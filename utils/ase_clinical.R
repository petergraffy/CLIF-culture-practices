# Clinical descriptive extensions; row-level objects are private intermediates only.
ase_early_phase <- function(hours) ifelse(hours<48,"First 48 ICU hours","After 48 ICU hours")

build_culture_antibiotic_timing <- function(con, events) {
  DBI::dbWriteTable(con,"clinical_culture_events",select(events,culture_event_id,hospitalization_id,collect_dttm),overwrite=TRUE)
  x<-DBI::dbGetQuery(con,"WITH doses AS (
    SELECT a.*, lag(admin_dttm) OVER(PARTITION BY a.hospitalization_id,med_category ORDER BY admin_dttm) previous_dose
    FROM antibiotics a JOIN hospitalizations h USING(hospitalization_id)
    WHERE is_iv_im=1 AND CAST(admin_dttm AS DATE)>=CAST(h.admission_dttm AS DATE)-2 AND admin_dttm<h.discharge_dttm
  ), starts AS (
    SELECT *, CASE WHEN previous_dose IS NULL OR date_diff('day',CAST(previous_dose AS DATE),CAST(admin_dttm AS DATE))>2 THEN 1 ELSE 0 END is_new FROM doses
  ), courses AS (
    SELECT *,sum(is_new) OVER(PARTITION BY hospitalization_id,med_category ORDER BY admin_dttm ROWS UNBOUNDED PRECEDING) course_id FROM starts
  ), course_doses AS (
    SELECT *,min(admin_dttm) OVER(PARTITION BY hospitalization_id,med_category,course_id) course_start FROM courses
  ), firsts AS (SELECT hospitalization_id,min(admin_dttm) first_dose FROM doses GROUP BY hospitalization_id)
  SELECT e.culture_event_id,(epoch(e.collect_dttm)-epoch(f.first_dose))/3600.0 hours_from_first_dose,
    f.first_dose IS NOT NULL has_recorded_parenteral_dose,
    EXISTS (SELECT 1 FROM course_doses d WHERE d.hospitalization_id=e.hospitalization_id AND d.admin_dttm<=e.collect_dttm AND d.admin_dttm>e.collect_dttm-INTERVAL '24 hours') dose_in_prior_24h,
    EXISTS (SELECT 1 FROM course_doses d WHERE d.hospitalization_id=e.hospitalization_id AND d.admin_dttm<=e.collect_dttm AND d.admin_dttm>e.collect_dttm-INTERVAL '24 hours' AND d.course_start>=e.collect_dttm-INTERVAL '48 hours') recent_new_course
  FROM clinical_culture_events e LEFT JOIN firsts f USING(hospitalization_id)")
  private<-left_join(events,x,by="culture_event_id") %>% mutate(antibiotic_relation=case_when(!has_recorded_parenteral_dose~"No recorded qualifying parenteral dose",hours_from_first_dose<0~"Before first dose",hours_from_first_dose==0~"Same recorded timestamp",TRUE~"After first dose"),
    recent_exposure=case_when(!dose_in_prior_24h~"No recorded dose in prior 24 hours",recent_new_course~"New drug course started within prior 48 hours",TRUE~"Ongoing drug course started more than 48 hours earlier"))
  summarize<-function(keys)private %>% group_by(across(all_of(keys))) %>% summarise(n_culture_events=n(),n_positive_culture_events=sum(any_positive_culture),
    n_hospitalizations=n_distinct(hospitalization_id),n_with_recorded_dose=sum(has_recorded_parenteral_dose),
    median_hours_from_first_dose=if(any(is.finite(hours_from_first_dose)))median(hours_from_first_dose,na.rm=TRUE) else NA_real_,.groups="drop") %>%
    group_by(ase_group,fluid_category) %>% mutate(percent_culture_events=100*n_culture_events/sum(n_culture_events),positivity=n_positive_culture_events/n_culture_events) %>% ungroup()
  private$fluid_category<-coalesce(private$fluid_category,"missing")
  list(private=private,relative=summarize(c("ase_group","fluid_category","antibiotic_relation")),exposure=summarize(c("ase_group","fluid_category","recent_exposure")))
}

build_icu_phase_analysis <- function(z, months) {
  phase_stays<-bind_rows(lapply(c("First 48 ICU hours","After 48 ICU hours"),function(phase) {
    s<-z$stays
    if(phase=="First 48 ICU hours")s$icu_out_dttm_clipped<-pmin(s$icu_out_dttm_clipped,s$icu_in_dttm+48*3600)
    else s$icu_in_dttm_clipped<-pmax(s$icu_in_dttm_clipped,s$icu_in_dttm+48*3600)
    s %>% filter(icu_out_dttm_clipped>icu_in_dttm_clipped) %>% mutate(icu_phase=phase,icu_los_days=as.numeric(difftime(icu_out_dttm_clipped,icu_in_dttm_clipped,units="days")))
  }))
  denominators<-phase_stays %>% group_by(ase_group,icu_phase) %>% summarise(n_icu_days=sum(icu_los_days),n_icu_admissions_exposed=n_distinct(icu_admission_id),.groups="drop")
  denominators<-tidyr::crossing(ase_group=ase_group_levels,icu_phase=c("First 48 ICU hours","After 48 ICU hours")) %>% left_join(denominators,by=c("ase_group","icu_phase")) %>% mutate(across(c(n_icu_days,n_icu_admissions_exposed),~coalesce(.x,0)))
  phase_events<-z$events %>% mutate(icu_phase=ase_early_phase(as.numeric(difftime(collect_dttm,icu_in_dttm,units="hours"))),fluid_category=coalesce(fluid_category,"missing"))
  phase_rows<-z$rows %>% mutate(icu_phase=ase_early_phase(as.numeric(difftime(collect_dttm,icu_in_dttm,units="hours"))))
  cultures<-phase_events %>% group_by(ase_group,icu_phase,fluid_category) %>% summarise(n_culture_events=n(),n_positive_culture_events=sum(any_positive_culture),.groups="drop")
  cultures<-tidyr::crossing(denominators,fluid_category=sort(unique(phase_events$fluid_category))) %>% left_join(cultures,by=c("ase_group","icu_phase","fluid_category")) %>%
    mutate(across(c(n_culture_events,n_positive_culture_events),~coalesce(.x,0L)),cultures_per_100_icu_days=if_else(n_icu_days>0,100*n_culture_events/n_icu_days,NA_real_),positivity=if_else(n_culture_events>0,n_positive_culture_events/n_culture_events,NA_real_))
  detections<-phase_rows %>% filter(positive_culture,!is.na(organism_category),!organism_category %in% c("na","unknown","missing")) %>% distinct(ase_group,icu_phase,icu_admission_id,collect_dttm,culture_event_id,organism_category)
  first<-detections %>% arrange(collect_dttm) %>% group_by(icu_admission_id,organism_category) %>% slice_head(n=1) %>% ungroup() %>% count(ase_group,icu_phase,organism_category,name="n_first_detection_icu_admissions")
  organisms<-detections %>% count(ase_group,icu_phase,organism_category,name="n_detection_events")
  organisms<-tidyr::crossing(denominators,organism_category=sort(unique(detections$organism_category))) %>% left_join(organisms,by=c("ase_group","icu_phase","organism_category")) %>% left_join(first,by=c("ase_group","icu_phase","organism_category")) %>%
    mutate(across(c(n_detection_events,n_first_detection_icu_admissions),~coalesce(.x,0L)),detections_per_100_icu_days=if_else(n_icu_days>0,100*n_detection_events/n_icu_days,NA_real_))
  list(stays=phase_stays,rows=phase_rows,denominators=denominators,cultures=cultures,organisms=organisms)
}

# Collapse simultaneous same-specimen collections; these are not validated blood-culture sets.
build_repeat_culture_yield <- function(z) {
  rows<-z$rows %>% mutate(fluid_category=coalesce(fluid_category,"missing"))
  collections<-rows %>% group_by(ase_group,hospitalization_id,icu_admission_id,fluid_category,collect_dttm) %>%
    summarise(result_status=collapse_result_status(result_status),organisms=list(sort(unique(organism_category[positive_culture & !is.na(organism_category) & !organism_category %in% c("na","unknown","missing")]))),
      result_available=if(n()>0 && all(!is.na(result_dttm)))max(result_dttm) else as.POSIXct(NA,tz="UTC"),.groups="drop") %>% arrange(icu_admission_id,fluid_category,collect_dttm)
  repeats<-list();n_index<-0L
  for(x in group_split(group_by(collections,icu_admission_id,fluid_category))) {
    if(!nrow(x))next
    index<-1L;n_index<-n_index+1L;seen<-x$organisms[[index]];unmapped_prior<-x$result_status[index]=="Positive" && !length(seen)
    for(i in seq_len(nrow(x))[-1L]) {
      hours<-as.numeric(difftime(x$collect_dttm[i],x$collect_dttm[index],units="hours"))
      if(hours>72){index<-i;n_index<-n_index+1L;seen<-x$organisms[[index]];unmapped_prior<-x$result_status[index]=="Positive" && !length(seen);next}
      novel<-setdiff(x$organisms[[i]],x$organisms[[index]])
      yield<-if(x$result_status[i]!="Positive")x$result_status[i] else if(!length(x$organisms[[i]]))"Positive organism unmapped" else if(x$result_status[index]=="Positive" && !length(x$organisms[[index]]))"Uncertain novelty: index organism unmapped" else if(length(novel))"New organism relative to index" else "Same organism only"
      new_to_episode<-setdiff(x$organisms[[i]],seen)
      incremental_yield<-if(x$result_status[i]!="Positive")x$result_status[i] else if(!length(x$organisms[[i]]))"Positive organism unmapped" else if(!length(new_to_episode))"Previously detected organism only" else if(unmapped_prior)"Uncertain novelty: prior positive organism unmapped" else "New organism not seen in prior episode collections"
      seen<-union(seen,x$organisms[[i]])
      unmapped_prior<-unmapped_prior || (x$result_status[i]=="Positive" && !length(x$organisms[[i]]))
      repeats[[length(repeats)+1L]]<-tibble(ase_group=x$ase_group[i],hospitalization_id=x$hospitalization_id[i],icu_admission_id=x$icu_admission_id[i],fluid_category=x$fluid_category[i],hours_since_index=hours,index_result=x$result_status[index],repeat_yield=yield,incremental_yield=incremental_yield,index_results_recorded_by_repeat=!is.na(x$result_available[index]) && x$result_available[index]<=x$collect_dttm[i])
    }
  }
  private<-if(length(repeats))bind_rows(repeats) else tibble(ase_group=character(),hospitalization_id=character(),icu_admission_id=integer(),fluid_category=character(),hours_since_index=numeric(),index_result=character(),repeat_yield=character(),incremental_yield=character(),index_results_recorded_by_repeat=logical())
  summary<-bind_rows(lapply(c(24,48,72),function(w)filter(private,hours_since_index<=w) %>% count(ase_group,fluid_category,index_result,repeat_yield,incremental_yield,index_results_recorded_by_repeat,name="n_repeat_collections") %>% mutate(window_hours=w))) %>%
    group_by(ase_group,fluid_category,index_result,window_hours) %>% mutate(n_repeat_collections_in_stratum=sum(n_repeat_collections),percent_repeat_collections=100*n_repeat_collections/n_repeat_collections_in_stratum) %>% ungroup()
  qc<-tibble(n_collection_episodes=n_index,n_distinct_collection_timestamps=nrow(collections),n_repeat_collections=nrow(private),max_index_window_hours=72,
    definition="same ICU admission and specimen; fixed 72-hour episode anchored at index collection; nested 24/48/72-hour windows")
  list(private=private,summary=summary,qc=qc)
}

build_ase_event_timing <- function(episodes, classification, full_stays, start, end, study_stays=NULL) {
  fields<-c("vasopressor_dttm","imv_dttm","aki_dttm","hyperbilirubinemia_dttm","thrombocytopenia_dttm")
  e<-filter(episodes,ase_without_lactate)
  for(f in c(fields,"blood_culture_dttm"))e[[f]]<-safe_ts(e[[f]])
  # Earliest qualifying dysfunction among all qualifying blood-culture anchors; not clinical recognition.
  times<-do.call(pmin,c(lapply(e[fields],as.numeric),list(na.rm=TRUE)))
  times[!is.finite(times)]<-NA_real_
  e$organ_dysfunction_proxy<-as.POSIXct(times,origin="1970-01-01",tz="UTC")
  primary<-e %>% arrange(organ_dysfunction_proxy,blood_culture_dttm,bc_id) %>% group_by(hospitalization_id) %>% slice_head(n=1) %>% ungroup() %>% transmute(hospitalization_id,timing_anchor="Earliest qualifying organ dysfunction",event_time=organ_dysfunction_proxy)
  blood<-e %>% arrange(blood_culture_dttm,bc_id) %>% group_by(hospitalization_id) %>% slice_head(n=1) %>% ungroup() %>% transmute(hospitalization_id,timing_anchor="First qualifying blood-culture anchor",event_time=blood_culture_dttm)
  full_stays<-full_stays %>% arrange(hospitalization_id,icu_in_dttm) %>% group_by(hospitalization_id) %>% mutate(icu_stay_sequence=row_number()) %>% ungroup()
  limits<-full_stays %>% group_by(hospitalization_id) %>% summarise(first_icu_in=min(icu_in_dttm),last_icu_out=max(icu_out_dttm),.groups="drop")
  x<-bind_rows(primary,blood) %>% left_join(limits,by="hospitalization_id")
  matched<-x %>% inner_join(select(full_stays,hospitalization_id,icu_stay_sequence,icu_in_dttm,icu_out_dttm),by="hospitalization_id",relationship="many-to-many") %>% filter(event_time>=icu_in_dttm,event_time<icu_out_dttm) %>% select(hospitalization_id,timing_anchor,icu_stay_sequence,icu_in_dttm)
  if(anyDuplicated(select(matched,hospitalization_id,timing_anchor)))stop("ASE event matches multiple ICU stays")
  x<-left_join(x,matched,by=c("hospitalization_id","timing_anchor")) %>% mutate(hours_from_first_icu_entry=as.numeric(difftime(event_time,first_icu_in,units="hours")),hours_into_matched_icu_stay=as.numeric(difftime(event_time,icu_in_dttm,units="hours")),
    within_study_window=event_time>=start & event_time<end,
    event_location=case_when(is.na(event_time)|is.na(first_icu_in)~"Unclassifiable timing",event_time<first_icu_in~"Before first ICU admission",!is.na(icu_in_dttm) & hours_into_matched_icu_stay<48~"During first 48 hours of an ICU stay",!is.na(icu_in_dttm) & hours_into_matched_icu_stay<=168~"During ICU stay: 48 hours to 7 days",!is.na(icu_in_dttm)~"During ICU stay: after 7 days",event_time>=last_icu_out~"After last ICU exit",TRUE~"Between ICU stays"))
  expected<-sum(classification$ase_group=="ASE")
  if(nrow(primary)!=expected || nrow(blood)!=expected)stop("ASE first-event timing does not cover all ASE hospitalizations")
  summary<-x %>% group_by(timing_anchor,event_location,within_study_window) %>% summarise(n_hospitalizations=n(),median_hours_from_first_icu_entry=if(any(is.finite(hours_from_first_icu_entry)))median(hours_from_first_icu_entry,na.rm=TRUE) else NA_real_,p25_hours_from_first_icu_entry=if(any(is.finite(hours_from_first_icu_entry)))unname(quantile(hours_from_first_icu_entry,.25,na.rm=TRUE)) else NA_real_,p75_hours_from_first_icu_entry=if(any(is.finite(hours_from_first_icu_entry)))unname(quantile(hours_from_first_icu_entry,.75,na.rm=TRUE)) else NA_real_,.groups="drop") %>% group_by(timing_anchor) %>% mutate(n_ase_hospitalizations=sum(n_hospitalizations),percent_ase_hospitalizations=100*n_hospitalizations/n_ase_hospitalizations) %>% ungroup()
  hours<-x %>% filter(!is.na(icu_in_dttm)) %>% group_by(timing_anchor) %>% summarise(n_hospitalizations=n(),n_during_first_icu_stay=sum(icu_stay_sequence==1),n_during_later_icu_stay=sum(icu_stay_sequence>1),median_hours_into_icu_stay=median(hours_into_matched_icu_stay),p25_hours=unname(quantile(hours_into_matched_icu_stay,.25)),p75_hours=unname(quantile(hours_into_matched_icu_stay,.75)),.groups="drop")
  day<-x %>% filter(!is.na(icu_in_dttm)) %>% mutate(icu_day=pmin(floor(hours_into_matched_icu_stay/24)+1,15)) %>% count(timing_anchor,icu_day,name="n_hospitalizations") %>% mutate(day_label=if_else(icu_day==15,"Day 15 or later",paste("Day",icu_day)),n_ase_hospitalizations=expected,percent_all_ase_hospitalizations=100*n_hospitalizations/expected)
  if(expected>0)day<-tidyr::crossing(timing_anchor=c("Earliest qualifying organ dysfunction","First qualifying blood-culture anchor"),icu_day=1:15) %>% left_join(select(day,timing_anchor,icu_day,n_hospitalizations),by=c("timing_anchor","icu_day")) %>% mutate(n_hospitalizations=coalesce(n_hospitalizations,0L),day_label=if_else(icu_day==15,"Day 15 or later",paste("Day",icu_day)),n_ase_hospitalizations=expected,percent_all_ase_hospitalizations=100*n_hospitalizations/expected)
  qc<-tibble(timing_anchor=c("Earliest qualifying organ dysfunction","First qualifying blood-culture anchor"),n_ase_hospitalizations=expected,n_first_anchors=c(nrow(primary),nrow(blood)),n_unclassifiable=c(sum(is.na(primary$event_time)),sum(is.na(blood$event_time))))
  admission_summary<-tibble(timing_anchor=character(),timing_status=character(),n_icu_stays=integer(),n_ase_icu_stays=integer(),percent_ase_icu_stays=numeric())
  if(!is.null(study_stays)) {
    admissions<-bind_rows(primary,blood) %>% inner_join(filter(study_stays,ase_group=="ASE"),by="hospitalization_id",relationship="many-to-many") %>%
      mutate(hours_from_icu_entry=as.numeric(difftime(event_time,icu_in_dttm,units="hours")),timing_status=case_when(is.na(event_time)~"Unclassifiable timing",event_time<=icu_in_dttm~"First ASE anchor at or before ICU entry",event_time>=icu_out_dttm~"First ASE anchor after this ICU stay",hours_from_icu_entry<48~"First ASE anchor within first 48 ICU hours",TRUE~"First ASE anchor later in this ICU stay"))
    admission_summary<-admissions %>% count(timing_anchor,timing_status,name="n_icu_stays") %>% group_by(timing_anchor) %>% mutate(n_ase_icu_stays=sum(n_icu_stays),percent_ase_icu_stays=100*n_icu_stays/n_ase_icu_stays) %>% ungroup()
  }
  list(private=x,summary=summary,hours=hours,day=day,admissions=admission_summary,qc=qc)
}
