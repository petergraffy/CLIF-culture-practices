# Additional retrospective hospitalization subgroups; whole-cohort outputs are unchanged.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(lubridate);library(readr);library(ggplot2)})
source("utils/clif_io.R");source("utils/preflight.R");source("utils/ase.R");source("utils/ase_characteristics.R");source("utils/ase_cumulative.R");source("utils/ase_culture_density.R");source("utils/ase_clinical.R");source("utils/trends.R");source("utils/susceptibility.R")
out<-project_output_dir("ase");stamp<-format(Sys.time(),"%Y%m%d_%H%M%S")
write_ase<-function(x,name)write_csv(mutate(x,site_name=clif_site_name),file.path(out,paste0(name,"_",clif_site_name,"_",stamp,".csv")))
availability<-function(status,detail=NA_character_)write_ase(tibble(analysis_status=status,detail=detail),"ase_availability")

ase_plot <- function(p) {
  colors<-setNames(c("#0072B2","#E69F00","#CC79A7"),ase_group_levels)
  labels<-function(x)stringr::str_wrap(x,24)
  p+scale_color_manual(values=colors,limits=ase_group_levels,labels=labels,name=NULL)+
    scale_fill_manual(values=colors,limits=ase_group_levels,labels=labels,name=NULL)
}

ase_fit_screens <- function(data,keys,outcomes) {
  models<-curves<-list();screens<-distinct(data,across(all_of(keys)))
  for(i in seq_len(nrow(screens))) {
    x<-semi_join(data,screens[i,,drop=FALSE],by=keys)
    for(j in seq_len(nrow(outcomes))) {
      spec<-outcomes[j,];selected<-x
      if(spec$proportion)selected<-filter(x,fraction_model_eligible)
      else if("rate_model_eligible" %in% names(x))selected<-filter(x,rate_model_eligible)
      fit<-fit_temporal_model(selected,spec$denominator,spec$count,spec$proportion)
      models[[length(models)+1L]]<-bind_cols(screens[i,,drop=FALSE],fit$summary) %>% mutate(outcome=spec$outcome,denominator=spec$denominator,n_usable_months=nrow(selected),first_usable_month=if(nrow(selected))min(selected$calendar_month) else as.POSIXct(NA),last_usable_month=if(nrow(selected))max(selected$calendar_month) else as.POSIXct(NA),effect_scale=if(spec$proportion)"annualized endpoint odds ratio" else "annualized endpoint rate ratio")
      if(nrow(fit$predictions))curves[[length(curves)+1L]]<-bind_cols(screens[rep(i,nrow(fit$predictions)),,drop=FALSE],fit$predictions) %>% mutate(outcome=spec$outcome,denominator=spec$denominator)
    }
  }
  models<-bind_rows(models)
  if(nrow(models))models<-models %>% group_by(across(all_of(c("ase_group","outcome","denominator",intersect("specimen_stratum",names(models)))))) %>% mutate(fdr_p_value=p.adjust(p_value,"BH"),direction=classify_temporal_direction(model_status,fdr_p_value,annual_percent_change,residual_dependence_flag)) %>% ungroup()
  list(models=models,curves=bind_rows(curves))
}

run_ase_analysis <- function() {
  check<-ase_source_check()
  if(check$status!="available") {availability(check$status,check$detail);message("ASE analysis: ",check$status," (",check$detail,")");return(invisible(NULL))}
  preflight_packages();if(!requireNamespace("duckdb",quietly=TRUE))stop("ASE requires DuckDB. Restore the updated renv.lock before running.")
  availability("analysis_started")
  start<-safe_ts(study_settings$study_start_date);end<-safe_ts(study_settings$study_end_date)+days(1)
  data<-read_culture_data(start,end)
  months<-seq(floor_date(start,"month"),floor_date(end-seconds(1),"month"),by="month")
  private<-project_intermediate_dir("ase")
  con<-ase_connect(file.path(private,"duckdb_tmp"));on.exit(DBI::dbDisconnect(con,shutdown=TRUE),add=TRUE)
  ase_register_sources(con,check$paths,unique(data$icu_admissions$hospitalization_id));ase_standardize_sources(con)
  # Validate observation labels before defining non-ASE from their absence.
  observed<-DBI::dbGetQuery(con,"SELECT 'antimicrobials' component,count(*) n FROM src_medication_admin_intermittent WHERE med_group='cms_sepsis_qualifying_antibiotics' AND mar_action_group='administered' UNION ALL SELECT 'vasopressors',count(*) FROM src_medication_admin_continuous WHERE med_group='vasoactives' AND mar_action_group='administered' UNION ALL SELECT 'labs',count(*) FROM src_labs WHERE lab_category IN ('creatinine','bilirubin_total','platelet_count') UNION ALL SELECT 'respiratory',count(*) FROM src_respiratory_support WHERE device_category IS NOT NULL")
  write_ase(observed,"ase_source_component_activity")
  if(any(observed$n==0)) {availability("skipped_unverified_component_coverage",paste(observed$component[observed$n==0],collapse="; "));return(invisible(NULL))}
  result<-compute_ase_hospitalizations(con,project_path("utils","ase","sql"))
  classification<-result$classification
  if(any(classification$ase_without_lactate & !classification$presumed_infection))stop("ASE requires presumed infection")
  if(anyNA(classification$ase_group)||any(!classification$ase_group %in% ase_group_levels))stop("Invalid subgroup classification")
  # This output contains identifiers and stays private.
  write_csv(classification,file.path(private,"ase_hospitalization_classification.csv"))
  if(nrow(result$episodes))write_csv(result$episodes,file.path(private,"ase_blood_culture_criteria.csv"))
  h<-data$hospitalization %>% semi_join(data$icu_admissions,by="hospitalization_id")
  qc<-tibble(metric=c("icu_hospitalizations","classified_adult_hospitalizations","excluded_or_unclassifiable_hospitalizations","ase_hospitalizations","presumed_infection_without_ase_hospitalizations","no_presumed_infection_hospitalizations","all_presumed_infection_hospitalizations","additional_ase_hospitalizations_with_lactate"),n=c(n_distinct(h$hospitalization_id),nrow(classification),n_distinct(h$hospitalization_id)-nrow(classification),sum(classification$ase_group=="ASE"),sum(classification$ase_group=="Presumed infection without ASE"),sum(classification$ase_group=="No presumed infection"),sum(classification$presumed_infection),sum(classification$ase_with_lactate & !classification$ase_without_lactate)))
  write_ase(qc,"ase_classification_qc")
  excluded<-h %>% anti_join(classification,by="hospitalization_id") %>% mutate(reason=case_when(is.na(age_at_admission)~"missing_age",age_at_admission<18~"not_adult",TRUE~"invalid_hospitalization_boundaries")) %>% count(reason,name="n_hospitalizations")
  write_ase(excluded,"ase_classification_exclusions")
  if(!nrow(classification)){availability("no_classifiable_adult_hospitalizations");return(invisible(NULL))}
  clinical_characteristics<-DBI::dbGetQuery(con,"SELECT h.hospitalization_id,
    CASE WHEN EXISTS (SELECT 1 FROM src_respiratory_support r WHERE r.hospitalization_id=h.hospitalization_id
      AND r.device_category='imv' AND r.recorded_dttm>=h.admission_dttm AND r.recorded_dttm<h.discharge_dttm) THEN 'yes' ELSE 'no' END recorded_imv,
    CASE WHEN EXISTS (SELECT 1 FROM med_continuous m WHERE m.hospitalization_id=h.hospitalization_id AND m.med_dose>0
      AND m.admin_dttm>=h.admission_dttm AND m.admin_dttm<h.discharge_dttm
      AND NOT EXISTS (SELECT 1 FROM src_adt a WHERE a.hospitalization_id=m.hospitalization_id
        AND a.location_category='procedural' AND m.admin_dttm>=a.in_dttm AND m.admin_dttm<a.out_dttm)) THEN 'yes' ELSE 'no' END recorded_vasopressor
    FROM src_hospitalization h")
  characteristics<-build_ase_characteristics(data,classification,read_tbl("patient"),clinical_characteristics)
  write_ase(characteristics$long,"ase_characteristics_long")
  write_ase(characteristics$wide,"ase_characteristics_table1")
  write_ase(characteristics$qc,"ase_characteristics_qc")
  z<-build_ase_group_aggregates(data,classification,months)
  eligible_stays<-semi_join(data$icu_admissions,classification,by="hospitalization_id")
  eligible_events<-semi_join(data$events,classification,by="hospitalization_id")
  reconciliation<-tibble(metric=c("hospitalizations","icu_stays","icu_days","culture_events","positive_culture_events"),
    expected=c(nrow(classification),nrow(eligible_stays),sum(eligible_stays$icu_los_days),nrow(eligible_events),sum(eligible_events$any_positive_culture)),
    subgroup_total=c(sum(z$cohort$n_hospitalizations),sum(z$cohort$n_icu_stays),sum(z$cohort$n_icu_days),sum(z$monthly$n_culture_events),sum(z$monthly$n_positive_culture_events))) %>%
    mutate(difference=subgroup_total-expected,reconciled=abs(difference)<1e-6)
  write_ase(reconciliation,"ase_cohort_reconciliation")
  if(any(!reconciliation$reconciled))stop("Three-group cohort denominators or culture counts failed reconciliation")
  write_ase(z$admission_monthly,"ase_monthly_admission_culture_proportions");write_ase(z$timing_bins,"ase_first_culture_timing_bins");write_ase(z$cohort,"ase_cohort_summary");write_ase(z$monthly,"ase_monthly_culture_rates");write_ase(z$specimen,"ase_monthly_specimen_culture_rates");write_ase(z$organisms,"ase_monthly_organism_detection_counts");write_ase(z$timing,"ase_first_culture_timing")
  write_ase(z$rows %>% filter(positive_culture) %>% group_by(ase_group,organism_group,organism_category) %>% summarise(n_detection_events=n_distinct(culture_event_id),n_hospitalizations=n_distinct(hospitalization_id),n_patients=n_distinct(patient_id),.groups="drop"),"ase_organism_distribution")
  antibiotics_timing<-build_culture_antibiotic_timing(con,z$events)
  write_ase(antibiotics_timing$relative,"ase_culture_antibiotic_timing")
  write_ase(antibiotics_timing$exposure,"ase_culture_recent_antibiotic_exposure")
  write_csv(antibiotics_timing$private,file.path(private,"culture_antibiotic_timing.csv"))
  phases<-build_icu_phase_analysis(z,months)
  write_ase(phases$denominators,"ase_early_later_icu_denominators")
  write_ase(phases$cultures,"ase_early_later_specimen_yield")
  write_ase(phases$organisms,"ase_early_later_organism_detection")
  repeats<-build_repeat_culture_yield(z)
  write_ase(repeats$summary,"ase_repeat_culture_yield")
  write_ase(repeats$qc,"ase_repeat_culture_qc")
  write_csv(repeats$private,file.path(private,"repeat_culture_comparisons.csv"))
  clinical_reconciliation<-tibble(metric=c("antibiotic_timing_culture_events","phase_culture_events","phase_icu_days","repeat_collection_accounting"),
    expected=c(nrow(z$events),nrow(z$events),sum(z$stays$icu_los_days),repeats$qc$n_distinct_collection_timestamps),
    observed=c(sum(antibiotics_timing$relative$n_culture_events),sum(phases$cultures$n_culture_events),sum(phases$denominators$n_icu_days),repeats$qc$n_collection_episodes+repeats$qc$n_repeat_collections)) %>% mutate(reconciled=abs(observed-expected)<1e-6)
  write_ase(clinical_reconciliation,"ase_clinical_reconciliation")
  if(any(!clinical_reconciliation$reconciled))stop("Clinical extensions failed cohort reconciliation")
  full_stays<-merge_icu_stays(DBI::dbGetQuery(con,"SELECT hospitalization_id,in_dttm,out_dttm,location_category FROM src_adt"),data$hospitalization)
  event_timing<-build_ase_event_timing(result$episodes,classification,full_stays,start,end,z$stays)
  write_ase(event_timing$qc,"ase_event_timing_qc")
  write_ase(event_timing$admissions,"ase_first_event_relative_to_icu_admissions")
  write_ase(event_timing$summary,"ase_event_location_timing")
  write_ase(event_timing$hours,"ase_event_hours_into_icu_stay")
  write_ase(event_timing$day,"ase_event_icu_day_distribution")
  write_csv(event_timing$private,file.path(private,"ase_first_event_timing.csv"))
  if(nrow(event_timing$summary)) {
  p<-ggplot(event_timing$summary,aes(stringr::str_wrap(event_location,25),percent_ase_hospitalizations,fill=within_study_window))+geom_col()+scale_x_discrete(limits=rev(stringr::str_wrap(c("Before first ICU admission","During first 48 hours of an ICU stay","During ICU stay: 48 hours to 7 days","During ICU stay: after 7 days","Between ICU stays","After last ICU exit","Unclassifiable timing"),25)))+facet_wrap(~timing_anchor,ncol=1)+coord_flip()+labs(x=NULL,y="ASE hospitalizations (%)",fill="Within study dates",caption=NULL)+theme_bw()
  ggsave(file.path(out,paste0("ase_event_location_timing_",clif_site_name,"_",stamp,".png")),p,width=12,height=9,dpi=200)
  }
  if(nrow(event_timing$day)) {
  p<-ggplot(event_timing$day,aes(icu_day,percent_all_ase_hospitalizations,color=timing_anchor))+geom_point()+geom_line()+labs(x="Day within matched ICU stay (15 = day 15 or later)",y="All ASE hospitalizations (%)",color=NULL,caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_event_icu_day_distribution_",clif_site_name,"_",stamp,".png")),p,width=12,height=7,dpi=200)
  }
  p<-ggplot(filter(antibiotics_timing$relative,fluid_category=="blood_buffy",n_culture_events>=20),aes(stringr::str_wrap(antibiotic_relation,20),100*positivity,color=ase_group))+geom_point(position=position_dodge(width=.5),size=2)+labs(x=NULL,y="Positive blood-culture events (%)",caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_blood_culture_antibiotic_timing_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=12,height=7,dpi=200)
  p<-ggplot(phases$cultures %>% group_by(ase_group,icu_phase) %>% summarise(n=sum(n_culture_events),days=first(n_icu_days),rate=if_else(days>0,100*n/days,NA_real_),.groups="drop"),aes(icu_phase,rate,fill=ase_group))+geom_col(position="dodge")+labs(x=NULL,y="Culture events per 100 phase-specific ICU-days",caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_early_later_culture_rates_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=12,height=7,dpi=200)
  top_phase_organisms<-phases$organisms %>% group_by(organism_category) %>% summarise(n=sum(n_detection_events),.groups="drop") %>% slice_max(n,n=10,with_ties=FALSE)
  if(nrow(top_phase_organisms)) {
    p<-ggplot(semi_join(phases$organisms,top_phase_organisms,by="organism_category"),aes(icu_phase,detections_per_100_icu_days,fill=ase_group))+geom_col(position="dodge")+facet_wrap(~organism_category,scales="free_y",ncol=3)+labs(x=NULL,y="Organism detections per 100 phase-specific ICU-days",caption=NULL)+theme_bw()+theme(legend.position="bottom",axis.text.x=element_text(angle=20,hjust=1))
    ggsave(file.path(out,paste0("ase_early_later_organism_rates_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=14,height=11,dpi=200)
  }
  repeat_plot<-repeats$summary %>% group_by(ase_group,window_hours,index_result) %>% summarise(n=sum(n_repeat_collections),n_new=sum(n_repeat_collections[incremental_yield=="New organism not seen in prior episode collections"]),percent_new=100*n_new/n,.groups="drop")
  repeat_plot<-filter(repeat_plot,n>=20)
  if(nrow(repeat_plot)) {
    p<-ggplot(repeat_plot,aes(factor(window_hours),percent_new,color=ase_group,group=ase_group))+geom_point()+geom_line()+facet_wrap(~index_result)+labs(x="Hours from index collection (nested windows)",y="Repeat collections detecting an organism new to the episode (%)",caption=NULL)+theme_bw()+theme(legend.position="bottom")
    ggsave(file.path(out,paste0("ase_repeat_culture_incremental_yield_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=12,height=8,dpi=200)
  }
  cumulative<-build_ase_cumulative(z,as.numeric(Sys.getenv("TIMING_MAX_ICU_DAY","14")),as.numeric(Sys.getenv("TIMING_MAX_ICU_HOUR","168")))
  write_ase(cumulative$first,"ase_cumulative_first_culture_by_icu_day")
  write_ase(cumulative$organisms,"ase_cumulative_organism_detection_by_icu_hour")
  write_ase(cumulative$cultures,"ase_cumulative_culture_events_by_specimen_icu_hour")
  write_ase(cumulative$observation,"ase_cumulative_observation_counts")
  write_ase(cumulative$qc,"ase_cumulative_qc")
  culture_density<-build_ase_culture_density(cumulative$cultures)
  write_ase(culture_density$curves,"ase_culture_timing_density")
  write_ase(culture_density$counts,"ase_culture_timing_counts")
  write_ase(culture_density$hourly_counts,"ase_culture_timing_histogram_1h")
  write_ase(culture_density$specimen_hourly_counts,"ase_specimen_histogram_1h")
  plot_ase_specimen_histograms(culture_density,out,clif_site_name,stamp)
  plot_ase_culture_density(culture_density,out,clif_site_name,stamp)
  cumulative_caption<-"Observed cumulative values; denominator includes uncultured new ICU admissions.\nRetrospective hospitalization groups; no censoring or competing-risk adjustment."
  p<-ggplot(cumulative$first,aes(icu_day,cumulative_percent,color=ase_group))+geom_step()+facet_wrap(~outcome)+
    labs(x="Days since ICU admission",y="ICU admissions with event (%)",caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_cumulative_first_culture_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=12,height=7,dpi=200)
  top_cumulative<-cumulative$organisms %>% distinct(ase_group,organism_category,n_total_detection_events) %>% group_by(organism_category) %>% summarise(n=sum(n_total_detection_events),.groups="drop") %>% slice_max(n,n=as.integer(Sys.getenv("TOP_N_TIMING_ORGANISMS","10")),with_ties=FALSE)
  if(nrow(top_cumulative)) {
    p<-ggplot(semi_join(cumulative$organisms,top_cumulative,by="organism_category"),aes(icu_hour,cumulative_percent,color=ase_group))+geom_step()+facet_wrap(~organism_category,scales="free_y",ncol=3)+
      labs(x="Hours since ICU admission",y="ICU admissions with organism detected (%)",caption=NULL)+theme_bw()+theme(legend.position="bottom")
    ggsave(file.path(out,paste0("ase_cumulative_organism_detection_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=14,height=11,dpi=200)
  }
  top_cumulative_specimens<-cumulative$cultures %>% distinct(ase_group,fluid_category,n_total_culture_events) %>% group_by(fluid_category) %>% summarise(n=sum(n_total_culture_events),.groups="drop") %>% slice_max(n,n=as.integer(Sys.getenv("TOP_N_CULTURE_TYPES","8")),with_ties=FALSE)
  if(nrow(top_cumulative_specimens)) {
    p<-ggplot(semi_join(cumulative$cultures,top_cumulative_specimens,by="fluid_category"),aes(icu_hour,cumulative_events_per_100_admissions,color=ase_group))+geom_step()+facet_wrap(~fluid_category,scales="free_y",ncol=2)+
      labs(x="Hours since ICU admission",y="Cumulative culture events per 100 ICU admissions",caption=NULL)+theme_bw()+theme(legend.position="bottom")
    ggsave(file.path(out,paste0("ase_cumulative_culture_events_by_specimen_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=13,height=10,dpi=200)
  }
  culture_outcomes<-tibble(outcome=c("culture_collection","culture_collection","positive_culture_detection","positive_culture_detection"),count=c("n_culture_events","n_culture_events","n_positive_culture_events","n_positive_culture_events"),denominator=rep(c("n_icu_days","n_icu_admissions"),2),proportion=FALSE)
  culture<-ase_fit_screens(z$monthly,"ase_group",culture_outcomes)
  write_ase(culture$models,"ase_culture_temporal_models");if(nrow(culture$curves))write_ase(culture$curves,"ase_culture_fitted_curves")
  detection<-ase_fit_screens(z$organisms,c("ase_group","organism_category"),tibble(outcome="organism_detection",count="n_detection_events",denominator=c("n_icu_days","n_icu_admissions"),proportion=FALSE))
  write_ase(detection$models,"ase_organism_temporal_models");if(nrow(detection$curves))write_ase(detection$curves,"ase_organism_fitted_curves")
  caption<-"Three mutually exclusive retrospective hospitalization groups; primary ASE excludes lactate.\nPoints: monthly observations. Curves: season-adjusted GAM means and pointwise 95% CIs."
  points<-z$monthly %>% select(ase_group,calendar_month,culture_events_per_100_icu_days,positive_events_per_100_icu_days) %>% pivot_longer(ends_with("icu_days"),names_to="outcome",values_to="observed") %>% mutate(outcome=if_else(outcome=="culture_events_per_100_icu_days","culture_collection","positive_culture_detection"))
  p<-ggplot(points,aes(calendar_month,observed,color=ase_group))+geom_point(alpha=.4)+facet_wrap(~outcome,scales="free_y",ncol=1)+labs(x=NULL,y="Events per 100 subgroup ICU-days",color=NULL,caption=NULL)+theme_bw()+theme(legend.position="bottom")
  if(nrow(culture$curves)){curves<-filter(culture$curves,denominator=="n_icu_days");p<-p+geom_ribbon(data=curves,aes(y=NULL,ymin=ci_low,ymax=ci_high,fill=ase_group),alpha=.12,color=NA)+geom_line(data=curves,aes(y=fitted_value))}
  ggsave(file.path(out,paste0("ase_culture_rates_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=12,height=8,dpi=200)
  p<-ggplot(z$monthly,aes(calendar_month,positivity,color=ase_group))+geom_point()+geom_line()+scale_y_continuous(labels=scales::percent)+labs(x=NULL,y="Positive fraction of subgroup culture events",color=NULL,caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_culture_positivity_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=12,height=6,dpi=200)
  top<-z$organisms %>% group_by(organism_category) %>% summarise(n=sum(n_detection_events),.groups="drop") %>% slice_max(n,n=12,with_ties=FALSE)
  if(nrow(top)) {
    points<-semi_join(z$organisms,top,by="organism_category")
    p<-ggplot(points,aes(calendar_month,detections_per_100_icu_days,color=ase_group))+geom_point(alpha=.35)+facet_wrap(~organism_category,scales="free_y",ncol=3)+labs(x=NULL,y="Detections per 100 subgroup ICU-days",color=NULL,caption=NULL)+theme_bw()+theme(legend.position="bottom")
    if(nrow(detection$curves)){curves<-detection$curves %>% filter(denominator=="n_icu_days") %>% semi_join(top,by="organism_category");p<-p+geom_ribbon(data=curves,aes(y=NULL,ymin=ci_low,ymax=ci_high,fill=ase_group),alpha=.12,color=NA)+geom_line(data=curves,aes(y=fitted_value))}
    ggsave(file.path(out,paste0("ase_organism_detection_rates_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=14,height=11,dpi=200)
  }
  p<-ggplot(z$timing,aes(ase_group,percent_with_icu_culture,fill=ase_group))+geom_col()+scale_x_discrete(limits=ase_group_levels,labels=function(x)stringr::str_wrap(x,20))+labs(x=NULL,y="ICU admissions with a culture (%)",caption=NULL)+theme_bw()+theme(legend.position="bottom")+theme(legend.position="none")
  ggsave(file.path(out,paste0("ase_admissions_cultured_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=7,height=5,dpi=200)
  top_specimens<-z$specimen %>% group_by(fluid_category) %>% summarise(n=sum(n_culture_events),.groups="drop") %>% slice_max(n,n=6,with_ties=FALSE)
  specimens<-semi_join(z$specimen,top_specimens,by="fluid_category")
  for(metric in c("cultures_per_100_icu_days","positivity")) {
    p<-ggplot(specimens,aes(calendar_month,.data[[metric]],color=ase_group))+geom_point(alpha=.5)+geom_line()+facet_wrap(~fluid_category,scales="free_y",ncol=2)+labs(x=NULL,y=if(metric=="positivity")"Positive fraction" else "Cultures per 100 subgroup ICU-days",color=NULL,caption=NULL)+theme_bw()+theme(legend.position="bottom")
    if(metric=="positivity")p<-p+scale_y_continuous(labels=scales::percent)
    ggsave(file.path(out,paste0("ase_specimen_",metric,"_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=13,height=10,dpi=200)
  }
  ast_path<-find_table_path("microbiology_susceptibility",required=FALSE)
  if(is.na(ast_path))write_ase(tibble(analysis_status="skipped_table_unavailable"),"ase_susceptibility_availability") else {
    ast<-read_tbl("microbiology_susceptibility");monthly_ast<-qc_ast<-list()
    coverage<-as.numeric(Sys.getenv("AST_MIN_TESTING_FRACTION","0.5"));linkage<-as.numeric(Sys.getenv("AST_MIN_LINKAGE_FRACTION","0.9"))
    if(any(!is.finite(c(coverage,linkage)))||any(c(coverage,linkage)<0|c(coverage,linkage)>1))stop("Invalid AST coverage/linkage thresholds")
    phase_ast<-phase_ast_qc<-list()
    for(g in ase_group_levels) for(phase in c("First 48 ICU hours","After 48 ICU hours")) {
      s<-filter(phases$stays,ase_group==g,icu_phase==phase)
      if(!nrow(s))next
      a_phase<-build_susceptibility_analysis(filter(phases$rows,ase_group==g,icu_phase==phase),ast,monthly_icu_denominators(s,months),project_path("config","mcide"),coverage,linkage)
      phase_ast_qc[[length(phase_ast_qc)+1L]]<-mutate(a_phase$qc,ase_group=g,icu_phase=phase)
      if(nrow(a_phase$monthly))phase_ast[[length(phase_ast)+1L]]<-mutate(a_phase$monthly,ase_group=g,icu_phase=phase)
    }
    write_ase(bind_rows(phase_ast_qc),"ase_early_later_susceptibility_qc")
    if(length(phase_ast))write_ase(bind_rows(phase_ast),"ase_early_later_monthly_susceptibility")
    for(g in ase_group_levels) {
      group_stays<-filter(z$stays,ase_group==g);group_rows<-filter(z$rows,ase_group==g)
      if(!nrow(group_stays))next
      a<-build_susceptibility_analysis(group_rows,ast,monthly_icu_denominators(group_stays,months),project_path("config","mcide"),coverage,linkage)
      qc_ast[[g]]<-mutate(a$qc,ase_group=g)
      if(nrow(a$monthly))monthly_ast[[g]]<-mutate(a$monthly,ase_group=g)
    }
    write_ase(bind_rows(qc_ast),"ase_susceptibility_qc")
    monthly_ast<-bind_rows(monthly_ast)
    if(nrow(monthly_ast)) {
      write_ase(monthly_ast,"ase_monthly_organism_antimicrobial_susceptibility")
      a<-ase_fit_screens(monthly_ast,c("ase_group","organism_category","antimicrobial_category","specimen_stratum"),tibble(outcome=c("susceptible_detection_rate","non_susceptible_detection_rate","non_susceptible_fraction"),count=c("n_susceptible","n_non_susceptible","n_non_susceptible"),denominator=c("n_icu_days","n_icu_days","n_interpretable"),proportion=c(FALSE,FALSE,TRUE)))
      write_ase(a$models,"ase_susceptibility_temporal_models");if(nrow(a$curves))write_ase(a$curves,"ase_susceptibility_fitted_curves")
      top_pairs<-monthly_ast %>% filter(specimen_stratum=="Overall") %>% group_by(organism_category,antimicrobial_category) %>% summarise(n=sum(n_interpretable),.groups="drop") %>% slice_max(n,n=6,with_ties=FALSE)
      points<-monthly_ast %>% filter(specimen_stratum=="Overall") %>% semi_join(top_pairs,by=c("organism_category","antimicrobial_category")) %>% mutate(pair=paste(organism_category,antimicrobial_category,sep=" / "))
      for(outcome in c("susceptible_detection_rate","non_susceptible_detection_rate","non_susceptible_fraction")) {
        fraction<-outcome=="non_susceptible_fraction"
        observed<-points %>% mutate(observed=if(fraction)non_susceptible_fraction else if(outcome=="susceptible_detection_rate")susceptible_per_100_icu_days else non_susceptible_per_100_icu_days,eligible=if(fraction)fraction_model_eligible else rate_model_eligible)
        p<-ggplot(observed,aes(calendar_month,observed,color=ase_group))+geom_point(aes(shape=eligible),alpha=.5,na.rm=TRUE)+facet_wrap(~pair,scales=if(fraction)"fixed" else "free_y",ncol=2)+labs(x=NULL,y=if(fraction)"Non-susceptible fraction among interpretable tests" else "Detections per 100 subgroup ICU-days",color=NULL,shape="Model eligible",caption=NULL)+theme_bw()+theme(legend.position="bottom")
        if(fraction)p<-p+scale_y_continuous(labels=scales::percent,limits=c(0,1))
        if(nrow(a$curves)) {
          fitted<-a$curves %>% filter(specimen_stratum=="Overall",.data$outcome==.env$outcome) %>% semi_join(top_pairs,by=c("organism_category","antimicrobial_category")) %>% mutate(pair=paste(organism_category,antimicrobial_category,sep=" / "))
          p<-p+geom_ribbon(data=fitted,aes(y=NULL,ymin=ci_low,ymax=ci_high,fill=ase_group),color=NA,alpha=.12)+geom_line(data=fitted,aes(y=fitted_value))
        }
        ggsave(file.path(out,paste0("ase_",outcome,"_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=13,height=10,dpi=200)
      }
      coverage_points<-points %>% pivot_longer(c(linkage_fraction,testing_fraction),names_to="coverage_metric",values_to="fraction")
      p<-ggplot(coverage_points,aes(calendar_month,fraction,color=ase_group,shape=coverage_metric))+geom_point(na.rm=TRUE)+facet_wrap(~pair,ncol=2)+scale_y_continuous(limits=c(0,1),labels=scales::percent)+geom_hline(yintercept=c(coverage,linkage),linetype=2,color="grey60")+labs(x=NULL,y="Linkage / interpretable testing coverage",color=NULL,shape=NULL)+theme_bw()+theme(legend.position="bottom")
      ggsave(file.path(out,paste0("ase_susceptibility_coverage_",clif_site_name,"_",stamp,".png")),ase_plot(p),width=13,height=10,dpi=200)
      write_ase(tibble(analysis_status="completed",n_estimated_models=sum(a$models$model_status=="estimated")),"ase_susceptibility_availability")
    } else write_ase(tibble(analysis_status="no_linked_interpretable_analysis"),"ase_susceptibility_availability")
  }
  upstream<-jsonlite::fromJSON(project_path("utils","ase","upstream.json"))
  jsonlite::write_json(list(classification="mutually exclusive full-hospitalization groups: ASE, presumed infection without ASE, no presumed infection",groups=as.list(ase_group_levels),classification_priority="ASE first; then presumed infection without ASE; otherwise no presumed infection",primary_include_lactate=FALSE,adult_minimum_age=18,calendar_timezone="UTC",denominators="subgroup ICU-days and new ICU admissions; includes uncultured stays",upstream=upstream,clinical_extensions=list(ase_timing="earliest qualifying organ dysfunction proxy and first qualifying blood culture; one per hospitalization per anchor",icu_phase_boundary_hours=48,repeat_index_episode_hours=72,repeat_windows_hours=list(24,48,72),antibiotics="recorded administered positive-dose CMS-sepsis-qualifying IV/IM drugs; first hospitalization dose and prior-24-hour exposure"),cumulative=list(estimator="observed cumulative proportion",time_origin="new ICU admission",carry_in_stays="excluded",denominator="all subgroup new ICU admissions, including uncultured",first_culture_day_horizon=as.numeric(Sys.getenv("TIMING_MAX_ICU_DAY","14")),hour_horizon=as.numeric(Sys.getenv("TIMING_MAX_ICU_HOUR","168")),competing_risk_adjusted=FALSE),characteristics=list(analysis_unit="hospitalization",categorical_denominator="all subgroup hospitalizations including missing",continuous_summary="median and quartiles among observed values",hospitalization_outcomes="full hospitalization",icu_and_culture_measures="study-window ICU time",support="recorded treatment within hospitalization; not baseline severity"),repeat_event_filter="not applicable to binary any-event hospitalization membership"),file.path(out,"ase_definition.json"),pretty=TRUE,auto_unbox=TRUE)
  availability("completed")
  message("ASE subgroup analysis completed.")
}
tryCatch(run_ase_analysis(),error=function(e){availability("failed");stop(e)})
