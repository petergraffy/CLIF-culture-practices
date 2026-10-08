# Retrospective first-ASE onset alignment. Hospital observation includes non-ICU time.
build_ase_onset_preview <- function(events, anchors, hospitalizations, defining, start, end) {
  a<-anchors %>% filter(timing_anchor=="Earliest qualifying organ dysfunction",!is.na(event_time)) %>% mutate(onset_stratum=case_when(event_time<=first_icu_in~"At/before first ICU entry",!is.na(icu_in_dttm)~"During ICU stay",TRUE~NA_character_))
  qc<-a %>% count(onset_stratum,name="n_ase_hospitalizations",.drop=FALSE)
  a<-a %>% filter(!is.na(onset_stratum)) %>% select(hospitalization_id,event_time,onset_stratum) %>% inner_join(hospitalizations %>% select(hospitalization_id,admission_dttm,discharge_dttm),by="hospitalization_id")
  if(anyDuplicated(a$hospitalization_id))stop("Onset alignment requires one anchor per hospitalization")
  strata<-c("At/before first ICU entry","During ICU stay")
  bins<-tibble(bin_start_hour=seq(-48,71,1),bin_end_hour=seq(-47,72,1))
  exposure<-crossing(a,bins) %>% mutate(hours=pmax(0,(pmin(as.numeric(discharge_dttm),as.numeric(end),as.numeric(event_time)+bin_end_hour*3600)-pmax(as.numeric(admission_dttm),as.numeric(start),as.numeric(event_time)+bin_start_hour*3600))/3600)) %>% group_by(onset_stratum,bin_start_hour,bin_end_hour) %>% summarise(observed_patient_hours=sum(hours),n_hospitalizations_observed=sum(hours>0),.groups="drop")
  n<-a %>% count(onset_stratum,name="n_anchor_hospitalizations")
  exposure<-crossing(onset_stratum=strata,bins) %>% left_join(exposure,by=c("onset_stratum","bin_start_hour","bin_end_hour")) %>% left_join(n,by="onset_stratum") %>% mutate(across(c(observed_patient_hours,n_hospitalizations_observed,n_anchor_hospitalizations),~coalesce(.x,0)))
  e<-events %>% inner_join(a,by="hospitalization_id") %>% mutate(relative_hour=as.numeric(difftime(collect_dttm,event_time,units="hours"))) %>% filter(relative_hour>=-48,relative_hour<72) %>% mutate(bin_start_hour=floor(relative_hour),bin_end_hour=bin_start_hour+1,fluid_category=coalesce(fluid_category,"missing")) %>% left_join(defining %>% distinct(hospitalization_id,collect_dttm) %>% mutate(defining_blood_collection=TRUE),by=c("hospitalization_id","collect_dttm")) %>% mutate(defining_blood_collection=coalesce(defining_blood_collection,FALSE)&fluid_category=="blood_buffy")
  e<-bind_rows(e %>% mutate(sensitivity="All collections"),e %>% filter(!defining_blood_collection) %>% mutate(sensitivity="Exclude linked ASE blood collection"))
  e<-bind_rows(e %>% mutate(specimen="All cultures"),e %>% mutate(specimen=fluid_category))
  counts<-e %>% count(onset_stratum,sensitivity,specimen,bin_start_hour,bin_end_hour,result_status,name="n_culture_events")
  result<-crossing(exposure,sensitivity=c("All collections","Exclude linked ASE blood collection"),specimen=unique(c("All cultures",events$fluid_category[!is.na(events$fluid_category)],"missing")),result_status=culture_result_levels) %>% left_join(counts,by=c("onset_stratum","sensitivity","specimen","bin_start_hour","bin_end_hour","result_status")) %>% mutate(n_culture_events=coalesce(n_culture_events,0L),collections_per_100_patient_hours=if_else(observed_patient_hours>0,100*n_culture_events/observed_patient_hours,NA_real_))
  yield<-result %>% group_by(onset_stratum,sensitivity,specimen,bin_start_hour,bin_end_hour) %>% summarise(n_positive=sum(n_culture_events[result_status=="Positive"]),n_culture_events=sum(n_culture_events),observed_patient_hours=first(observed_patient_hours),n_hospitalizations_observed=first(n_hospitalizations_observed),.groups="drop") %>% mutate(percent_positive=if_else(n_culture_events>0,100*n_positive/n_culture_events,NA_real_))
  list(results=result,yield=yield,exposure=exposure,qc=qc)
}

plot_ase_onset_preview <- function(x,out,site) {
  color<-c("Negative/no growth"="#80B1D3","Positive"="#D55E00","Mixed/contaminated"="#CC79A7","Indeterminate"="#BBBBBB")
  strata<-c("At/before first ICU entry","During ICU stay")
  label<-function(v)stringr::str_wrap(gsub("_"," ",v),25)
  for(sensitivity in unique(x$results$sensitivity)) {
    suffix<-if(sensitivity=="All collections")"all" else "exclude_linked_blood"
    d<-x$results %>% filter(specimen=="All cultures",.data$sensitivity==.env$sensitivity) %>% mutate(onset_stratum=factor(onset_stratum,levels=strata),result_status=factor(result_status,levels=names(color)))
    p<-ggplot(d,aes((bin_start_hour+bin_end_hour)/2,collections_per_100_patient_hours,fill=result_status))+geom_col(width=1)+geom_vline(xintercept=0,linetype=2)+facet_wrap(~onset_stratum,ncol=2)+scale_fill_manual(values=color,drop=FALSE,name="Eventual culture result")+labs(x="Hours relative to earliest qualifying organ-dysfunction proxy",y="Collections per 100 observed hospital patient-hours",title=paste(site,"—",sensitivity),caption=NULL)+theme_bw()+theme(legend.position="bottom")
    ggsave(file.path(out,paste0("ase_onset_collection_results_",site,"_",suffix,".png")),p,width=13,height=7,dpi=200)
    d<-x$yield %>% filter(specimen %in% c("All cultures","blood_buffy","respiratory_tract","genito_urinary_tract"),.data$sensitivity==.env$sensitivity) %>% mutate(onset_stratum=factor(onset_stratum,levels=strata),specimen=factor(specimen,levels=c("All cultures","blood_buffy","respiratory_tract","genito_urinary_tract")))
    p<-ggplot(d %>% mutate(percent_positive=if_else(n_culture_events>=20,percent_positive,NA_real_)),aes((bin_start_hour+bin_end_hour)/2,percent_positive))+geom_line(na.rm=TRUE,color="#D55E00")+geom_point(na.rm=TRUE,color="#D55E00")+geom_vline(xintercept=0,linetype=2)+facet_grid(specimen~onset_stratum,labeller=labeller(specimen=label))+labs(x="Hours relative to earliest qualifying organ-dysfunction proxy",y="Collected cultures with a positive result (%)",title=paste(site,"—",sensitivity),caption=NULL)+theme_bw()
    ggsave(file.path(out,paste0("ase_onset_positivity_",site,"_",suffix,".png")),p,width=13,height=10,dpi=200)
  }
}
