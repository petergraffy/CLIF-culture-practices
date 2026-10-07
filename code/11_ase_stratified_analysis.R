# Additional retrospective hospitalization subgroups; whole-cohort outputs are unchanged.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(lubridate);library(readr);library(ggplot2)})
source("utils/clif_io.R");source("utils/preflight.R");source("utils/ase.R");source("utils/ase_characteristics.R");source("utils/trends.R");source("utils/susceptibility.R")
out<-project_output_dir("ase");stamp<-format(Sys.time(),"%Y%m%d_%H%M%S")
write_ase<-function(x,name)write_csv(mutate(x,site_name=clif_site_name),file.path(out,paste0(name,"_",clif_site_name,"_",stamp,".csv")))
availability<-function(status,detail=NA_character_)write_ase(tibble(analysis_status=status,detail=detail),"ase_availability")

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
  # This output contains identifiers and stays private.
  write_csv(classification,file.path(private,"ase_hospitalization_classification.csv"))
  if(nrow(result$episodes))write_csv(result$episodes,file.path(private,"ase_blood_culture_criteria.csv"))
  h<-data$hospitalization %>% semi_join(data$icu_admissions,by="hospitalization_id")
  qc<-tibble(metric=c("icu_hospitalizations","classified_adult_hospitalizations","excluded_or_unclassifiable_hospitalizations","ase_hospitalizations","non_ase_hospitalizations","additional_ase_hospitalizations_with_lactate"),n=c(n_distinct(h$hospitalization_id),nrow(classification),n_distinct(h$hospitalization_id)-nrow(classification),sum(classification$ase_group=="ASE"),sum(classification$ase_group=="Non-ASE"),sum(classification$ase_with_lactate & !classification$ase_without_lactate)))
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
  write_ase(z$admission_monthly,"ase_monthly_admission_culture_proportions");write_ase(z$timing_bins,"ase_first_culture_timing_bins");write_ase(z$cohort,"ase_cohort_summary");write_ase(z$monthly,"ase_monthly_culture_rates");write_ase(z$specimen,"ase_monthly_specimen_culture_rates");write_ase(z$organisms,"ase_monthly_organism_detection_counts");write_ase(z$timing,"ase_first_culture_timing")
  write_ase(z$rows %>% filter(positive_culture) %>% group_by(ase_group,organism_group,organism_category) %>% summarise(n_detection_events=n_distinct(culture_event_id),n_hospitalizations=n_distinct(hospitalization_id),n_patients=n_distinct(patient_id),.groups="drop"),"ase_organism_distribution")
  culture_outcomes<-tibble(outcome=c("culture_collection","culture_collection","positive_culture_detection","positive_culture_detection"),count=c("n_culture_events","n_culture_events","n_positive_culture_events","n_positive_culture_events"),denominator=rep(c("n_icu_days","n_icu_admissions"),2),proportion=FALSE)
  culture<-ase_fit_screens(z$monthly,"ase_group",culture_outcomes)
  write_ase(culture$models,"ase_culture_temporal_models");if(nrow(culture$curves))write_ase(culture$curves,"ase_culture_fitted_curves")
  detection<-ase_fit_screens(z$organisms,c("ase_group","organism_category"),tibble(outcome="organism_detection",count="n_detection_events",denominator=c("n_icu_days","n_icu_admissions"),proportion=FALSE))
  write_ase(detection$models,"ase_organism_temporal_models");if(nrow(detection$curves))write_ase(detection$curves,"ase_organism_fitted_curves")
  caption<-"Retrospective hospitalization classification; primary ASE excludes lactate.\nPoints: monthly observations. Curves: season-adjusted GAM means and pointwise 95% CIs."
  points<-z$monthly %>% select(ase_group,calendar_month,culture_events_per_100_icu_days,positive_events_per_100_icu_days) %>% pivot_longer(ends_with("icu_days"),names_to="outcome",values_to="observed") %>% mutate(outcome=if_else(outcome=="culture_events_per_100_icu_days","culture_collection","positive_culture_detection"))
  p<-ggplot(points,aes(calendar_month,observed,color=ase_group))+geom_point(alpha=.4)+facet_wrap(~outcome,scales="free_y",ncol=1)+labs(x=NULL,y="Events per 100 subgroup ICU-days",color=NULL,caption=caption)+theme_bw()
  if(nrow(culture$curves)){curves<-filter(culture$curves,denominator=="n_icu_days");p<-p+geom_ribbon(data=curves,aes(y=NULL,ymin=ci_low,ymax=ci_high,fill=ase_group),alpha=.12,color=NA)+geom_line(data=curves,aes(y=fitted_value))}
  ggsave(file.path(out,paste0("ase_culture_rates_",clif_site_name,"_",stamp,".png")),p,width=12,height=8,dpi=200)
  p<-ggplot(z$monthly,aes(calendar_month,positivity,color=ase_group))+geom_point()+geom_line()+scale_y_continuous(labels=scales::percent)+labs(x=NULL,y="Positive fraction of subgroup culture events",color=NULL,caption="Observed monthly fractions; retrospective hospitalization classification, excluding lactate.")+theme_bw()
  ggsave(file.path(out,paste0("ase_culture_positivity_",clif_site_name,"_",stamp,".png")),p,width=12,height=6,dpi=200)
  top<-z$organisms %>% group_by(organism_category) %>% summarise(n=sum(n_detection_events),.groups="drop") %>% slice_max(n,n=12,with_ties=FALSE)
  if(nrow(top)) {
    points<-semi_join(z$organisms,top,by="organism_category")
    p<-ggplot(points,aes(calendar_month,detections_per_100_icu_days,color=ase_group))+geom_point(alpha=.35)+facet_wrap(~organism_category,scales="free_y",ncol=3)+labs(x=NULL,y="Detections per 100 subgroup ICU-days",color=NULL,caption=caption)+theme_bw()
    if(nrow(detection$curves)){curves<-detection$curves %>% filter(denominator=="n_icu_days") %>% semi_join(top,by="organism_category");p<-p+geom_ribbon(data=curves,aes(y=NULL,ymin=ci_low,ymax=ci_high,fill=ase_group),alpha=.12,color=NA)+geom_line(data=curves,aes(y=fitted_value))}
    ggsave(file.path(out,paste0("ase_organism_detection_rates_",clif_site_name,"_",stamp,".png")),p,width=14,height=11,dpi=200)
  }
  p<-ggplot(z$timing,aes(ase_group,percent_with_icu_culture,fill=ase_group))+geom_col()+labs(x=NULL,y="ICU admissions with a culture (%)",caption="New ICU admissions only; includes admissions with no ICU culture.")+theme_bw()+theme(legend.position="none")
  ggsave(file.path(out,paste0("ase_admissions_cultured_",clif_site_name,"_",stamp,".png")),p,width=7,height=5,dpi=200)
  top_specimens<-z$specimen %>% group_by(fluid_category) %>% summarise(n=sum(n_culture_events),.groups="drop") %>% slice_max(n,n=6,with_ties=FALSE)
  specimens<-semi_join(z$specimen,top_specimens,by="fluid_category")
  for(metric in c("cultures_per_100_icu_days","positivity")) {
    p<-ggplot(specimens,aes(calendar_month,.data[[metric]],color=ase_group))+geom_point(alpha=.5)+geom_line()+facet_wrap(~fluid_category,scales="free_y",ncol=2)+labs(x=NULL,y=if(metric=="positivity")"Positive fraction" else "Cultures per 100 subgroup ICU-days",color=NULL,caption="Observed monthly values; identical specimen categories selected across groups.")+theme_bw()
    if(metric=="positivity")p<-p+scale_y_continuous(labels=scales::percent)
    ggsave(file.path(out,paste0("ase_specimen_",metric,"_",clif_site_name,"_",stamp,".png")),p,width=13,height=10,dpi=200)
  }
  ast_path<-find_table_path("microbiology_susceptibility",required=FALSE)
  if(is.na(ast_path))write_ase(tibble(analysis_status="skipped_table_unavailable"),"ase_susceptibility_availability") else {
    ast<-read_tbl("microbiology_susceptibility");monthly_ast<-qc_ast<-list()
    coverage<-as.numeric(Sys.getenv("AST_MIN_TESTING_FRACTION","0.5"));linkage<-as.numeric(Sys.getenv("AST_MIN_LINKAGE_FRACTION","0.9"))
    if(any(!is.finite(c(coverage,linkage)))||any(c(coverage,linkage)<0|c(coverage,linkage)>1))stop("Invalid AST coverage/linkage thresholds")
    for(g in c("ASE","Non-ASE")) {
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
        p<-ggplot(observed,aes(calendar_month,observed,color=ase_group))+geom_point(aes(shape=eligible),alpha=.5,na.rm=TRUE)+facet_wrap(~pair,scales=if(fraction)"fixed" else "free_y",ncol=2)+labs(x=NULL,y=if(fraction)"Non-susceptible fraction among interpretable tests" else "Detections per 100 subgroup ICU-days",color=NULL,shape="Model eligible",caption=caption)+theme_bw()
        if(fraction)p<-p+scale_y_continuous(labels=scales::percent,limits=c(0,1))
        if(nrow(a$curves)) {
          fitted<-a$curves %>% filter(specimen_stratum=="Overall",.data$outcome==.env$outcome) %>% semi_join(top_pairs,by=c("organism_category","antimicrobial_category")) %>% mutate(pair=paste(organism_category,antimicrobial_category,sep=" / "))
          p<-p+geom_ribbon(data=fitted,aes(y=NULL,ymin=ci_low,ymax=ci_high,fill=ase_group),color=NA,alpha=.12)+geom_line(data=fitted,aes(y=fitted_value))
        }
        ggsave(file.path(out,paste0("ase_",outcome,"_",clif_site_name,"_",stamp,".png")),p,width=13,height=10,dpi=200)
      }
      coverage_points<-points %>% pivot_longer(c(linkage_fraction,testing_fraction),names_to="coverage_metric",values_to="fraction")
      p<-ggplot(coverage_points,aes(calendar_month,fraction,color=ase_group,shape=coverage_metric))+geom_point(na.rm=TRUE)+facet_wrap(~pair,ncol=2)+scale_y_continuous(limits=c(0,1),labels=scales::percent)+geom_hline(yintercept=c(coverage,linkage),linetype=2,color="grey60")+labs(x=NULL,y="Linkage / interpretable testing coverage",color=NULL,shape=NULL)+theme_bw()
      ggsave(file.path(out,paste0("ase_susceptibility_coverage_",clif_site_name,"_",stamp,".png")),p,width=13,height=10,dpi=200)
      write_ase(tibble(analysis_status="completed",n_estimated_models=sum(a$models$model_status=="estimated")),"ase_susceptibility_availability")
    } else write_ase(tibble(analysis_status="no_linked_interpretable_analysis"),"ase_susceptibility_availability")
  }
  upstream<-jsonlite::fromJSON(project_path("utils","ase","upstream.json"))
  jsonlite::write_json(list(classification="any qualifying ASE episode during full hospitalization",primary_include_lactate=FALSE,adult_minimum_age=18,calendar_timezone="UTC",denominators="subgroup ICU-days and new ICU admissions; includes uncultured stays",upstream=upstream,characteristics=list(analysis_unit="hospitalization",categorical_denominator="all subgroup hospitalizations including missing",continuous_summary="median and quartiles among observed values",hospitalization_outcomes="full hospitalization",icu_and_culture_measures="study-window ICU time",support="recorded treatment within hospitalization; not baseline severity"),repeat_event_filter="not applicable to binary any-event hospitalization membership"),file.path(out,"ase_definition.json"),pretty=TRUE,auto_unbox=TRUE)
  availability("completed")
  message("ASE subgroup analysis completed.")
}
tryCatch(run_ase_analysis(),error=function(e){availability("failed");stop(e)})
