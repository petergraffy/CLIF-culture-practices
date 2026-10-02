# Optional standardized susceptibility analysis; no raw-text resistance substitution.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(lubridate);library(ggplot2);library(glue)})
source("utils/clif_io.R");source("utils/susceptibility.R");source("utils/trends.R")
out_dir <- project_output_dir("susceptibility");stamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
path <- find_table_path("microbiology_susceptibility",required=FALSE)
write_availability <- function(status,n_estimated=0L) {
  write_csv(tibble(site_name=clif_site_name,susceptibility_table_available=!is.na(path),analysis_status=status,n_estimated_models=n_estimated),file.path(out_dir,glue("susceptibility_availability_{clif_site_name}_{stamp}.csv")))
}
run_susceptibility <- function() {
  if(is.na(path)){write_availability("skipped_table_unavailable");message("Susceptibility table unavailable; analysis skipped.");return(invisible(NULL))}
  write_availability("analysis_started")
  coverage_min <- as.numeric(Sys.getenv("AST_MIN_TESTING_FRACTION","0.5"));linkage_min <- as.numeric(Sys.getenv("AST_MIN_LINKAGE_FRACTION","0.9"))
  if(any(!is.finite(c(coverage_min,linkage_min))) || any(c(coverage_min,linkage_min)<0 | c(coverage_min,linkage_min)>1))stop("AST coverage/linkage thresholds must be between 0 and 1.")
  start <- safe_ts(study_settings$study_start_date);end <- safe_ts(study_settings$study_end_date)+days(1)
  rows <- read_clif_csv(latest_project_intermediate_file("^icu_culture_rows_.*\\.csv$","cohort"),show_col_types=FALSE) %>% mutate(collect_dttm=safe_ts(collect_dttm),organism_id=as.character(organism_id))
  stays <- read_culture_data(start,end)$icu_admissions
  month_seq <- seq(floor_date(min(stays$icu_in_dttm_clipped),"month"),floor_date(max(stays$icu_out_dttm_clipped-seconds(1)),"month"),by="month")
  result <- build_susceptibility_analysis(rows,read_tbl("microbiology_susceptibility"),monthly_icu_denominators(stays,month_seq),project_path("config","mcide"),coverage_min,linkage_min)
  write_csv(mutate(result$qc,site_name=clif_site_name),file.path(out_dir,glue("susceptibility_qc_{clif_site_name}_{stamp}.csv")))
  write_csv(mutate(result$linkage_qc,site_name=clif_site_name),file.path(out_dir,glue("monthly_susceptibility_linkage_qc_{clif_site_name}_{stamp}.csv")))
  if(!nrow(result$monthly)){write_availability(result$analysis_status);message("Susceptibility: ",result$analysis_status);return(invisible(NULL))}
  monthly <- result$monthly %>% mutate(site_name=clif_site_name)
  write_csv(monthly,file.path(out_dir,glue("monthly_organism_antimicrobial_susceptibility_{clif_site_name}_{stamp}.csv")))
  # Fit each outcome once; summaries and plotted predictions come from the same fit.
  keys <- c("organism_category","antimicrobial_category","specimen_stratum")
  models <- curves <- list()
  screens <- distinct(monthly,across(all_of(keys)))
  for(i in seq_len(nrow(screens))) {
    dat <- semi_join(monthly,screens[i,],by=keys)
    for(outcome in c("susceptible_detection_rate","non_susceptible_detection_rate","non_susceptible_fraction")) {
      fraction <- outcome=="non_susceptible_fraction"
      selected <- if(fraction)filter(dat,fraction_model_eligible) else filter(dat,rate_model_eligible)
      count <- if(outcome=="susceptible_detection_rate")"n_susceptible" else "n_non_susceptible"
      fit <- fit_temporal_model(selected,if(fraction)"n_interpretable" else "n_icu_days",count,fraction)
      models[[length(models)+1]] <- bind_cols(screens[i,],fit$summary) %>% mutate(outcome=outcome,n_usable_months=nrow(selected),first_usable_month=if(nrow(selected))min(selected$calendar_month) else as.POSIXct(NA),last_usable_month=if(nrow(selected))max(selected$calendar_month) else as.POSIXct(NA),minimum_testing_fraction_for_rate=coverage_min,minimum_linkage_fraction=linkage_min)
      if(nrow(fit$predictions))curves[[length(curves)+1]]<-bind_cols(screens[rep(i,nrow(fit$predictions)),],fit$predictions) %>% mutate(outcome=outcome)
    }
  }
  fits <- bind_rows(models) %>% group_by(outcome,specimen_stratum) %>% mutate(fdr_p_value=p.adjust(p_value,"BH"),direction=classify_temporal_direction(model_status,fdr_p_value,annual_percent_change,residual_dependence_flag)) %>% ungroup() %>% mutate(site_name=clif_site_name,effect_scale=if_else(outcome=="non_susceptible_fraction","annualized endpoint odds ratio","annualized endpoint rate ratio"))
  curves <- bind_rows(curves)
  write_csv(fits,file.path(out_dir,glue("susceptibility_temporal_models_{clif_site_name}_{stamp}.csv")))
  if(nrow(curves))write_csv(curves,file.path(out_dir,glue("susceptibility_gam_fitted_curves_{clif_site_name}_{stamp}.csv")))
  n_estimated<-sum(fits$model_status=="estimated")
  write_availability(if(result$analysis_status=="no_interpretable_tests")result$analysis_status else if(n_estimated)"completed" else "insufficient_model_data",n_estimated)
  top_pairs <- monthly %>% filter(specimen_stratum=="Overall") %>% group_by(organism_category,antimicrobial_category) %>% summarise(n_tested=sum(n_interpretable),.groups="drop") %>% slice_max(n_tested,n=12,with_ties=FALSE)
  points <- monthly %>% filter(specimen_stratum=="Overall") %>% semi_join(top_pairs,by=c("organism_category","antimicrobial_category")) %>% mutate(pair=gsub("_"," ",paste(organism_category,antimicrobial_category,sep=" / ")))
  for(outcome in c("susceptible_detection_rate","non_susceptible_detection_rate","non_susceptible_fraction")) {
    fraction<-outcome=="non_susceptible_fraction"
    observed<-points %>% mutate(observed_value=if(fraction)non_susceptible_fraction else if(outcome=="susceptible_detection_rate")susceptible_per_100_icu_days else non_susceptible_per_100_icu_days,model_eligible=if(fraction)fraction_model_eligible else rate_model_eligible)
    p<-ggplot(observed,aes(calendar_month,observed_value))+geom_point(aes(shape=model_eligible),alpha=0.6,na.rm=TRUE)+facet_wrap(~pair,scales=if(fraction)"fixed" else "free_y",ncol=3)+labs(x=NULL,y=if(fraction)"Non-susceptible fraction among interpretable tests" else paste(if(outcome=="susceptible_detection_rate")"Susceptible" else "Non-susceptible","detections per 100 ICU days"),shape="Model eligible",caption="Points: observed values. Line and ribbon: season-adjusted long-term mean and pointwise 95% CI.\nValidated zero-detection months retained; unavailable testing remains unknown. See coverage figure.")+theme_bw()+theme(legend.position="bottom")
    if(fraction)p<-p+scale_y_continuous(limits=c(0,1),labels=scales::percent)
    if(nrow(curves)) {
      fit_data<-curves %>% filter(specimen_stratum=="Overall",.data$outcome==.env$outcome) %>% semi_join(top_pairs,by=c("organism_category","antimicrobial_category")) %>% mutate(pair=gsub("_"," ",paste(organism_category,antimicrobial_category,sep=" / ")))
      p<-p+geom_ribbon(data=fit_data,mapping=aes(x=calendar_month,ymin=ci_low,ymax=ci_high),inherit.aes=FALSE,alpha=0.15)+geom_line(data=fit_data,aes(y=fitted_value),color="#9C3333",na.rm=TRUE)
    }
    ggsave(file.path(out_dir,glue("{outcome}_{clif_site_name}_{stamp}.png")),p,width=14,height=10,dpi=200)
  }
  coverage<-points %>% pivot_longer(c(linkage_fraction,testing_fraction,testing_fraction_linkable),names_to="metric",values_to="fraction")
  p<-ggplot(coverage,aes(calendar_month,fraction,color=metric))+geom_point(na.rm=TRUE)+facet_wrap(~pair,ncol=3)+geom_hline(yintercept=coverage_min,linetype=2,color="grey50")+geom_hline(yintercept=linkage_min,linetype=3,color="grey50")+scale_y_continuous(limits=c(0,1),labels=scales::percent)+labs(x=NULL,y="Isolate linkage / interpretable testing coverage",color=NULL,caption=glue("Testing uses all observed positive isolates; linkable-only coverage shown for comparison.\nDashed: minimum testing {coverage_min}; dotted: minimum linkage {linkage_min}. Missing values are undefined, not zero."))+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out_dir,glue("susceptibility_testing_coverage_{clif_site_name}_{stamp}.png")),p,width=14,height=10,dpi=200)
}
tryCatch(run_susceptibility(),error=function(e){write_availability("failed");stop(e)})
