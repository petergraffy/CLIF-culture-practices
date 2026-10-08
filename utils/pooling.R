# Aggregate-only cross-site analysis. No patient/isolate identifiers are accepted.
pooling_forbidden <- function(x) {
  forbidden <- c("patient_id", "hospitalization_id", "organism_id", "icu_stay_id", "culture_event_id", "detection_event_id", "culture_row_id")
  if (any(names(x) %in% forbidden) || any(grepl("_dttm$", names(x)))) stop("Pooling input contains private identifiers or exact clinical timestamps.")
}

one_export <- function(run_dir, folder, pattern, required = TRUE) {
  paths <- list.files(file.path(run_dir, folder), pattern = pattern, full.names = TRUE)
  if (length(paths) > 1) stop("Ambiguous aggregate export: ", pattern)
  if (!length(paths)) {
    if (required) stop("Missing aggregate export: ", pattern, "; rerun the site's current pipeline.")
    return(NULL)
  }
  x <- readr::read_csv(paths, show_col_types = FALSE)
  pooling_forbidden(x)
  x
}

read_pooling_sites <- function(registry) {
  required <- c("site_name", "run_dir", "validated_start_date", "validated_end_date", "culture_qc_pass", "ast_qc_pass")
  if (!all(required %in% names(registry))) stop("Pooling registry missing required columns.")
  if (!nrow(registry)) stop("No sites registered. Complete the pooling registry with reviewed aggregate runs.")
  if (anyNA(registry[required]) || anyDuplicated(registry$site_name) || anyDuplicated(registry$run_dir)) stop("Registry has missing values or duplicate sites/runs.")
  if (any(!nzchar(registry$site_name)) || any(!nzchar(registry$run_dir))) stop("Registry site names and run directories must be nonempty.")
  if (any(!registry$culture_qc_pass %in% c(TRUE,FALSE)) || any(!registry$ast_qc_pass %in% c(TRUE,FALSE))) stop("QC fields must be TRUE/FALSE.")
  culture <- ast <- denominators <- list(); audit <- hashes <- list(); reference_hashes <- NULL; source_paths <- character()
  for (i in seq_len(nrow(registry))) {
    r <- registry[i, ]
    manifest_path <- file.path(r$run_dir, "provenance", "run_manifest.json")
    m <- jsonlite::fromJSON(manifest_path)
    if (!identical(m$analysis_status, "completed") || !identical(m$site_name, r$site_name)) stop("Run incomplete or site does not match manifest: ", r$site_name)
    if (is.null(m$code_and_mapping_md5) || is.null(m$study_start_date) || is.null(m$study_end_date)) stop("Run lacks required provenance.")
    comparable <- c("utils/study_settings.R","utils/culture_core.R","utils/trends.R","utils/susceptibility.R","code/08_organism_trends.R","code/09_susceptibility_trends.R","config/mcide/clif_microbiology_susceptibility_category.csv","config/mcide/clif_microbiology_susceptibility_antibiotics_category.csv")
    source_hashes <- unlist(m$code_and_mapping_md5[comparable])
    if(length(source_hashes)!=length(comparable) || anyNA(source_hashes))stop("Run lacks comparable source/schema hashes: ",r$site_name)
    if(r$culture_qc_pass) {
      if (!identical(unname(source_hashes), unname(tools::md5sum(comparable)))) stop("Central and site analysis/schema code differs; rerun sites with current code.")
      if(is.null(reference_hashes))reference_hashes<-source_hashes else if(!identical(source_hashes,reference_hashes))stop("Site runs use different analysis or mCIDE versions; rerun with matching code.")
    }
    start <- as.Date(r$validated_start_date); end <- as.Date(r$validated_end_date)
    if (is.na(start) || is.na(end) || start > end || start < as.Date(m$study_start_date) || end > as.Date(m$study_end_date)) stop("Validated interval lies outside run's study interval: ", r$site_name)
    # Retain only complete calendar months, never partial-month denominators.
    first <- lubridate::ceiling_date(start, "month", change_on_boundary = FALSE)
    last <- lubridate::floor_date(end + 1, "month") - lubridate::period(month=1)
    export_audit <- one_export(r$run_dir,"site_exports","^site_export_privacy_audit_.*[.]csv$")
    if (nrow(export_audit)) stop("Run export audit is not clean: ", r$site_name)
    ast_paths <- list.files(file.path(r$run_dir,"susceptibility"),pattern="^monthly_organism_antimicrobial_susceptibility_.*[.]csv$",full.names=TRUE)
    audit[[i]] <- tibble::tibble(site_name=r$site_name,run_id=m$run_id,culture_qc_pass=r$culture_qc_pass,ast_qc_pass=r$ast_qc_pass,ast_available=length(ast_paths)>0,first_full_month=first,last_full_month=last)
    input_paths <- c(manifest_path,list.files(file.path(r$run_dir,"organism_trends"),pattern="^monthly_(all_organism_detection_counts|icu_denominators_for_organism_trends)_.*[.]csv$",full.names=TRUE),ast_paths,list.files(file.path(r$run_dir,"site_exports"),pattern="^site_export_privacy_audit_.*[.]csv$",full.names=TRUE))
    source_paths <- c(source_paths,input_paths)
    hashes[[i]]<-tibble::tibble(site_name=r$site_name,run_id=m$run_id,file_basename=basename(input_paths),md5=unname(tools::md5sum(input_paths)))
    if (!r$culture_qc_pass) next
    den <- one_export(r$run_dir,"organism_trends","^monthly_icu_denominators_for_organism_trends_.*[.]csv$")
    if(!all(c("calendar_month","n_icu_days","n_icu_admissions","n_observed_culture_events") %in% names(den)))stop("Monthly denominator export schema invalid.")
    den <- den %>% dplyr::mutate(calendar_month=as.Date(calendar_month),site_name=r$site_name) %>% dplyr::filter(calendar_month>=first,calendar_month<=last)
    if (anyDuplicated(den$calendar_month) || anyNA(den$calendar_month) || any(lubridate::day(den$calendar_month)!=1)) stop("Invalid or duplicate calendar months.")
    denominators[[r$site_name]] <- den
    x <- one_export(r$run_dir,"organism_trends","^monthly_all_organism_detection_counts_.*[.]csv$")
    if (!all(c("organism_category","n_detection_events","calendar_month") %in% names(x))) stop("All-organism export schema invalid.")
    if ("site_name" %in% names(x) && any(is.na(x$site_name) | x$site_name!=r$site_name)) stop("Export site does not match registry.")
    x <- x %>% dplyr::mutate(calendar_month=as.Date(calendar_month),site_name=r$site_name)
    if(anyNA(x$calendar_month) || any(lubridate::day(x$calendar_month)!=1))stop("Invalid calendar months.")
    x <- x %>% dplyr::filter(calendar_month>=first,calendar_month<=last)
    if (anyDuplicated(x[c("calendar_month","organism_category")]) || anyNA(x$organism_category)) stop("Duplicate or missing organism category.")
    culture[[r$site_name]] <- x %>% dplyr::select(site_name,calendar_month,organism_category,n_detection_events)
    if (r$ast_qc_pass) {
      x <- one_export(r$run_dir,"susceptibility","^monthly_organism_antimicrobial_susceptibility_.*[.]csv$",FALSE)
      if (!is.null(x)) {
        if(!all(c("organism_category","antimicrobial_category","specimen_stratum","n_susceptible","n_non_susceptible","n_interpretable","testing_fraction","calendar_month","n_icu_days","n_culture_isolates","n_linkable_culture_isolates","n_observed_culture_events","n_positive_rows_missing_organism_category") %in% names(x)))stop("Susceptibility export schema invalid.")
        if("site_name" %in% names(x) && any(is.na(x$site_name) | x$site_name!=r$site_name))stop("AST export site does not match registry.")
        x<-x %>% dplyr::mutate(calendar_month=as.Date(calendar_month),site_name=r$site_name)
        if(anyNA(x$calendar_month) || any(lubridate::day(x$calendar_month)!=1))stop("Invalid AST calendar months.")
        if(any(is.finite(x$testing_fraction) & (x$testing_fraction<0 | x$testing_fraction>1)))stop("Invalid AST testing coverage.")
        x<-x %>% dplyr::filter(calendar_month>=first,calendar_month<=last)
        if(anyNA(x[c("organism_category","antimicrobial_category","specimen_stratum")]))stop("Missing AST screen keys.")
        if(any(!is.finite(x$n_susceptible)|!is.finite(x$n_non_susceptible)|!is.finite(x$n_interpretable)|x$n_susceptible+x$n_non_susceptible!=x$n_interpretable))stop("Invalid susceptibility counts.")
        matched <- x %>% dplyr::left_join(dplyr::select(den,calendar_month,validated_icu_days=n_icu_days),by="calendar_month")
        if(any(!is.finite(matched$validated_icu_days) | !is.finite(matched$n_icu_days) | abs(matched$n_icu_days-matched$validated_icu_days)>1e-8))stop("AST denominators do not match culture denominators.")
        source_check <- x %>% dplyr::left_join(dplyr::select(den,calendar_month,total_culture_events=n_observed_culture_events),by="calendar_month")
        if(any(!is.finite(source_check$n_observed_culture_events) | source_check$n_observed_culture_events<0 | source_check$n_observed_culture_events>source_check$total_culture_events))stop("Invalid specimen-specific source coverage.")
        ast[[r$site_name]]<-x
      }
    }
  }
  den <- dplyr::bind_rows(denominators)
  if (!nrow(den)) stop("No culture-QC-approved sites with complete months.")
  if (any(!is.finite(den$n_icu_days) | den$n_icu_days < 0 | !is.finite(den$n_icu_admissions) | den$n_icu_admissions < 0 | !is.finite(den$n_observed_culture_events) | den$n_observed_culture_events < 0)) stop("Invalid monthly denominators or source coverage counts.")
  counts <- dplyr::bind_rows(culture)
  # Absent categories are true zero only because these exports enumerate ALL detections.
  grid <- tidyr::crossing(den,organism_category=sort(unique(counts$organism_category))) %>% dplyr::left_join(counts,by=c("site_name","calendar_month","organism_category")) %>% dplyr::mutate(n_detection_events=dplyr::coalesce(n_detection_events,0L))
  list(culture=grid,ast=dplyr::bind_rows(ast),audit=dplyr::bind_rows(audit),hashes=dplyr::bind_rows(hashes),source_paths=source_paths)
}

common_pooling_months <- function(data, denominator, count, proportion=FALSE, coverage_min=0.5,linkage_min=0.9) {
  pooling_forbidden(data)
  if (anyDuplicated(data[c("site_name","calendar_month")])) stop("Duplicate site/month rows within pooling screen.")
  if (any(!is.finite(data[[count]]) | data[[count]] < 0 | data[[count]] != round(data[[count]]))) stop("Invalid monthly counts.")
  if (proportion && any(data[[count]]>data[[denominator]],na.rm=TRUE)) stop("Non-susceptible count exceeds tested denominator.")
  sites <- unique(data$site_name)
  x <- data %>% dplyr::filter(is.finite(.data[[denominator]]),.data[[denominator]]>0)
  if ("n_observed_culture_events" %in% names(x)) x <- dplyr::filter(x,n_observed_culture_events>0)
  if ("testing_fraction" %in% names(x)) {
    x <- screen_susceptibility_months(x,coverage_min,linkage_min)
    x <- if(proportion) dplyr::filter(x,fraction_model_eligible) else dplyr::filter(x,rate_model_eligible)
  }
  months <- x %>% dplyr::count(calendar_month) %>% dplyr::filter(n==length(sites)) %>% dplyr::pull(calendar_month)
  x %>% dplyr::filter(calendar_month %in% months) %>% dplyr::arrange(site_name,calendar_month)
}

meta_pool <- function(effects) {
  valid <- effects %>% dplyr::filter(model_status=="estimated",is.finite(log_annual_effect),is.finite(log_annual_se),log_annual_se>0)
  base <- tibble::tibble(n_sites=nrow(valid),model_status="insufficient_sites",annual_effect=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_value=NA_real_,tau2=NA_real_,I2=NA_real_,Q_p_value=NA_real_,prediction_low=NA_real_,prediction_high=NA_real_,inference_status="not_estimated",model_warning=NA_character_)
  if (nrow(valid)<2) return(base)
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(metafor::rma.uni(yi=valid$log_annual_effect,sei=valid$log_annual_se,method="REML",test="adhoc"),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
  if (inherits(fit,"error")) {base$model_status<-paste0("failed: ",conditionMessage(fit));return(base)}
  pi <- predict(fit,transf=exp)
  base %>% dplyr::mutate(model_status="estimated",annual_effect=exp(as.numeric(fit$b)),ci_low=exp(fit$ci.lb),ci_high=exp(fit$ci.ub),p_value=fit$pval,tau2=fit$tau2,I2=fit$I2,Q_p_value=fit$QEp,prediction_low=if(nrow(valid)>=3)pi$pi.lb else NA_real_,prediction_high=if(nrow(valid)>=3)pi$pi.ub else NA_real_,inference_status=dplyr::case_when(any(is.na(valid$residual_dependence_flag) | valid$residual_dependence_flag)~"exploratory_residual_dependence",nrow(valid)<3~"exploratory_few_sites",TRUE~"eligible"),model_warning=if(length(warnings))paste(unique(warnings),collapse="; ") else NA_character_)
}

fit_joint_pool <- function(data,denominator,count,proportion=FALSE) {
  empty <- function(status) list(summary=tibble::tibble(model_status=status,n_sites=dplyr::n_distinct(data$site_name),theta=NA_real_,long_term_edf=NA_real_,residual_dependence_flag=NA,model_warning=NA_character_),curves=tibble::tibble(),diagnostics=tibble::tibble())
  if (dplyr::n_distinct(data$site_name)<3) return(empty("insufficient_sites"))
  dat <- data %>% dplyr::mutate(site=factor(site_name),y=.data[[count]],denominator=.data[[denominator]],time_years=as.numeric(as.Date(calendar_month)-min(as.Date(calendar_month)))/365.25,season_sin=sin(2*pi*(lubridate::month(calendar_month)-1)/12),season_cos=cos(2*pi*(lubridate::month(calendar_month)-1)/12),log_denominator=log(denominator))
  k <- min(8L,max(4L,floor(dplyr::n_distinct(dat$calendar_month)/12)))
  s <- mgcv::s
  formula <- if(proportion) cbind(y,denominator-y) ~ site+s(time_years,k=k,bs="cr")+s(time_years,site,k=k,bs="sz")+site:season_sin+site:season_cos else y ~ site+s(time_years,k=k,bs="cr")+s(time_years,site,k=k,bs="sz")+site:season_sin+site:season_cos+offset(log_denominator)
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(mgcv::gam(formula,data=dat,family=if(proportion)quasibinomial() else mgcv::nb(),method="REML"),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
  if(inherits(fit,"error"))return(empty(paste0("failed: ",conditionMessage(fit))))
  if(!isTRUE(fit$converged))return(empty("not_converged"))
  dat$residual <- as.numeric(residuals(fit,type="pearson"))
  diagnostics <- dat %>% dplyr::group_by(site_name) %>% dplyr::group_modify(function(.x,.y){
    z <- dplyr::arrange(.x,calendar_month);adj <- which(diff(lubridate::year(z$calendar_month)*12+lubridate::month(z$calendar_month))==1)
    r <- if(length(adj)>=12 && sd(z$residual[adj])>0 && sd(z$residual[adj+1])>0)cor(z$residual[adj],z$residual[adj+1]) else NA_real_
    tibble::tibble(residual_lag1=r,residual_dependence_flag=if(is.finite(r))abs(r)>0.3 else NA)
  }) %>% dplyr::ungroup()
  pred <- dat;pred$log_denominator<-0;pred$season_sin<-0;pred$season_cos<-0
  xp <- predict(fit,pred,type="lpmatrix");vc <- vcov(fit,unconditional=TRUE);eta <- as.numeric(xp%*%coef(fit));se <- sqrt(pmax(0,rowSums((xp%*%vc)*xp)))
  inverse <- if(proportion)plogis else function(z)100*exp(z)
  curves <- tibble::tibble(site_name=dat$site_name,calendar_month=dat$calendar_month,fitted_value=inverse(eta),ci_low=inverse(eta-1.96*se),ci_high=inverse(eta+1.96*se),curve_type="site")
  # Arithmetic response-scale mean with equal, fixed site weights. Delta-method covariance
  # uses the full joint coefficient covariance; sites are not assumed independent here.
  pooled <- lapply(sort(unique(dat$calendar_month)),function(month){
    i <- which(dat$calendar_month==month);values <- inverse(eta[i]);mean_value<-mean(values)
    deriv <- if(proportion)values*(1-values) else values
    gradient <- colMeans(xp[i,,drop=FALSE]*deriv)
    response_se <- sqrt(as.numeric(gradient%*%vc%*%gradient))
    link_mean <- if(proportion)qlogis(mean_value) else log(mean_value)
    link_se <- response_se/if(proportion)(mean_value*(1-mean_value)) else mean_value
    trans <- if(proportion)plogis else exp
    tibble::tibble(site_name="Equal-site mean",calendar_month=month,fitted_value=mean_value,ci_low=trans(link_mean-1.96*link_se),ci_high=trans(link_mean+1.96*link_se),curve_type="pooled_equal_site")
  })
  list(summary=tibble::tibble(model_status="estimated",n_sites=nlevels(dat$site),theta=if(proportion)NA_real_ else fit$family$getTheta(TRUE),long_term_edf=unname(summary(fit)$s.table[1,"edf"]),residual_dependence_flag=if(anyNA(diagnostics$residual_dependence_flag))NA else any(diagnostics$residual_dependence_flag),model_warning=if(length(warnings))paste(unique(warnings),collapse="; ") else NA_character_),curves=dplyr::bind_rows(curves,dplyr::bind_rows(pooled)),diagnostics=diagnostics)
}

pool_screen <- function(data,denominator,count,proportion=FALSE,coverage_min=0.5,linkage_min=0.9) {
  common <- common_pooling_months(data,denominator,count,proportion,coverage_min,linkage_min)
  site_names <- sort(unique(data$site_name))
  effects <- dplyr::bind_rows(lapply(site_names,function(site){
    x <- dplyr::filter(common,site_name==site)
    fit <- fit_temporal_model(x,denominator,count,proportion)$summary
    fit %>% dplyr::mutate(site_name=site,n_common_months=dplyr::n_distinct(common$calendar_month),common_start=if(nrow(common))min(common$calendar_month) else as.Date(NA),common_end=if(nrow(common))max(common$calendar_month) else as.Date(NA))
  }))
  eligible <- effects %>% dplyr::filter(model_status=="estimated") %>% dplyr::pull(site_name)
  joint_data <- common %>% dplyr::filter(site_name %in% eligible)
  joint <- fit_joint_pool(joint_data,denominator,count,proportion)
  window <- tibble::tibble(n_registered_screen_sites=length(site_names),n_common_months=dplyr::n_distinct(common$calendar_month),common_start=if(nrow(common))min(common$calendar_month) else as.Date(NA),common_end=if(nrow(common))max(common$calendar_month) else as.Date(NA))
  list(effects=effects,meta=dplyr::bind_cols(window,meta_pool(effects)),joint=dplyr::bind_cols(window,joint$summary),curves=joint$curves,diagnostics=joint$diagnostics)
}

run_pooling <- function(inputs,coverage_min=0.5,linkage_min=0.9) {
  outputs <- list(effects=list(),meta=list(),joint=list(),curves=list(),diagnostics=list())
  add <- function(data,keys,denominator,count,proportion=FALSE) {
    result <- pool_screen(data,denominator,count,proportion,coverage_min,linkage_min)
    for(name in names(outputs)) {
      x <- result[[name]]
      if(nrow(x)) outputs[[name]][[length(outputs[[name]])+1]] <<- dplyr::bind_cols(keys[rep(1,nrow(x)),,drop=FALSE],x)
    }
  }
  for(organism in unique(inputs$culture$organism_category)) for(den in c("n_icu_days","n_icu_admissions")) {
    add(dplyr::filter(inputs$culture,organism_category==organism),tibble::tibble(screen="organism",organism_category=organism,antimicrobial_category=NA_character_,specimen_stratum="Overall",outcome="organism_detection_rate",denominator=den,effect_scale="annualized endpoint rate ratio"),den,"n_detection_events")
  }
  if(nrow(inputs$ast)) {
    keys <- inputs$ast %>% dplyr::distinct(organism_category,antimicrobial_category,specimen_stratum)
    for(i in seq_len(nrow(keys))) {
      dat <- dplyr::semi_join(inputs$ast,keys[i,],by=names(keys))
      for(outcome in c("susceptible_detection_rate","non_susceptible_detection_rate","non_susceptible_fraction")) {
        fraction <- outcome=="non_susceptible_fraction";den<-if(fraction)"n_interpretable" else "n_icu_days";count<-if(outcome=="susceptible_detection_rate")"n_susceptible" else "n_non_susceptible"
        add(dat,keys[i,] %>% dplyr::mutate(screen="susceptibility",outcome=outcome,denominator=den,effect_scale=if(fraction)"annualized endpoint odds ratio" else "annualized endpoint rate ratio"),den,count,fraction)
      }
    }
  }
  outputs <- lapply(outputs,dplyr::bind_rows)
  if(nrow(outputs$meta)) outputs$meta <- outputs$meta %>% dplyr::group_by(screen,outcome,denominator,specimen_stratum) %>% dplyr::mutate(fdr_p_value=p.adjust(p_value,"BH"),direction=dplyr::case_when(model_status!="estimated"~"Not estimated",inference_status!="eligible"~"Exploratory",fdr_p_value<0.05 & annual_effect>1~"Increasing",fdr_p_value<0.05 & annual_effect<1~"Decreasing",TRUE~"No clear net change")) %>% dplyr::ungroup()
  outputs
}

plot_pooling <- function(result,out_dir) {
  keys <- c("screen","organism_category","antimicrobial_category","specimen_stratum","outcome","denominator")
  screens <- dplyr::distinct(result$meta,dplyr::across(dplyr::all_of(keys)))
  for(i in seq_len(nrow(screens))) {
    key <- screens[i,]; filename<-sprintf("screen_%04d",i)
    title <- gsub("_"," ",paste(na.omit(unlist(key[c("organism_category","antimicrobial_category","specimen_stratum","outcome","denominator")])),collapse=" | "))
    effects <- dplyr::semi_join(result$effects,key,by=keys) %>% dplyr::filter(model_status=="estimated") %>% dplyr::transmute(label=site_name,estimate=annual_irr,low=annual_ci_low,high=annual_ci_high)
    meta <- dplyr::semi_join(result$meta,key,by=keys) %>% dplyr::filter(model_status=="estimated")
    if(nrow(meta))effects<-dplyr::bind_rows(effects,tibble::tibble(label="Random-effects pooled",estimate=meta$annual_effect,low=meta$ci_low,high=meta$ci_high))
    if(nrow(effects)) {
      p <- ggplot2::ggplot(effects,ggplot2::aes(estimate,factor(label,levels=rev(effects$label))))+ggplot2::geom_vline(xintercept=1,linetype=2,color="grey60")+ggplot2::geom_segment(ggplot2::aes(x=low,xend=high,yend=factor(label,levels=rev(effects$label))))+ggplot2::geom_point(size=2)+ggplot2::scale_x_log10()+ggplot2::labs(title=title,x=if(key$outcome=="non_susceptible_fraction")"Annualized endpoint odds ratio (95% CI)" else "Annualized endpoint rate ratio (95% CI)",y=NULL,caption=NULL)+ggplot2::theme_bw(base_size=10)
      ggplot2::ggsave(file.path(out_dir,paste0(filename,"_forest.png")),p,width=11,height=max(4,1+nrow(effects)*0.4),dpi=160)
    }
    if(nrow(result$curves)) {
      curves<-dplyr::semi_join(result$curves,key,by=keys)
      if(nrow(curves)) {
        site_names<-sort(unique(curves$site_name[curves$curve_type=="site"]))
        palette<-c("Equal-site mean"="black",setNames(grDevices::hcl.colors(length(site_names),"Dark 3"),site_names))
        p<-ggplot2::ggplot(curves,ggplot2::aes(calendar_month,fitted_value,color=site_name,group=site_name))+ggplot2::geom_line(alpha=0.7)+ggplot2::geom_ribbon(data=dplyr::filter(curves,curve_type=="pooled_equal_site"),ggplot2::aes(ymin=ci_low,ymax=ci_high),fill="grey50",color=NA,alpha=0.2)+ggplot2::geom_line(data=dplyr::filter(curves,curve_type=="pooled_equal_site"),color="black",linewidth=1)+ggplot2::scale_color_manual(values=palette)+ggplot2::labs(title=title,x=NULL,y=if(key$outcome=="non_susceptible_fraction")"Non-susceptible tested fraction" else if(key$denominator=="n_icu_days")"Detections per 100 ICU days" else "Detections per 100 ICU admissions",color="Site",caption=NULL)+ggplot2::theme_bw(base_size=10)+ggplot2::theme(legend.position="bottom")
        ggplot2::ggsave(file.path(out_dir,paste0(filename,"_joint.png")),p,width=11,height=6,dpi=160)
      }
    }
  }
  readr::write_csv(dplyr::mutate(screens,figure_prefix=sprintf("screen_%04d",seq_len(nrow(screens)))),file.path(out_dir,"figure_index.csv"))
}
