# Descriptive Table 1: one observation per classified hospitalization with study ICU time.
build_ase_characteristics <- function(data, classification, patient, clinical=NULL) {
  patient_fields <- c("sex_category", "race_category", "ethnicity_category")
  p <- patient %>% select(any_of(c("patient_id", patient_fields))) %>% distinct()
  if (anyDuplicated(p$patient_id)) stop("Conflicting patient demographics for the same patient_id")
  h <- data$hospitalization %>% semi_join(data$icu_admissions, by="hospitalization_id") %>%
    inner_join(select(classification, hospitalization_id, ase_group), by="hospitalization_id") %>%
    left_join(p, by="patient_id")
  if (anyDuplicated(h$hospitalization_id)) stop("Characteristics require one row per hospitalization")
  if (!is.null(clinical)) {
    if (anyDuplicated(clinical$hospitalization_id)) stop("Clinical characteristics must be unique by hospitalization")
    h <- left_join(h,clinical,by="hospitalization_id")
  }
  available <- names(h)
  for (field in c(patient_fields, "admission_type_category", "discharge_category")) {
    if (!field %in% names(h)) h[[field]] <- NA_character_
    h[[field]] <- trimws(tolower(as.character(h[[field]])))
    h[[field]][is.na(h[[field]]) | h[[field]] %in% c("", "unknown", "missing", "na", "n/a")] <- NA_character_
  }
  for (field in c("recorded_imv","recorded_vasopressor")) {
    if (!field %in% names(h)) h[[field]] <- NA_character_
  }
  h$in_hospital_death <- ifelse(is.na(h$discharge_category),NA_character_,ifelse(h$discharge_category=="expired","yes","no"))
  stays <- data$icu_admissions %>% group_by(hospitalization_id) %>%
    summarise(icu_days=sum(icu_los_days), icu_stays=n(), .groups="drop")
  cultures <- data$events %>% group_by(hospitalization_id) %>%
    summarise(culture_events=n_distinct(culture_event_id), positive_events=n_distinct(culture_event_id[any_positive_culture]), .groups="drop")
  h <- h %>% left_join(stays, by="hospitalization_id") %>% left_join(cultures, by="hospitalization_id") %>%
    mutate(culture_events=coalesce(culture_events,0L), positive_events=coalesce(positive_events,0L),
           hospital_days=as.numeric(difftime(safe_ts(discharge_dttm),safe_ts(admission_dttm),units="days")),
           age_at_admission=suppressWarnings(as.numeric(age_at_admission)),
           any_icu_culture=if_else(culture_events>0,"yes","no"),
           any_positive_icu_culture=if_else(positive_events>0,"yes","no"))
  continuous <- c(age_at_admission="Age at admission (years)", hospital_days="Hospital length of stay (days)",
    icu_days="ICU days within study window per hospitalization", icu_stays="ICU stays overlapping study window per hospitalization",
    culture_events="ICU culture events per hospitalization", positive_events="Positive ICU culture events per hospitalization")
  categorical <- c(sex_category="Sex", race_category="Race", ethnicity_category="Ethnicity",
    admission_type_category="Admission type", discharge_category="Discharge disposition",
    in_hospital_death="In-hospital death (expired discharge)",
    recorded_imv="Recorded invasive mechanical ventilation during hospitalization",
    recorded_vasopressor="Recorded nonprocedural vasopressor during hospitalization",
    any_icu_culture="Any ICU culture in study window", any_positive_icu_culture="Any positive ICU culture in study window")
  rows <- list()
  for (g in ase_group_levels) {
    x <- filter(h, ase_group==g); total <- nrow(x)
    add <- function(field,label,level,type,n,observed,missing,percent=NA_real_,median=NA_real_,p25=NA_real_,p75=NA_real_,display=NA_character_) {
      rows[[length(rows)+1L]] <<- tibble(ase_group=g, characteristic=field, label=label, level=level,
        summary_type=type, n_hospitalizations=total, n_observed=observed, n_missing=missing,
        n=n, percent=percent, median=median, p25=p25, p75=p75, display=display,
        availability=if(total==0)"empty_subgroup" else if(field %in% c(names(continuous),"any_icu_culture","any_positive_icu_culture","hospitalizations","unique_patients","in_hospital_death") || field %in% available) "available" else "source_column_unavailable")
    }
    for (field in names(continuous)) {
      v <- x[[field]]; v[!is.finite(v) | v<0] <- NA_real_
      observed <- sum(!is.na(v)); q <- if(observed) unname(quantile(v,c(.25,.5,.75),na.rm=TRUE)) else rep(NA_real_,3)
      add(field,continuous[[field]],"", "median_iqr",NA_integer_,observed,total-observed,
        median=q[2],p25=q[1],p75=q[3],display=if(observed)sprintf("%.1f [%.1f, %.1f]",q[2],q[1],q[3]) else "Unavailable")
    }
    for (field in names(categorical)) {
      v <- x[[field]]; observed <- sum(!is.na(v)); levels <- sort(unique(h[[field]][!is.na(h[[field]])]))
      for (level in c(levels,"Missing / unknown")) {
        n <- if(level=="Missing / unknown")sum(is.na(v)) else sum(v==level,na.rm=TRUE)
        percent <- if(total)100*n/total else NA_real_
        add(field,categorical[[field]],level,"n_percent",n,observed,total-observed,percent,
          display=if(total)sprintf("%d (%.1f%%)",n,percent) else "Unavailable")
      }
    }
    add("hospitalizations","Hospitalizations","","count",total,total,0,display=as.character(total))
    add("unique_patients","Unique patients within subgroup","","count",n_distinct(x$patient_id),total,0,display=as.character(n_distinct(x$patient_id)))
  }
  long <- bind_rows(rows)
  wide <- long %>% select(characteristic,label,level,summary_type,ase_group,display,n_observed,n_missing,availability) %>%
    tidyr::pivot_wider(names_from=ase_group,values_from=c(display,n_observed,n_missing,availability))
  overlap <- h %>% distinct(patient_id,ase_group) %>% count(patient_id) %>% summarise(n=sum(n>1))
  qc <- tibble(metric=c("classified_hospitalizations","unique_patients_overall","patients_in_multiple_subgroups"),
    n=c(nrow(h),n_distinct(h$patient_id),overlap$n))
  list(long=long,wide=wide,qc=qc)
}
