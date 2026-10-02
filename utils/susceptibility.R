# Canonical CLIF categories are authoritative; raw MICs/text never infer status.
normalize_susceptibility <- function(ast, mcide_dir = "config/mcide") {
  required <- c("organism_id", "antimicrobial_category", "susceptibility_category")
  missing <- setdiff(required, names(ast))
  if (length(missing)) stop("microbiology_susceptibility missing mCIDE fields: ", paste(missing, collapse = ", "))
  drugs <- readr::read_csv(file.path(mcide_dir, "clif_microbiology_susceptibility_antibiotics_category.csv"), show_col_types = FALSE)$antimicrobial_category
  ast %>% dplyr::mutate(organism_id = dplyr::na_if(as.character(organism_id), ""),
    antimicrobial_category_raw = clean_micro_label(antimicrobial_category),
    antimicrobial_category = stringr::str_replace_all(antimicrobial_category_raw, "[ -]+", "_"),
    antimicrobial_mapped = antimicrobial_category %in% drugs,
    susceptibility_category_raw = clean_micro_label(susceptibility_category),
    susceptibility_category = stringr::str_replace_all(susceptibility_category_raw, "[ -]+", "_"),
    unmapped_result = !is.na(susceptibility_category) & !susceptibility_category %in% c("susceptible", "non_susceptible", "indeterminate", "na"),
    susceptibility_category = dplyr::case_when(susceptibility_category %in% c("susceptible", "non_susceptible", "indeterminate") ~ susceptibility_category, TRUE ~ "unavailable"))
}

# Shared by local and pooled models. Culture capture is assumed complete.
# Zero organisms is distinct from untested organisms; source activity is still required.
screen_susceptibility_months <- function(x, coverage_min = 0.5, linkage_min = 0.9) {
  if (any(!is.finite(c(coverage_min,linkage_min))) || any(c(coverage_min,linkage_min)<0 | c(coverage_min,linkage_min)>1)) stop("AST coverage/linkage thresholds must be between 0 and 1.")
  required <- c("n_culture_isolates","n_linkable_culture_isolates","n_interpretable","n_observed_culture_events","n_positive_rows_missing_organism_category","n_icu_days")
  if (!all(required %in% names(x))) stop("Susceptibility export lacks observation/linkage fields; rerun the site pipeline.")
  if (any(!is.finite(x$n_culture_isolates) | !is.finite(x$n_linkable_culture_isolates) | !is.finite(x$n_interpretable) | x$n_culture_isolates<0 | x$n_linkable_culture_isolates<0 | x$n_interpretable<0 | x$n_linkable_culture_isolates>x$n_culture_isolates | x$n_interpretable>x$n_linkable_culture_isolates)) stop("Invalid AST isolate/test counts.")
  x %>% dplyr::mutate(
    linkage_fraction = dplyr::if_else(n_culture_isolates>0,n_linkable_culture_isolates/n_culture_isolates,NA_real_),
    testing_fraction = dplyr::if_else(n_culture_isolates>0,n_interpretable/n_culture_isolates,NA_real_),
    testing_fraction_linkable = dplyr::if_else(n_linkable_culture_isolates>0,n_interpretable/n_linkable_culture_isolates,NA_real_),
    source_observed = is.finite(n_observed_culture_events) & n_observed_culture_events>0,
    zero_detection_validated = n_culture_isolates==0 & source_observed & n_positive_rows_missing_organism_category==0,
    linkage_eligible = n_culture_isolates>0 & dplyr::coalesce(linkage_fraction>=linkage_min,FALSE),
    observation_status = dplyr::case_when(!is.finite(n_icu_days) | n_icu_days<=0 ~ "no_icu_exposure", !source_observed ~ "source_unavailable",
      n_culture_isolates==0 & n_positive_rows_missing_organism_category>0 ~ "organism_identification_incomplete",
      zero_detection_validated ~ "organism_not_detected",
      !linkage_eligible ~ "insufficient_linkage", n_interpretable==0 ~ "no_interpretable_tests", testing_fraction<coverage_min ~ "insufficient_testing_coverage", TRUE ~ "adequate_testing"),
    rate_model_eligible = observation_status %in% c("organism_not_detected","adequate_testing"),
    fraction_model_eligible = source_observed & linkage_eligible & n_interpretable>0,
    minimum_testing_fraction_for_rate=coverage_min, minimum_linkage_fraction=linkage_min)
}

build_susceptibility_analysis <- function(rows, ast, denominators, mcide_dir = "config/mcide", coverage_min = 0.5, linkage_min = 0.9) {
  ast <- normalize_susceptibility(ast, mcide_dir)
  qc <- tibble::tibble(metric=c("ast_rows","missing_organism_id_rows","unmapped_antimicrobial_rows","unmapped_result_rows"),n=c(nrow(ast),sum(is.na(ast$organism_id)),sum(!ast$antimicrobial_mapped),sum(ast$unmapped_result)))
  # Retain every positive isolate in coverage denominators, including missing linkage IDs.
  positives <- rows %>% dplyr::filter(positive_culture) %>% dplyr::mutate(organism_id=dplyr::na_if(as.character(organism_id),""),organism_category=dplyr::na_if(organism_category,""))
  isolates <- positives %>% dplyr::filter(!is.na(organism_category)) %>% dplyr::distinct(organism_id,culture_event_id,patient_id,hospitalization_id,organism_category,fluid_category,collect_dttm,dplyr::across(dplyr::any_of(c("organism_name","organism_group")))) %>% dplyr::mutate(calendar_month=lubridate::floor_date(collect_dttm,"month"))
  linkable <- isolates %>% dplyr::filter(!is.na(organism_id))
  if (anyDuplicated(linkable$organism_id)) stop("organism_id maps to multiple ICU isolate events or categories; fix linkage before susceptibility analysis.")
  qc <- dplyr::bind_rows(qc,tibble::tibble(metric=c("positive_culture_rows_missing_organism_id","positive_culture_rows_missing_organism_category","linkable_icu_isolates","ast_rows_unlinked_to_icu_isolate"),n=c(sum(is.na(positives$organism_id)),sum(is.na(positives$organism_category)),nrow(linkable),sum(!ast$organism_id %in% linkable$organism_id))))
  tests <- ast %>% dplyr::filter(!is.na(organism_id),antimicrobial_mapped) %>% dplyr::group_by(organism_id,antimicrobial_category) %>%
    dplyr::summarise(n_source_rows=dplyr::n(),conflicting_result=all(c("susceptible","non_susceptible") %in% susceptibility_category),
      susceptibility_status=if(conflicting_result)"indeterminate" else if(any(susceptibility_category=="non_susceptible"))"non_susceptible" else if(any(susceptibility_category=="susceptible"))"susceptible" else if(any(susceptibility_category=="indeterminate"))"indeterminate" else "unavailable",.groups="drop")
  linked <- dplyr::inner_join(linkable,tests,by="organism_id",relationship="one-to-many")
  qc <- dplyr::bind_rows(qc,tibble::tibble(metric=c("conflicting_isolate_antimicrobial_pairs","linked_icu_isolate_antimicrobial_pairs","interpretable_icu_isolate_antimicrobial_pairs"),n=c(sum(tests$conflicting_result),nrow(linked),sum(linked$susceptibility_status %in% c("susceptible","non_susceptible")))))
  strata <- function(x) dplyr::bind_rows(dplyr::mutate(x,specimen_stratum="Overall"),dplyr::mutate(x,specimen_stratum=dplyr::coalesce(fluid_category,"missing")))
  totals <- strata(isolates) %>% dplyr::group_by(calendar_month,organism_category,specimen_stratum) %>% dplyr::summarise(n_culture_isolates=dplyr::n(),n_linkable_culture_isolates=sum(!is.na(organism_id)),n_missing_organism_id=sum(is.na(organism_id)),.groups="drop")
  source_counts <- strata(rows %>% dplyr::mutate(calendar_month=lubridate::floor_date(collect_dttm,"month"))) %>% dplyr::group_by(calendar_month,specimen_stratum) %>% dplyr::summarise(n_observed_culture_events=dplyr::n_distinct(culture_event_id),n_positive_rows_missing_organism_category=sum(positive_culture & (is.na(organism_category) | organism_category=="")),.groups="drop")
  linkage_qc <- totals %>% tidyr::complete(calendar_month=denominators$calendar_month,tidyr::nesting(organism_category,specimen_stratum),fill=list(n_culture_isolates=0L,n_linkable_culture_isolates=0L,n_missing_organism_id=0L)) %>% dplyr::mutate(linkage_fraction=dplyr::if_else(n_culture_isolates>0,n_linkable_culture_isolates/n_culture_isolates,NA_real_),linkage_eligible=linkage_fraction>=linkage_min)
  status <- if(!nrow(linkable))"no_linkable_icu_isolates" else if(!nrow(linked))"no_linked_icu_tests" else if(!any(linked$susceptibility_status %in% c("susceptible","non_susceptible")))"no_interpretable_tests" else "ready_for_models"
  if(!nrow(linked)) return(list(monthly=tibble::tibble(),qc=qc,linkage_qc=linkage_qc,analysis_status=status))
  monthly <- strata(linked) %>% dplyr::group_by(calendar_month,organism_category,antimicrobial_category,specimen_stratum) %>%
    dplyr::summarise(n_susceptible=sum(susceptibility_status=="susceptible"),n_non_susceptible=sum(susceptibility_status=="non_susceptible"),n_indeterminate=sum(susceptibility_status=="indeterminate"),n_unavailable_reported=sum(susceptibility_status=="unavailable"),n_conflicting=sum(conflicting_result),.groups="drop") %>%
    tidyr::complete(calendar_month=denominators$calendar_month,tidyr::nesting(organism_category,antimicrobial_category,specimen_stratum),fill=list(n_susceptible=0L,n_non_susceptible=0L,n_indeterminate=0L,n_unavailable_reported=0L,n_conflicting=0L)) %>%
    dplyr::left_join(totals,by=c("calendar_month","organism_category","specimen_stratum")) %>% dplyr::left_join(source_counts,by=c("calendar_month","specimen_stratum")) %>% dplyr::left_join(dplyr::select(denominators,-dplyr::any_of("n_observed_culture_events")),by="calendar_month") %>%
    dplyr::mutate(dplyr::across(c(n_culture_isolates,n_linkable_culture_isolates,n_missing_organism_id,n_observed_culture_events,n_positive_rows_missing_organism_category),~dplyr::coalesce(.x,0L)),
      n_interpretable=n_susceptible+n_non_susceptible,n_without_reported_test=n_culture_isolates-n_interpretable-n_indeterminate-n_unavailable_reported,
      non_susceptible_fraction=dplyr::if_else(n_interpretable>0,n_non_susceptible/n_interpretable,NA_real_)) %>% screen_susceptibility_months(coverage_min,linkage_min) %>%
    dplyr::mutate(susceptible_per_100_icu_days=dplyr::if_else(n_icu_days>0 & (n_interpretable>0 | zero_detection_validated),100*n_susceptible/n_icu_days,NA_real_),non_susceptible_per_100_icu_days=dplyr::if_else(n_icu_days>0 & (n_interpretable>0 | zero_detection_validated),100*n_non_susceptible/n_icu_days,NA_real_))
  list(monthly=monthly,qc=qc,linkage_qc=linkage_qc,analysis_status=status)
}
