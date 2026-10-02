# CLIF 2.1 mCIDE fields are authoritative; never reinterpret MICs or unavailable tests.
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
build_susceptibility_analysis <- function(rows, ast, denominators, mcide_dir = "config/mcide") {
  ast <- normalize_susceptibility(ast, mcide_dir)
  qc <- tibble::tibble(metric = c("ast_rows", "missing_organism_id_rows", "unmapped_antimicrobial_rows", "unmapped_result_rows"), n = c(nrow(ast), sum(is.na(ast$organism_id)), sum(!ast$antimicrobial_mapped), sum(ast$unmapped_result)))
  culture_isolates <- rows %>% dplyr::filter(positive_culture, !is.na(organism_id), organism_id != "", !is.na(organism_category)) %>%
    dplyr::distinct(organism_id, culture_event_id, patient_id, hospitalization_id, organism_category, fluid_category, collect_dttm) %>%
    dplyr::mutate(calendar_month = lubridate::floor_date(collect_dttm, "month"))
  if (anyDuplicated(culture_isolates$organism_id)) stop("organism_id maps to multiple ICU isolate events or categories; fix linkage before susceptibility analysis.")
  qc <- dplyr::bind_rows(qc, tibble::tibble(metric = c("positive_culture_rows_missing_organism_id", "ast_rows_unlinked_to_icu_isolate"), n = c(sum(rows$positive_culture & (is.na(rows$organism_id) | rows$organism_id == "")), sum(!ast$organism_id %in% culture_isolates$organism_id))))
  # No AST result timestamps are specified by CLIF. Conflicting S/NS reports are indeterminate,
  # rather than arbitrarily taking a row as the final result or assuming resistance.
  tests <- ast %>% dplyr::filter(!is.na(organism_id), antimicrobial_mapped) %>%
    dplyr::group_by(organism_id, antimicrobial_category) %>%
    dplyr::summarise(n_source_rows = dplyr::n(), conflicting_result = all(c("susceptible", "non_susceptible") %in% susceptibility_category),
      susceptibility_status = if (conflicting_result) "indeterminate" else if (any(susceptibility_category == "non_susceptible")) "non_susceptible" else if (any(susceptibility_category == "susceptible")) "susceptible" else if (any(susceptibility_category == "indeterminate")) "indeterminate" else "unavailable", .groups = "drop")
  qc <- dplyr::bind_rows(qc, tibble::tibble(metric = "conflicting_isolate_antimicrobial_pairs", n = sum(tests$conflicting_result)))
  linked <- dplyr::inner_join(culture_isolates, tests, by = "organism_id", relationship = "one-to-many")
  if (!nrow(linked)) return(list(monthly = tibble::tibble(), qc = qc))
  # Overall and specimen-specific analyses are separate strata.
  linked <- dplyr::bind_rows(dplyr::mutate(linked, specimen_stratum = "Overall"), dplyr::mutate(linked, specimen_stratum = dplyr::coalesce(fluid_category, "missing")))
  all_isolates <- dplyr::bind_rows(dplyr::mutate(culture_isolates, specimen_stratum = "Overall"), dplyr::mutate(culture_isolates, specimen_stratum = dplyr::coalesce(fluid_category, "missing"))) %>%
    dplyr::count(calendar_month, organism_category, specimen_stratum, name = "n_culture_isolates")
  monthly <- linked %>% dplyr::group_by(calendar_month, organism_category, antimicrobial_category, specimen_stratum) %>%
    dplyr::summarise(n_susceptible = sum(susceptibility_status == "susceptible"), n_non_susceptible = sum(susceptibility_status == "non_susceptible"), n_indeterminate = sum(susceptibility_status == "indeterminate"), n_unavailable_reported = sum(susceptibility_status == "unavailable"), n_conflicting = sum(conflicting_result), .groups = "drop") %>%
    tidyr::complete(calendar_month = denominators$calendar_month, tidyr::nesting(organism_category, antimicrobial_category, specimen_stratum), fill = list(n_susceptible = 0L, n_non_susceptible = 0L, n_indeterminate = 0L, n_unavailable_reported = 0L, n_conflicting = 0L)) %>%
    dplyr::left_join(all_isolates, by = c("calendar_month", "organism_category", "specimen_stratum")) %>% dplyr::left_join(denominators, by = "calendar_month") %>%
    dplyr::mutate(n_culture_isolates = dplyr::coalesce(n_culture_isolates, 0L), n_interpretable = n_susceptible + n_non_susceptible,
      n_without_reported_test = pmax(0L, n_culture_isolates - n_interpretable - n_indeterminate - n_unavailable_reported),
      testing_fraction = dplyr::if_else(n_culture_isolates > 0, n_interpretable / n_culture_isolates, NA_real_),
      non_susceptible_fraction = dplyr::if_else(n_interpretable > 0, n_non_susceptible / n_interpretable, NA_real_),
      susceptible_per_100_icu_days = dplyr::if_else(n_icu_days > 0 & n_interpretable > 0, 100 * n_susceptible / n_icu_days, NA_real_),
      non_susceptible_per_100_icu_days = dplyr::if_else(n_icu_days > 0 & n_interpretable > 0, 100 * n_non_susceptible / n_icu_days, NA_real_))
  list(monthly = monthly, qc = qc)
}
