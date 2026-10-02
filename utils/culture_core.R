# Shared definitions: half-open ICU stays, culture events, and result categories.
# Functions are pure where possible so synthetic fixtures can exercise site-independent logic.
read_clif_csv <- function(path, ...) {
  header <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  types <- readr::cols(.default = readr::col_guess())
  for (id in intersect(header, c("patient_id", "hospitalization_id", "organism_id"))) types$cols[[id]] <- readr::col_character()
  readr::read_csv(path, col_types = types, ...)
}
safe_ts <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = tz))
  if (is.numeric(x)) return(as.POSIXct(ifelse(x > 1e12, x / 1000, x), origin = "1970-01-01", tz = tz))
  suppressWarnings(lubridate::parse_date_time(x, orders = c("ymd_HMS", "ymd_HM", "ymd", "ymdTz", "ymdT", "mdy_HMS", "mdy_HM", "mdy"), tz = tz, quiet = TRUE))
}
clean_micro_label <- function(x) {
  x <- stringr::str_to_lower(stringr::str_squish(as.character(x)))
  dplyr::na_if(x, "")
}
culture_result_levels <- c("Negative/no growth", "Positive", "Mixed/contaminated", "Indeterminate")
classify_culture_result <- function(organism_group, organism_category, organism_name) {
  fields <- lapply(list(organism_group, organism_category, organism_name), clean_micro_label)
  label <- dplyr::coalesce(fields[[2]], fields[[1]], fields[[3]])
  text <- paste(dplyr::coalesce(fields[[1]], ""), dplyr::coalesce(fields[[2]], ""), dplyr::coalesce(fields[[3]], ""))
  negative <- stringr::str_detect(text, "no[_ ]growth|negative for|\\b(no|none|not) .*isolated|no .*detected")
  mixed <- stringr::str_detect(text, "mixed[_ ]([a-z]+[_ ])*flora|normal[_ ]([a-z]+[_ ])*flora|commensal[_ ]flora|contaminat")
  unknown <- is.na(label) | label %in% c("na", "unknown", "missing", "other", "other_unspecified", "indeterminate", "pending", "not_reported", "not_available") |
    stringr::str_detect(text, "pending|cancelled|canceled|not performed|insufficient|indeterminate")
  dplyr::case_when(mixed ~ "Mixed/contaminated", negative ~ "Negative/no growth", unknown ~ "Indeterminate", TRUE ~ "Positive")
}
collapse_result_status <- function(x) {
  # A named isolate takes precedence; unresolved reports cannot establish no growth.
  if (any(x == "Positive")) "Positive" else if (any(x == "Mixed/contaminated")) "Mixed/contaminated" else if (any(x == "Indeterminate")) "Indeterminate" else "Negative/no growth"
}
merge_icu_stays <- function(adt, hospitalization) {
  h <- hospitalization %>% dplyr::mutate(dplyr::across(c(admission_dttm, discharge_dttm), safe_ts))
  if (anyDuplicated(h$hospitalization_id)) stop("hospitalization_id is not unique in hospitalization.")
  adt %>% dplyr::transmute(hospitalization_id, icu_in_dttm = safe_ts(in_dttm), icu_out_dttm_raw = safe_ts(out_dttm), location_category = clean_micro_label(location_category)) %>%
    dplyr::filter(location_category == "icu", !is.na(icu_in_dttm)) %>%
    dplyr::left_join(h %>% dplyr::select(patient_id, hospitalization_id, admission_dttm, discharge_dttm), by = "hospitalization_id") %>%
    dplyr::mutate(icu_interval_missing_out = is.na(icu_out_dttm_raw), icu_out_dttm = dplyr::coalesce(icu_out_dttm_raw, discharge_dttm)) %>%
    dplyr::filter(!is.na(patient_id), !is.na(icu_out_dttm), icu_out_dttm > icu_in_dttm) %>%
    dplyr::arrange(patient_id, hospitalization_id, icu_in_dttm, icu_out_dttm) %>%
    dplyr::group_by(patient_id, hospitalization_id) %>%
    dplyr::mutate(prior_max = dplyr::lag(cummax(as.numeric(icu_out_dttm))), seq = cumsum(is.na(prior_max) | as.numeric(icu_in_dttm) > prior_max)) %>%
    dplyr::group_by(patient_id, hospitalization_id, seq) %>%
    dplyr::summarise(admission_dttm = dplyr::first(admission_dttm), discharge_dttm = dplyr::first(discharge_dttm), icu_in_dttm = min(icu_in_dttm), icu_out_dttm = max(icu_out_dttm), icu_interval_missing_out = any(icu_interval_missing_out), n_icu_adt_rows = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(patient_id, hospitalization_id, icu_in_dttm) %>%
    dplyr::mutate(icu_admission_id = dplyr::row_number(), icu_interval_id = icu_admission_id, icu_admission_month = lubridate::floor_date(icu_in_dttm, "month")) %>% dplyr::select(-seq)
}
apply_specimen_overrides <- function(rows, overrides) {
  rows$fluid_category_original <- rows$fluid_category
  rows$specimen_mapping_overridden <- FALSE
  if (!nrow(overrides)) return(rows)
  for (i in seq_len(nrow(overrides))) {
    rule <- overrides[i, ]
    hit <- !is.na(rows$fluid_name) & rows$fluid_name == clean_micro_label(rule$fluid_name) &
      !is.na(rows$fluid_category_original) & rows$fluid_category_original == clean_micro_label(rule$original_category) &
      as.Date(rows$collect_dttm) >= as.Date(rule$start_date) & as.Date(rows$collect_dttm) <= as.Date(rule$end_date)
    if (any(hit & rows$specimen_mapping_overridden)) stop("Overlapping specimen override rules.")
    if (is.na(rule$reason) || !nzchar(rule$reason)) stop("Specimen overrides require a documented reason.")
    rows$fluid_category[hit] <- clean_micro_label(rule$replacement_category)
    rows$specimen_mapping_overridden[hit] <- TRUE
  }
  rows
}
build_culture_data <- function(hospitalization, adt, micro, start = as.POSIXct(NA), end = as.POSIXct(NA), overrides = NULL) {
  stays <- merge_icu_stays(adt, hospitalization)
  stays$icu_in_dttm_clipped <- if (!is.na(start)) pmax(stays$icu_in_dttm, start) else stays$icu_in_dttm
  # end is exclusive, including when a calendar study_end_date is supplied.
  stays$icu_out_dttm_clipped <- if (!is.na(end)) pmin(stays$icu_out_dttm, end) else stays$icu_out_dttm
  stays <- stays %>% dplyr::filter(icu_out_dttm_clipped > icu_in_dttm_clipped) %>%
    dplyr::mutate(icu_los_days = as.numeric(difftime(icu_out_dttm_clipped, icu_in_dttm_clipped, units = "days")))
  for (col in c("organism_id", "organism_name", "organism_category", "organism_group", "order_dttm", "result_dttm", "fluid_name", "fluid_category", "method_name")) if (!col %in% names(micro)) micro[[col]] <- NA_character_
  rows <- micro %>% dplyr::mutate(microbiology_row_id = dplyr::row_number(), dplyr::across(c(order_dttm, collect_dttm, result_dttm), safe_ts),
    dplyr::across(c(fluid_name, fluid_category, method_name, method_category, organism_name, organism_category, organism_group), clean_micro_label), organism_id = as.character(organism_id)) %>%
    dplyr::filter(method_category == "culture", !is.na(collect_dttm))
  if (!is.na(start)) rows <- dplyr::filter(rows, collect_dttm >= start)
  if (!is.na(end)) rows <- dplyr::filter(rows, collect_dttm < end)
  if (!is.null(overrides)) rows <- apply_specimen_overrides(rows, overrides)
  rows <- rows %>% dplyr::inner_join(stays %>% dplyr::select(patient_id, hospitalization_id, icu_admission_id, icu_interval_id, icu_in_dttm, icu_out_dttm, icu_interval_missing_out, icu_admission_month), by = c("patient_id", "hospitalization_id"), relationship = "many-to-many") %>%
    dplyr::filter(collect_dttm >= icu_in_dttm, collect_dttm < icu_out_dttm) %>%
    dplyr::group_by(patient_id, hospitalization_id, icu_admission_id, order_dttm, collect_dttm, fluid_name, method_name) %>%
    dplyr::mutate(culture_event_id = dplyr::cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::mutate(isolate_key = dplyr::coalesce(dplyr::na_if(organism_id, ""), paste(organism_category, organism_group, organism_name, sep = "|"))) %>%
    dplyr::arrange(!is.na(result_dttm), result_dttm) %>%
    dplyr::group_by(culture_event_id, isolate_key) %>% dplyr::slice_tail(n = 1) %>% dplyr::ungroup() %>%
    dplyr::mutate(result_status = classify_culture_result(organism_group, organism_category, organism_name), positive_culture = result_status == "Positive", no_growth = result_status == "Negative/no growth") %>%
    dplyr::select(-isolate_key)
  events <- rows %>% dplyr::group_by(culture_event_id, patient_id, hospitalization_id, icu_admission_id, icu_interval_id, icu_in_dttm, icu_out_dttm, icu_admission_month, order_dttm, collect_dttm, fluid_name, method_name) %>%
    dplyr::summarise(fluid_category = dplyr::first(fluid_category), method_category = "culture", n_culture_rows = dplyr::n(), result_status = collapse_result_status(result_status),
      organism_groups = paste(sort(unique(stats::na.omit(organism_group))), collapse = "; "), .groups = "drop") %>%
    dplyr::mutate(any_positive_culture = result_status == "Positive")
  list(hospitalization = hospitalization, icu_admissions = stays, rows = rows, events = events)
}
read_culture_data <- function(start = as.POSIXct(NA), end = as.POSIXct(NA)) {
  path <- project_path("config", "specimen_category_overrides.csv")
  overrides <- if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else NULL
  build_culture_data(read_tbl("hospitalization"), read_tbl("adt"), read_tbl("microbiology_culture"), start, end, overrides)
}
monthly_icu_denominators <- function(stays, month_seq) {
  admissions <- stays %>% dplyr::filter(icu_in_dttm >= icu_in_dttm_clipped) %>% dplyr::count(calendar_month = icu_admission_month, name = "n_icu_admissions")
  windows <- tibble::tibble(calendar_month = month_seq, month_end = month_seq %m+% lubridate::period(month = 1))
  days <- tidyr::crossing(stays, windows) %>% dplyr::mutate(overlap_start = pmax(icu_in_dttm_clipped, calendar_month), overlap_end = pmin(icu_out_dttm_clipped, month_end), overlap_days = as.numeric(difftime(overlap_end, overlap_start, units = "days"))) %>%
    dplyr::filter(overlap_days > 0) %>% dplyr::group_by(calendar_month) %>% dplyr::summarise(n_icu_days = sum(overlap_days), .groups = "drop")
  tibble::tibble(calendar_month = month_seq) %>% dplyr::left_join(admissions, by = "calendar_month") %>% dplyr::left_join(days, by = "calendar_month") %>%
    dplyr::mutate(n_icu_admissions = dplyr::coalesce(n_icu_admissions, 0L), n_icu_days = dplyr::coalesce(n_icu_days, 0))
}

# Shared taxonomy fallback; bacterial species names can contain viral substrings.
classify_microbe_taxonomy <- function(x) {
  x_clean <- stringr::str_to_lower(dplyr::coalesce(x, ""))
  dplyr::case_when(
    stringr::str_detect(x_clean, "candida|yeast|fung|aspergillus|cryptococcus|mold|mould|saccharomyces|fusarium|mucor|rhizopus|pneumocystis|torulopsis") ~ "Fungi/yeast",
    stringr::str_detect(x_clean, "mycobacter|tuberculosis|\\bafb\\b|acid fast") ~ "Mycobacteria/AFB",
    stringr::str_detect(x_clean, "anaerob|bacteroides|clostrid|prevotella|fusobacter|cutibacter|propionibacter|propionbacterium|leptotrichia|clostridioides") ~ "Anaerobes",
    stringr::str_detect(x_clean, "amebiasis|cryptosporidium|echinoco|giardia|protozo|toxoplasma|trichomonas|parasite") ~ "Parasite/protozoa",
    stringr::str_detect(x_clean, "haemophilus|gram[_ ]negative") ~ "Gram negative bacteria",
    stringr::str_detect(x_clean, "gram[_ ]positive") ~ "Gram positive bacteria",
    stringr::str_detect(x_clean, "adenovirus|cytomegalovirus|enterovirus|epstein|hepatitis|herpes|hhv|hiv|influenza|measles|mumps|papovavirus|parainfluenza|polyomavirus|respiratory_syncytial|rsv|rhinovirus|rotavirus|rubella|virus|viral|covid|sars|cmv") ~ "Virus",
    stringr::str_detect(x_clean, "borrelia|chlamydia|coxiella|leptospira|mycoplasma|rickettsia|treponema") ~ "Atypical/other bacteria",
    stringr::str_detect(x_clean, "staphylococcus|streptococcus|enterococcus|bacillus|corynebacter|lactobacillus|listeria|leuconostoc|micrococcus|nocardia|rhodococcus|stomatococcus|mrsa|vre|gram_positive|gram positive|gpc|coag_pos|coag_neg|coagneg") ~ "Gram positive bacteria",
    stringr::str_detect(x_clean, "acinetobacter|agrobacterium|alcaligenes|branhamelia|moraxella|pseudomonas|stenotrophomonas|xanthomonas|klebsiella|enterobacter|escherichia|serratia|haemophilus|citrobacter|proteus|neisseria|salmonella|shigella|campylobacter|burkholderia|cepacia|legionella|flavimonas|flavobacterium|helicobacter|methylobacterium|vibrio|esbl|cre|carbapenem|gram_negative|gram negative|gnr") ~ "Gram negative bacteria",
    stringr::str_detect(x_clean, "no_growth|no growth") ~ "No growth/negative",
    stringr::str_detect(x_clean, "bacteria|gram|cocci|bacilli|rods|flora") ~ "Other bacteria",
    TRUE ~ "Other/unspecified"
  )
}
