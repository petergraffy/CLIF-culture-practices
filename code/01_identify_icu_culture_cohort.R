# ================================================================================================
# Identify ICU Culture Cohort
#
# Cohort:
#   Patients/hospitalizations with at least one microbiology culture collected during an ICU stay.
#
# Export:
#   Every microbiology culture row collected during an ICU interval for qualifying hospitalizations,
#   across all specimen/fluid types.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(janitor)
  library(lubridate)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("utils/clif_io.R")

site_name <- clif_site_name
tables_path <- clif_tables_path

safe_ts <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = tz))
  if (is.numeric(x)) {
    x2 <- ifelse(x > 1e12, x / 1000, x)
    return(as.POSIXct(x2, origin = "1970-01-01", tz = tz))
  }
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c("ymd_HMS", "ymd_HM", "ymd", "ymdTz", "ymdT", "mdy_HMS", "mdy_HM", "mdy"),
    tz = tz,
    quiet = TRUE
  ))
}

count_nonmissing <- function(x) sum(!is.na(x))

message("Using CLIF tables: ", tables_path)

hospitalization <- read_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm),
    admission_year = year(admission_dttm),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission))
  )

adt <- read_tbl("adt") %>%
  transmute(
    hospitalization_id,
    icu_in_dttm = safe_ts(in_dttm),
    icu_out_dttm_raw = safe_ts(out_dttm),
    location_category = str_to_lower(str_trim(as.character(location_category)))
  )

icu_intervals <- adt %>%
  filter(location_category == "icu", !is.na(icu_in_dttm)) %>%
  left_join(
    hospitalization %>% select(patient_id, hospitalization_id, admission_dttm, discharge_dttm),
    by = "hospitalization_id"
  ) %>%
  mutate(
    icu_out_dttm = coalesce(icu_out_dttm_raw, discharge_dttm),
    icu_interval_missing_out = is.na(icu_out_dttm_raw),
    icu_interval_id = row_number()
  ) %>%
  filter(!is.na(patient_id), !is.na(icu_out_dttm), icu_out_dttm > icu_in_dttm) %>%
  select(
    patient_id,
    hospitalization_id,
    icu_interval_id,
    admission_dttm,
    discharge_dttm,
    icu_in_dttm,
    icu_out_dttm,
    icu_interval_missing_out
  )

microbiology_culture <- read_tbl("microbiology_culture") %>%
  mutate(microbiology_row_id = row_number(), .before = 1) %>%
  transmute(
    microbiology_row_id,
    patient_id,
    hospitalization_id,
    organism_id = if ("organism_id" %in% names(.)) organism_id else NA,
    order_dttm = safe_ts(order_dttm),
    collect_dttm = safe_ts(collect_dttm),
    result_dttm = safe_ts(result_dttm),
    fluid_name = str_to_lower(str_trim(as.character(fluid_name))),
    fluid_category = str_to_lower(str_trim(as.character(fluid_category))),
    method_name = str_to_lower(str_trim(as.character(method_name))),
    method_category = str_to_lower(str_trim(as.character(method_category))),
    organism_name = str_to_lower(str_trim(as.character(organism_name))),
    organism_category = str_to_lower(str_trim(as.character(organism_category))),
    organism_group = str_to_lower(str_trim(as.character(organism_group)))
  ) %>%
  mutate(
    organism_group = coalesce(na_if(organism_group, ""), organism_category),
    no_growth = organism_group %in% c("no_growth", "no growth"),
    positive_culture = !is.na(organism_group) & !no_growth
  ) %>%
  filter(method_category == "culture", !is.na(collect_dttm))

icu_culture_rows <- microbiology_culture %>%
  inner_join(
    icu_intervals %>%
      select(patient_id, hospitalization_id, icu_interval_id, icu_in_dttm, icu_out_dttm, icu_interval_missing_out),
    by = c("patient_id", "hospitalization_id")
  ) %>%
  filter(collect_dttm >= icu_in_dttm, collect_dttm <= icu_out_dttm) %>%
  distinct(microbiology_row_id, icu_interval_id, .keep_all = TRUE) %>%
  arrange(patient_id, hospitalization_id, icu_in_dttm, collect_dttm, microbiology_row_id)

cohort_hospitalizations <- icu_culture_rows %>%
  distinct(patient_id, hospitalization_id) %>%
  left_join(hospitalization, by = c("patient_id", "hospitalization_id")) %>%
  arrange(patient_id, admission_dttm, hospitalization_id)

cohort_icu_intervals <- icu_culture_rows %>%
  distinct(patient_id, hospitalization_id, icu_interval_id, icu_in_dttm, icu_out_dttm, icu_interval_missing_out) %>%
  arrange(patient_id, hospitalization_id, icu_in_dttm)

culture_event_summary <- icu_culture_rows %>%
  group_by(
    patient_id,
    hospitalization_id,
    icu_interval_id,
    order_dttm,
    collect_dttm,
    fluid_name,
    fluid_category,
    method_name,
    method_category
  ) %>%
  summarise(
    n_culture_rows = n(),
    any_positive_culture = any(positive_culture, na.rm = TRUE),
    organism_groups = paste(sort(unique(na.omit(organism_group))), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(patient_id, hospitalization_id, collect_dttm, fluid_category, fluid_name)

fluid_summary <- icu_culture_rows %>%
  group_by(fluid_category, fluid_name) %>%
  summarise(
    n_culture_rows = n(),
    n_culture_events = n_distinct(paste(patient_id, hospitalization_id, icu_interval_id, collect_dttm, fluid_name, method_name)),
    n_positive_rows = sum(positive_culture, na.rm = TRUE),
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_culture_rows), fluid_category, fluid_name)

cohort_summary <- tibble(
  site_name = site_name,
  tables_path = tables_path,
  cohort_definition = "hospitalizations with at least one microbiology culture collected during an ICU interval",
  culture_export_definition = "all microbiology_culture rows with method_category == culture and collect_dttm within an ICU interval for qualifying hospitalizations",
  n_patients = n_distinct(cohort_hospitalizations$patient_id),
  n_hospitalizations = n_distinct(cohort_hospitalizations$hospitalization_id),
  n_icu_intervals_with_culture = n_distinct(cohort_icu_intervals$icu_interval_id),
  n_culture_rows = nrow(icu_culture_rows),
  n_culture_events = nrow(culture_event_summary),
  n_positive_culture_rows = sum(icu_culture_rows$positive_culture, na.rm = TRUE),
  n_positive_culture_events = sum(culture_event_summary$any_positive_culture, na.rm = TRUE),
  n_fluid_categories = n_distinct(icu_culture_rows$fluid_category, na.rm = TRUE),
  n_organism_groups = n_distinct(icu_culture_rows$organism_group, na.rm = TRUE),
  n_rows_with_collect_dttm = count_nonmissing(icu_culture_rows$collect_dttm),
  min_collect_dttm = suppressWarnings(min(icu_culture_rows$collect_dttm, na.rm = TRUE)),
  max_collect_dttm = suppressWarnings(max(icu_culture_rows$collect_dttm, na.rm = TRUE))
)

out_dir <- file.path("output", "cohort")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

culture_path <- file.path(out_dir, glue("icu_culture_rows_{site_name}_{stamp}.csv"))
event_path <- file.path(out_dir, glue("icu_culture_events_{site_name}_{stamp}.csv"))
hospitalization_path <- file.path(out_dir, glue("icu_culture_cohort_hospitalizations_{site_name}_{stamp}.csv"))
icu_interval_path <- file.path(out_dir, glue("icu_culture_cohort_icu_intervals_{site_name}_{stamp}.csv"))
summary_path <- file.path(out_dir, glue("icu_culture_cohort_summary_{site_name}_{stamp}.csv"))
fluid_path <- file.path(out_dir, glue("icu_culture_fluid_summary_{site_name}_{stamp}.csv"))

readr::write_csv(icu_culture_rows, culture_path)
readr::write_csv(culture_event_summary, event_path)
readr::write_csv(cohort_hospitalizations, hospitalization_path)
readr::write_csv(cohort_icu_intervals, icu_interval_path)
readr::write_csv(cohort_summary, summary_path)
readr::write_csv(fluid_summary, fluid_path)

message("ICU culture cohort summary:")
print(cohort_summary)
message("")
message("Top culture fluid categories:")
print(fluid_summary %>% select(fluid_category, fluid_name, n_culture_rows, n_positive_rows, n_hospitalizations, n_patients) %>% head(25), n = 25)
message("")
message("Wrote ICU culture rows: ", culture_path)
message("Wrote ICU culture events: ", event_path)
message("Wrote cohort hospitalizations: ", hospitalization_path)
message("Wrote cohort ICU intervals: ", icu_interval_path)
message("Wrote cohort summary: ", summary_path)
message("Wrote fluid summary: ", fluid_path)
