# ================================================================================================
# Identify ICU Culture Cohort
#
# Cohort:
#   Patients/hospitalizations with at least one microbiology culture collected during an ICU stay.
#
# Export:
#   Aggregate, non-PHI cohort and fluid summaries under output/cohort/.
#   Row-level culture extracts are private intermediates and are written only to data/intermediate/
#   when WRITE_ROW_LEVEL_INTERMEDIATES=true.
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
study_start_date <- study_settings$study_start_date
study_end_date <- study_settings$study_end_date
write_row_level_intermediates <- tolower(Sys.getenv("WRITE_ROW_LEVEL_INTERMEDIATES", unset = "true")) %in% c("true", "1", "yes", "y")


count_nonmissing <- function(x) sum(!is.na(x))

message("Using CLIF tables: ", tables_path)

study_start_dttm <- if (!is.na(study_start_date) && nzchar(study_start_date)) {
  safe_ts(study_start_date)
} else {
  as.POSIXct(NA)
}
study_end_dttm <- if (!is.na(study_end_date) && nzchar(study_end_date)) {
  safe_ts(study_end_date) + days(1)
} else {
  as.POSIXct(NA)
}

if (!is.na(study_start_dttm)) message("Study start: ", study_start_dttm)
if (!is.na(study_end_dttm)) message("Study end: ", study_end_dttm)

culture_data <- read_culture_data(study_start_dttm, study_end_dttm)
hospitalization <- culture_data$hospitalization %>% mutate(across(c(admission_dttm, discharge_dttm), safe_ts), admission_year = year(admission_dttm))
icu_intervals <- culture_data$icu_admissions
icu_culture_rows <- culture_data$rows

cohort_hospitalizations <- icu_culture_rows %>%
  distinct(patient_id, hospitalization_id) %>%
  left_join(hospitalization, by = c("patient_id", "hospitalization_id")) %>%
  arrange(patient_id, admission_dttm, hospitalization_id)

cohort_icu_intervals <- icu_culture_rows %>%
  distinct(patient_id, hospitalization_id, icu_interval_id, icu_in_dttm, icu_out_dttm, icu_interval_missing_out) %>%
  arrange(patient_id, hospitalization_id, icu_in_dttm)

culture_event_summary <- culture_data$events

fluid_summary <- icu_culture_rows %>%
  group_by(fluid_category, fluid_name) %>%
  summarise(
    n_culture_rows = n(),
    n_culture_events = n_distinct(culture_event_id),
    n_positive_rows = sum(positive_culture, na.rm = TRUE),
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_culture_rows), fluid_category, fluid_name)

cohort_summary <- tibble(
  site_name = site_name,
  study_start_date = if_else(is.na(study_start_dttm), NA_character_, as.character(as.Date(study_start_dttm))),
  study_end_date = if_else(is.na(study_end_dttm), NA_character_, as.character(as.Date(study_end_dttm - days(1)))),
  cohort_definition = "hospitalizations with at least one microbiology culture collected during an ICU interval",
  culture_event_definition = "unique patient/hospitalization/merged ICU stay/order and collection time/source specimen/source method culture events with method_category == culture collected during ICU time",
  n_patients = n_distinct(cohort_hospitalizations$patient_id),
  n_hospitalizations = n_distinct(cohort_hospitalizations$hospitalization_id),
  n_icu_intervals_with_culture = n_distinct(cohort_icu_intervals$icu_interval_id),
  n_culture_rows = nrow(icu_culture_rows),
  n_culture_events = nrow(culture_event_summary),
  n_positive_culture_rows = sum(icu_culture_rows$positive_culture, na.rm = TRUE),
  n_positive_culture_events = sum(culture_event_summary$any_positive_culture, na.rm = TRUE),
  n_fluid_categories = n_distinct(icu_culture_rows$fluid_category, na.rm = TRUE),
  n_organism_groups = n_distinct(icu_culture_rows$organism_group, na.rm = TRUE),
  n_rows_with_collection_time = count_nonmissing(icu_culture_rows$collect_dttm),
  first_collect_date = as.Date(suppressWarnings(min(icu_culture_rows$collect_dttm, na.rm = TRUE))),
  last_collect_date = as.Date(suppressWarnings(max(icu_culture_rows$collect_dttm, na.rm = TRUE)))
)

out_dir <- project_output_dir("cohort")
intermediate_dir <- project_intermediate_path("cohort")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

summary_path <- file.path(out_dir, glue("icu_culture_cohort_summary_{site_name}_{stamp}.csv"))
fluid_path <- file.path(out_dir, glue("icu_culture_fluid_summary_{site_name}_{stamp}.csv"))

readr::write_csv(cohort_summary, summary_path)
readr::write_csv(fluid_summary, fluid_path)

if (write_row_level_intermediates) {
  intermediate_dir <- project_intermediate_dir("cohort")
  culture_path <- file.path(intermediate_dir, glue("icu_culture_rows_{site_name}_{stamp}.csv"))
  event_path <- file.path(intermediate_dir, glue("icu_culture_events_{site_name}_{stamp}.csv"))
  hospitalization_path <- file.path(intermediate_dir, glue("icu_culture_cohort_hospitalizations_{site_name}_{stamp}.csv"))
  icu_interval_path <- file.path(intermediate_dir, glue("icu_culture_cohort_icu_intervals_{site_name}_{stamp}.csv"))

  readr::write_csv(icu_culture_rows, culture_path)
  readr::write_csv(culture_event_summary, event_path)
  readr::write_csv(cohort_hospitalizations, hospitalization_path)
  readr::write_csv(cohort_icu_intervals, icu_interval_path)
}

message("ICU culture cohort summary:")
print(cohort_summary)
message("")
message("Top culture fluid categories:")
print(fluid_summary %>% select(fluid_category, fluid_name, n_culture_rows, n_positive_rows, n_hospitalizations, n_patients) %>% head(25), n = 25)
message("")
message("Wrote cohort summary: ", summary_path)
message("Wrote fluid summary: ", fluid_path)
if (write_row_level_intermediates) {
  message("")
  message("Wrote private row-level intermediates under ignored path: ", intermediate_dir)
  message("Do not share or copy data/intermediate/ into pooled site exports.")
}
