# ================================================================================================
# Diagnose Microbiology Fluid Category Shifts
#
# Goal:
#   Check whether apparent drops in selected culture types reflect ICU cohort logic or upstream
#   changes in microbiology_culture fluid_category assignment.
# ================================================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(glue)
  library(lubridate)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("utils/clif_io.R")

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

site_name <- clif_site_name
diagnosis_start <- safe_ts(Sys.getenv("DIAGNOSIS_START_DATE", unset = "2022-01-01"))
diagnosis_end <- safe_ts(Sys.getenv("DIAGNOSIS_END_DATE", unset = "2024-12-31")) + days(1) - seconds(1)
shift_date <- safe_ts(Sys.getenv("CATEGORY_SHIFT_DATE", unset = "2023-05-01"))
event_path <- Sys.getenv("ICU_CULTURE_EVENTS_PATH", unset = NA_character_)

target_categories <- str_split(
  Sys.getenv(
    "TARGET_FLUID_CATEGORIES",
    unset = "pleural_cavity_fluid,respiratory_tract_lower,respiratory_tract,other_unspecified,peritoneum,woundsite"
  ),
  ",",
  simplify = TRUE
) %>%
  as.character() %>%
  str_trim()

message("Using CLIF tables: ", clif_tables_path)
message("Diagnosis window: ", diagnosis_start, " to ", diagnosis_end)
message("Shift date: ", shift_date)

micro <- read_tbl("microbiology_culture") %>%
  transmute(
    patient_id,
    hospitalization_id,
    collect_dttm = safe_ts(collect_dttm),
    fluid_name = str_to_lower(str_trim(as.character(fluid_name))),
    fluid_category = str_to_lower(str_trim(as.character(fluid_category))),
    method_category = str_to_lower(str_trim(as.character(method_category)))
  ) %>%
  filter(
    method_category == "culture",
    !is.na(collect_dttm),
    collect_dttm >= diagnosis_start,
    collect_dttm <= diagnosis_end
  ) %>%
  mutate(
    month = floor_date(collect_dttm, "month"),
    period = if_else(collect_dttm < shift_date, "pre_shift", "post_shift")
  )

raw_category_ranges <- micro %>%
  filter(fluid_category %in% target_categories) %>%
  group_by(fluid_category) %>%
  summarise(
    n_rows = n(),
    first_collect_dttm = min(collect_dttm),
    last_collect_dttm = max(collect_dttm),
    n_months = n_distinct(month),
    .groups = "drop"
  ) %>%
  arrange(fluid_category)

raw_monthly_target_categories <- micro %>%
  filter(fluid_category %in% target_categories) %>%
  count(month, fluid_category) %>%
  arrange(month, fluid_category)

raw_period_category_counts <- micro %>%
  count(period, fluid_category, sort = TRUE) %>%
  group_by(period) %>%
  slice_max(n, n = 25, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(period, desc(n))

post_shift_other_names <- micro %>%
  filter(period == "post_shift", fluid_category == "other_unspecified") %>%
  count(fluid_name, sort = TRUE)

out_dir <- file.path("output", "diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

range_path <- file.path(out_dir, glue("raw_fluid_category_ranges_{site_name}_{stamp}.csv"))
monthly_path <- file.path(out_dir, glue("raw_monthly_target_fluid_categories_{site_name}_{stamp}.csv"))
period_path <- file.path(out_dir, glue("raw_period_fluid_category_counts_{site_name}_{stamp}.csv"))
other_name_path <- file.path(out_dir, glue("post_shift_other_unspecified_fluid_names_{site_name}_{stamp}.csv"))

write_csv(raw_category_ranges, range_path)
write_csv(raw_monthly_target_categories, monthly_path)
write_csv(raw_period_category_counts, period_path)
write_csv(post_shift_other_names, other_name_path)

message("Raw selected category ranges:")
print(raw_category_ranges, width = Inf)
message("")
message("Top post-shift other_unspecified fluid names:")
print(head(post_shift_other_names, 20), n = 20, width = Inf)
message("")
message("Wrote diagnostics:")
message(range_path)
message(monthly_path)
message(period_path)
message(other_name_path)
