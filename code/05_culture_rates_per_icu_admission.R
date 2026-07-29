# ================================================================================================
# Monthly ICU Culture Collection Rates
#
# Question:
#   How often are cultures collected across sites, care settings, specimen types, and calendar time?
#
# Denominator:
#   All ICU admissions, defined as merged ICU ADT intervals. Back-to-back or overlapping ICU ADT
#   rows within the same hospitalization are counted as one ICU admission.
#
# Numerator:
#   ICU culture events collected during ICU time, collapsed from microbiology culture rows by
#   patient/hospitalization/ICU admission/collection time/specimen/method.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(glue)
  library(lubridate)
  library(readr)
  library(scales)
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

clean_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_sentence()
}

month_bar_width <- 25 * 24 * 60 * 60

site_name <- clif_site_name
tables_path <- clif_tables_path
study_start_date <- Sys.getenv("STUDY_START_DATE", unset = Sys.getenv("PLOT_START_DATE", unset = NA_character_))
study_end_date <- Sys.getenv("STUDY_END_DATE", unset = Sys.getenv("PLOT_END_DATE", unset = NA_character_))
top_n_types <- as.integer(Sys.getenv("TOP_N_CULTURE_TYPES", unset = "8"))

study_start_dttm <- if (!is.na(study_start_date) && nzchar(study_start_date)) safe_ts(study_start_date) else as.POSIXct(NA)
study_end_dttm <- if (!is.na(study_end_date) && nzchar(study_end_date)) safe_ts(study_end_date) + days(1) - seconds(1) else as.POSIXct(NA)

message("Using CLIF tables: ", tables_path)
if (!is.na(study_start_dttm)) message("Study start: ", study_start_dttm)
if (!is.na(study_end_dttm)) message("Study end: ", study_end_dttm)

hospitalization <- read_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm)
  )

adt <- read_tbl("adt") %>%
  transmute(
    hospitalization_id,
    icu_in_dttm = safe_ts(in_dttm),
    icu_out_dttm_raw = safe_ts(out_dttm),
    location_category = str_to_lower(str_trim(as.character(location_category)))
  )

icu_adt_intervals <- adt %>%
  filter(location_category == "icu", !is.na(icu_in_dttm)) %>%
  left_join(hospitalization, by = "hospitalization_id") %>%
  mutate(
    icu_out_dttm = coalesce(icu_out_dttm_raw, discharge_dttm),
    icu_interval_missing_out = is.na(icu_out_dttm_raw)
  ) %>%
  filter(!is.na(patient_id), !is.na(icu_out_dttm), icu_out_dttm > icu_in_dttm) %>%
  arrange(patient_id, hospitalization_id, icu_in_dttm, icu_out_dttm)

icu_admissions <- icu_adt_intervals %>%
  group_by(patient_id, hospitalization_id) %>%
  mutate(
    prior_max_icu_out_num = lag(cummax(as.numeric(icu_out_dttm))),
    new_icu_admission = is.na(prior_max_icu_out_num) | as.numeric(icu_in_dttm) > prior_max_icu_out_num,
    icu_admission_seq = cumsum(new_icu_admission)
  ) %>%
  group_by(patient_id, hospitalization_id, icu_admission_seq) %>%
  summarise(
    admission_dttm = first(admission_dttm),
    discharge_dttm = first(discharge_dttm),
    icu_in_dttm = min(icu_in_dttm, na.rm = TRUE),
    icu_out_dttm = max(icu_out_dttm, na.rm = TRUE),
    icu_interval_missing_out = any(icu_interval_missing_out, na.rm = TRUE),
    n_icu_adt_rows = n(),
    .groups = "drop"
  ) %>%
  arrange(patient_id, hospitalization_id, icu_in_dttm) %>%
  mutate(
    icu_admission_id = row_number(),
    icu_admission_month = floor_date(icu_in_dttm, "month")
  ) %>%
  filter(is.na(study_start_dttm) | icu_in_dttm >= study_start_dttm) %>%
  filter(is.na(study_end_dttm) | icu_in_dttm <= study_end_dttm)

microbiology_culture <- read_tbl("microbiology_culture") %>%
  mutate(microbiology_row_id = row_number(), .before = 1) %>%
  transmute(
    microbiology_row_id,
    patient_id,
    hospitalization_id,
    order_dttm = safe_ts(order_dttm),
    collect_dttm = safe_ts(collect_dttm),
    fluid_name = str_to_lower(str_trim(as.character(fluid_name))),
    fluid_category = str_to_lower(str_trim(as.character(fluid_category))),
    method_name = str_to_lower(str_trim(as.character(method_name))),
    method_category = str_to_lower(str_trim(as.character(method_category))),
    organism_category = if ("organism_category" %in% names(.)) {
      str_to_lower(str_trim(as.character(organism_category)))
    } else {
      NA_character_
    },
    organism_group = if ("organism_group" %in% names(.)) {
      str_to_lower(str_trim(as.character(organism_group)))
    } else {
      NA_character_
    }
  ) %>%
  mutate(
    organism_group = coalesce(na_if(organism_group, ""), organism_category),
    no_growth = organism_group %in% c("no_growth", "no growth"),
    positive_culture = !is.na(organism_group) & !no_growth
  ) %>%
  filter(method_category == "culture", !is.na(collect_dttm)) %>%
  filter(is.na(study_start_dttm) | collect_dttm >= study_start_dttm) %>%
  filter(is.na(study_end_dttm) | collect_dttm <= study_end_dttm)

icu_culture_rows <- microbiology_culture %>%
  inner_join(
    icu_admissions %>%
      select(patient_id, hospitalization_id, icu_admission_id, icu_in_dttm, icu_out_dttm),
    by = c("patient_id", "hospitalization_id"),
    relationship = "many-to-many"
  ) %>%
  filter(collect_dttm >= icu_in_dttm, collect_dttm <= icu_out_dttm) %>%
  distinct(microbiology_row_id, icu_admission_id, .keep_all = TRUE)

icu_culture_events <- icu_culture_rows %>%
  group_by(
    patient_id,
    hospitalization_id,
    icu_admission_id,
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
    .groups = "drop"
  ) %>%
  mutate(
    culture_month = floor_date(collect_dttm, "month"),
    culture_type = clean_label(coalesce(na_if(fluid_category, ""), "missing")),
    culture_type = if_else(culture_type %in% c("Other", "Other unspecified"), "Other", culture_type)
  )

if (nrow(icu_admissions) == 0) stop("No ICU admissions after filters.")

month_min <- min(c(icu_admissions$icu_admission_month, icu_culture_events$culture_month), na.rm = TRUE)
month_max <- max(c(icu_admissions$icu_admission_month, icu_culture_events$culture_month), na.rm = TRUE)
month_seq <- seq(month_min, month_max, by = "month")

monthly_icu_admissions <- icu_admissions %>%
  count(icu_admission_month, name = "n_icu_admissions") %>%
  complete(icu_admission_month = month_seq, fill = list(n_icu_admissions = 0L)) %>%
  rename(calendar_month = icu_admission_month)

top_types <- icu_culture_events %>%
  count(culture_type, sort = TRUE) %>%
  slice_head(n = top_n_types) %>%
  pull(culture_type)

monthly_culture_events_by_type <- icu_culture_events %>%
  mutate(culture_type_plot = if_else(culture_type %in% top_types, culture_type, "Other")) %>%
  group_by(calendar_month = culture_month, specimen_type = culture_type_plot) %>%
  summarise(
    n_culture_events = n(),
    n_icu_admissions_with_culture_type = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  group_by(calendar_month, specimen_type) %>%
  summarise(
    n_culture_events = sum(n_culture_events),
    n_icu_admissions_with_culture_type = sum(n_icu_admissions_with_culture_type),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    specimen_type,
    fill = list(n_culture_events = 0L, n_icu_admissions_with_culture_type = 0L)
  ) %>%
  left_join(monthly_icu_admissions, by = "calendar_month") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_events_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_culture_events / n_icu_admissions,
      NA_real_
    ),
    icu_admissions_with_culture_type_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_icu_admissions_with_culture_type / n_icu_admissions,
      NA_real_
    ),
    specimen_type = fct_relevel(factor(specimen_type), "Other", after = Inf)
  ) %>%
  arrange(calendar_month, specimen_type)

monthly_overall_rates <- icu_culture_events %>%
  group_by(calendar_month = culture_month) %>%
  summarise(
    n_culture_events = n(),
    n_icu_admissions_with_any_culture = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    fill = list(n_culture_events = 0L, n_icu_admissions_with_any_culture = 0L)
  ) %>%
  left_join(monthly_icu_admissions, by = "calendar_month") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_events_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_culture_events / n_icu_admissions,
      NA_real_
    ),
    icu_admissions_with_any_culture_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_icu_admissions_with_any_culture / n_icu_admissions,
      NA_real_
    )
  ) %>%
  arrange(calendar_month)

monthly_result_status_rates <- icu_culture_events %>%
  mutate(culture_result = if_else(any_positive_culture, "Positive", "Negative")) %>%
  group_by(calendar_month = culture_month, culture_result) %>%
  summarise(
    n_culture_events = n(),
    n_icu_admissions_with_culture_result = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    culture_result = c("Negative", "Positive"),
    fill = list(n_culture_events = 0L, n_icu_admissions_with_culture_result = 0L)
  ) %>%
  left_join(monthly_icu_admissions, by = "calendar_month") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_result = factor(culture_result, levels = c("Negative", "Positive")),
    culture_events_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_culture_events / n_icu_admissions,
      NA_real_
    ),
    icu_admissions_with_culture_result_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_icu_admissions_with_culture_result / n_icu_admissions,
      NA_real_
    )
  ) %>%
  arrange(calendar_month, culture_result)

culture_type_palette <- c(
  "Blood buffy" = "#E41A1C",
  "Respiratory tract" = "#377EB8",
  "Genito urinary tract" = "#4DAF4A",
  "Meninges csf" = "#984EA3",
  "Pleural cavity fluid" = "#FFFF33",
  "Respiratory tract lower" = "#A65628",
  "Woundsite" = "#F781BF",
  "Catheter tip" = "#FF7F00",
  "Other" = "#999999"
)
specimen_type_levels <- levels(monthly_culture_events_by_type$specimen_type)
extra_specimen_types <- setdiff(specimen_type_levels, names(culture_type_palette))
extra_palette <- if (length(extra_specimen_types) > 0) {
  setNames(hue_pal()(length(extra_specimen_types)), extra_specimen_types)
} else {
  character()
}
available_palette <- c(culture_type_palette, extra_palette)[specimen_type_levels]
result_status_palette <- c("Negative" = "#8F8F8F", "Positive" = "#B44E4E")

plot_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = grid::unit(3, "pt"),
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.caption.position = "plot"
  )

p_rate_stacked <- ggplot(
  monthly_culture_events_by_type,
  aes(calendar_month, culture_events_per_100_icu_admissions, fill = specimen_type)
) +
  geom_col(width = month_bar_width, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = available_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Monthly ICU Culture Collection Rates by Specimen Type",
    subtitle = "Culture events per 100 ICU admissions; Other includes less common and unspecified specimen types",
    x = NULL,
    y = "Culture events per 100 ICU admissions",
    fill = NULL
  ) +
  plot_theme

p_rate_lines <- monthly_culture_events_by_type %>%
  filter(specimen_type != "Other") %>%
  ggplot(aes(calendar_month, culture_events_per_100_icu_admissions, color = specimen_type)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = available_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Monthly ICU Culture Collection Rates for Major Specimen Types",
    x = NULL,
    y = "Culture events per 100 ICU admissions",
    color = NULL
  ) +
  plot_theme

p_result_status_rate <- ggplot(
  monthly_result_status_rates,
  aes(calendar_month, culture_events_per_100_icu_admissions, fill = culture_result)
) +
  geom_col(width = month_bar_width, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = result_status_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Monthly ICU Culture Collection Rates by Result Status",
    subtitle = "Positive and negative culture events per 100 ICU admissions",
    x = NULL,
    y = "Culture events per 100 ICU admissions",
    fill = NULL
  ) +
  plot_theme

out_dir <- file.path("output", "rates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

paths <- c(
  monthly_icu_admissions = file.path(out_dir, glue("monthly_icu_admissions_{site_name}_{stamp}.csv")),
  monthly_overall_rates = file.path(out_dir, glue("monthly_overall_culture_rates_per_100_icu_admissions_{site_name}_{stamp}.csv")),
  monthly_result_status_rates = file.path(out_dir, glue("monthly_result_status_culture_rates_per_100_icu_admissions_{site_name}_{stamp}.csv")),
  monthly_type_rates = file.path(out_dir, glue("monthly_specimen_type_culture_rates_per_100_icu_admissions_{site_name}_{stamp}.csv")),
  stacked_rate_plot = file.path(out_dir, glue("monthly_specimen_type_culture_rates_stacked_per_100_icu_admissions_{site_name}_{stamp}.png")),
  line_rate_plot = file.path(out_dir, glue("monthly_specimen_type_culture_rates_lines_per_100_icu_admissions_{site_name}_{stamp}.png")),
  result_status_rate_plot = file.path(out_dir, glue("monthly_result_status_culture_rates_stacked_per_100_icu_admissions_{site_name}_{stamp}.png"))
)

write_csv(monthly_icu_admissions, paths[["monthly_icu_admissions"]])
write_csv(monthly_overall_rates, paths[["monthly_overall_rates"]])
write_csv(monthly_result_status_rates, paths[["monthly_result_status_rates"]])
write_csv(monthly_culture_events_by_type, paths[["monthly_type_rates"]])

ggsave(paths[["stacked_rate_plot"]], p_rate_stacked, width = 12, height = 7, dpi = 300)
ggsave(paths[["line_rate_plot"]], p_rate_lines, width = 12, height = 7, dpi = 300)
ggsave(paths[["result_status_rate_plot"]], p_result_status_rate, width = 12, height = 7, dpi = 300)

message("Monthly ICU admission denominator summary:")
print(monthly_icu_admissions %>% summarise(
  first_month = min(calendar_month),
  last_month = max(calendar_month),
  total_icu_admissions = sum(n_icu_admissions),
  median_monthly_icu_admissions = median(n_icu_admissions),
  min_monthly_icu_admissions = min(n_icu_admissions),
  max_monthly_icu_admissions = max(n_icu_admissions)
), width = Inf)

message("")
message("Overall culture collection rate summary:")
print(monthly_overall_rates %>% summarise(
  total_culture_events = sum(n_culture_events),
  median_monthly_events_per_100_icu_admissions = median(culture_events_per_100_icu_admissions, na.rm = TRUE),
  min_monthly_events_per_100_icu_admissions = min(culture_events_per_100_icu_admissions, na.rm = TRUE),
  max_monthly_events_per_100_icu_admissions = max(culture_events_per_100_icu_admissions, na.rm = TRUE)
), width = Inf)

message("")
message("Culture result status rate summary:")
print(monthly_result_status_rates %>%
  group_by(culture_result) %>%
  summarise(
    total_culture_events = sum(n_culture_events),
    median_monthly_events_per_100_icu_admissions = median(culture_events_per_100_icu_admissions, na.rm = TRUE),
    min_monthly_events_per_100_icu_admissions = min(culture_events_per_100_icu_admissions, na.rm = TRUE),
    max_monthly_events_per_100_icu_admissions = max(culture_events_per_100_icu_admissions, na.rm = TRUE),
    .groups = "drop"
  ), width = Inf)

message("")
message("Specimen types displayed separately:")
print(tibble(specimen_type = top_types), n = Inf)
message("")
message("Wrote outputs:")
print(paths)
