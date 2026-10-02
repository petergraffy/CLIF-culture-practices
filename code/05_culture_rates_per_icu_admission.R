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


clean_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_sentence()
}

month_bar_width <- 25 * 24 * 60 * 60

site_name <- clif_site_name
tables_path <- clif_tables_path
study_start_date <- study_settings$study_start_date
study_end_date <- study_settings$study_end_date
top_n_types <- as.integer(Sys.getenv("TOP_N_CULTURE_TYPES", unset = "8"))

study_start_dttm <- if (!is.na(study_start_date) && nzchar(study_start_date)) safe_ts(study_start_date) else as.POSIXct(NA)
study_end_dttm <- if (!is.na(study_end_date) && nzchar(study_end_date)) safe_ts(study_end_date) + days(1) else as.POSIXct(NA)

message("Using CLIF tables: ", tables_path)
if (!is.na(study_start_dttm)) message("Study start: ", study_start_dttm)
if (!is.na(study_end_dttm)) message("Study end: ", study_end_dttm)

culture_data <- read_culture_data(study_start_dttm, study_end_dttm)
icu_admissions <- culture_data$icu_admissions
icu_culture_rows <- culture_data$rows
icu_culture_events <- culture_data$events %>% mutate(culture_month = floor_date(collect_dttm, "month"), culture_type = clean_label(coalesce(fluid_category, "missing")), culture_type = if_else(culture_type %in% c("Other", "Other unspecified"), "Other", culture_type))

if (nrow(icu_admissions) == 0) stop("No ICU admissions after filters.")

month_min <- floor_date(min(icu_admissions$icu_in_dttm_clipped), "month")
month_max <- floor_date(max(icu_admissions$icu_out_dttm_clipped - seconds(1)), "month")
month_seq <- seq(month_min, month_max, by = "month")

monthly_icu_admissions <- monthly_icu_denominators(icu_admissions, month_seq) %>%
  select(calendar_month, n_icu_admissions)

top_types <- icu_culture_events %>%
  count(culture_type, sort = TRUE) %>%
  slice_head(n = top_n_types) %>%
  pull(culture_type)

monthly_culture_events_by_type <- icu_culture_events %>%
  mutate(culture_type_plot = if_else(culture_type %in% top_types, culture_type, "Other")) %>%
  group_by(calendar_month = culture_month, specimen_type = culture_type_plot) %>%
  summarise(
    n_culture_events = n(),
    n_active_icu_admissions_with_culture_type = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  group_by(calendar_month, specimen_type) %>%
  summarise(
    n_culture_events = sum(n_culture_events),
    n_active_icu_admissions_with_culture_type = sum(n_active_icu_admissions_with_culture_type),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    specimen_type,
    fill = list(n_culture_events = 0L, n_active_icu_admissions_with_culture_type = 0L)
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
    specimen_type = fct_relevel(factor(specimen_type), "Other", after = Inf)
  ) %>%
  arrange(calendar_month, specimen_type)

monthly_overall_rates <- icu_culture_events %>%
  group_by(calendar_month = culture_month) %>%
  summarise(
    n_culture_events = n(),
    n_active_icu_admissions_with_any_culture = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    fill = list(n_culture_events = 0L, n_active_icu_admissions_with_any_culture = 0L)
  ) %>%
  left_join(monthly_icu_admissions, by = "calendar_month") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_events_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_culture_events / n_icu_admissions,
      NA_real_
    )
  ) %>%
  arrange(calendar_month)

monthly_result_status_rates <- icu_culture_events %>%
  mutate(culture_result = result_status) %>%
  group_by(calendar_month = culture_month, culture_result) %>%
  summarise(
    n_culture_events = n(),
    n_active_icu_admissions_with_culture_result = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    culture_result = culture_result_levels,
    fill = list(n_culture_events = 0L, n_active_icu_admissions_with_culture_result = 0L)
  ) %>%
  left_join(monthly_icu_admissions, by = "calendar_month") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_result = factor(culture_result, levels = culture_result_levels),
    culture_events_per_100_icu_admissions = if_else(
      n_icu_admissions > 0,
      100 * n_culture_events / n_icu_admissions,
      NA_real_
    )
  ) %>%
  arrange(calendar_month, culture_result)

# Admission-cohort proportions: assign every observed culture to its stay's entry month.
admission_cohort <- icu_admissions %>%
  filter(is.na(study_start_dttm) | icu_in_dttm >= study_start_dttm) %>%
  filter(is.na(study_end_dttm) | icu_in_dttm < study_end_dttm) %>%
  left_join(icu_culture_events %>% group_by(icu_admission_id) %>% summarise(n_observed_culture_events = n(), .groups = "drop"), by = "icu_admission_id") %>%
  mutate(n_observed_culture_events = coalesce(n_observed_culture_events, 0L), followup_truncated = !is.na(study_end_dttm) & icu_out_dttm > study_end_dttm) %>%
  group_by(calendar_month = icu_admission_month) %>%
  summarise(n_icu_admissions = n(), n_icu_admissions_with_any_culture = sum(n_observed_culture_events > 0), n_followup_truncated = sum(followup_truncated), .groups = "drop") %>%
  mutate(site_name = site_name, proportion_icu_admissions_cultured = n_icu_admissions_with_any_culture / n_icu_admissions)

specimen_type_levels <- levels(monthly_culture_events_by_type$specimen_type)
available_palette <- clif_complete_specimen_palette(specimen_type_levels)
result_status_palette <- c("Negative/no growth" = "#8F8F8F", "Positive" = "#B44E4E", "Mixed/contaminated" = "#D49B35", "Indeterminate" = "#6688AA")

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
    subtitle = "Culture events by classified result per 100 ICU admissions",
    x = NULL,
    y = "Culture events per 100 ICU admissions",
    fill = NULL
  ) +
  plot_theme

out_dir <- project_output_dir("rates")
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

write_csv(admission_cohort, file.path(out_dir, glue("monthly_admission_cohort_culture_proportions_{site_name}_{stamp}.csv")))
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
