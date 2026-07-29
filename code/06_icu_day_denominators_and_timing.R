# ================================================================================================
# ICU-Day Culture Denominators and Timing
#
# Questions:
#   1. How often are cultures collected per ICU day?
#   2. When during an ICU admission are cultures first collected?
#
# Denominator:
#   All ICU admissions, defined as merged ICU ADT intervals. ICU days are allocated to calendar
#   months using the overlap between each ICU interval and each month.
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
  library(purrr)
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

drop_negative_organism_name <- function(x) {
  str_detect(
    coalesce(x, ""),
    regex("^(no|none|not) .*isolated|no growth|no .* detected|negative for", ignore_case = TRUE)
  )
}

clip_ts <- function(x, lower = as.POSIXct(NA), upper = as.POSIXct(NA)) {
  if (!is.na(lower)) x <- pmax(x, lower)
  if (!is.na(upper)) x <- pmin(x, upper)
  x
}

month_bar_width <- 25 * 24 * 60 * 60

site_name <- clif_site_name
tables_path <- clif_tables_path
study_start_date <- Sys.getenv("STUDY_START_DATE", unset = Sys.getenv("PLOT_START_DATE", unset = NA_character_))
study_end_date <- Sys.getenv("STUDY_END_DATE", unset = Sys.getenv("PLOT_END_DATE", unset = NA_character_))
top_n_types <- as.integer(Sys.getenv("TOP_N_CULTURE_TYPES", unset = "8"))
timing_max_day <- as.integer(Sys.getenv("TIMING_MAX_ICU_DAY", unset = "14"))
timing_max_hour <- as.integer(Sys.getenv("TIMING_MAX_ICU_HOUR", unset = "168"))
top_n_organisms <- as.integer(Sys.getenv("TOP_N_TIMING_ORGANISMS", unset = "10"))

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
    icu_admission_month = floor_date(icu_in_dttm, "month"),
    icu_in_dttm_clipped = clip_ts(icu_in_dttm, study_start_dttm, study_end_dttm),
    icu_out_dttm_clipped = clip_ts(icu_out_dttm, study_start_dttm, study_end_dttm)
  ) %>%
  filter(icu_out_dttm_clipped > icu_in_dttm_clipped) %>%
  mutate(icu_los_days = as.numeric(difftime(icu_out_dttm_clipped, icu_in_dttm_clipped, units = "days")))

if (nrow(icu_admissions) == 0) stop("No ICU admissions after filters.")

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
    },
    organism_name = if ("organism_name" %in% names(.)) {
      str_to_lower(str_trim(as.character(organism_name)))
    } else {
      NA_character_
    }
  ) %>%
  mutate(
    organism_group = coalesce(na_if(organism_group, ""), organism_category),
    organism_name = coalesce(na_if(organism_name, ""), organism_group),
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
    icu_in_dttm,
    icu_out_dttm,
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
    calendar_month = floor_date(collect_dttm, "month"),
    specimen_type = clean_label(coalesce(na_if(fluid_category, ""), "missing")),
    specimen_type = if_else(specimen_type %in% c("Other", "Other unspecified"), "Other", specimen_type),
    hours_since_icu_admit = as.numeric(difftime(collect_dttm, icu_in_dttm, units = "hours")),
    icu_day = floor(hours_since_icu_admit / 24) + 1L
  )

month_min <- floor_date(min(icu_admissions$icu_in_dttm_clipped, na.rm = TRUE), "month")
month_max <- floor_date(max(icu_admissions$icu_out_dttm_clipped, na.rm = TRUE), "month")
month_seq <- seq(month_min, month_max, by = "month")

monthly_icu_days <- icu_admissions %>%
  mutate(month_start = map(icu_in_dttm_clipped, ~ seq(floor_date(.x, "month"), month_max, by = "month"))) %>%
  select(icu_admission_id, month_start, icu_in_dttm_clipped, icu_out_dttm_clipped) %>%
  unnest(month_start) %>%
  mutate(
    month_end = month_start %m+% months(1),
    overlap_start = pmax(icu_in_dttm_clipped, month_start),
    overlap_end = pmin(icu_out_dttm_clipped, month_end),
    icu_days = as.numeric(difftime(overlap_end, overlap_start, units = "days"))
  ) %>%
  filter(icu_days > 0) %>%
  group_by(calendar_month = month_start) %>%
  summarise(
    n_icu_days = sum(icu_days),
    n_icu_admissions_contributing_days = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    fill = list(n_icu_days = 0, n_icu_admissions_contributing_days = 0L)
  )

monthly_icu_admissions <- icu_admissions %>%
  count(icu_admission_month, name = "n_icu_admissions") %>%
  complete(icu_admission_month = month_seq, fill = list(n_icu_admissions = 0L)) %>%
  rename(calendar_month = icu_admission_month)

top_types <- icu_culture_events %>%
  count(specimen_type, sort = TRUE) %>%
  slice_head(n = top_n_types) %>%
  pull(specimen_type)

monthly_events_by_type_per_icu_day <- icu_culture_events %>%
  mutate(specimen_type_plot = if_else(specimen_type %in% top_types, specimen_type, "Other")) %>%
  group_by(calendar_month, specimen_type = specimen_type_plot) %>%
  summarise(
    n_culture_events = n(),
    n_positive_culture_events = sum(any_positive_culture, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(calendar_month, specimen_type) %>%
  summarise(
    n_culture_events = sum(n_culture_events),
    n_positive_culture_events = sum(n_positive_culture_events),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    specimen_type,
    fill = list(n_culture_events = 0L, n_positive_culture_events = 0L)
  ) %>%
  left_join(monthly_icu_days, by = "calendar_month") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_events_per_100_icu_days = if_else(n_icu_days > 0, 100 * n_culture_events / n_icu_days, NA_real_),
    positive_culture_events_per_100_icu_days = if_else(n_icu_days > 0, 100 * n_positive_culture_events / n_icu_days, NA_real_),
    specimen_type = fct_relevel(factor(specimen_type), "Other", after = Inf)
  ) %>%
  arrange(calendar_month, specimen_type)

monthly_overall_per_icu_day <- icu_culture_events %>%
  group_by(calendar_month) %>%
  summarise(
    n_culture_events = n(),
    n_positive_culture_events = sum(any_positive_culture, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    calendar_month = month_seq,
    fill = list(n_culture_events = 0L, n_positive_culture_events = 0L)
  ) %>%
  left_join(monthly_icu_days, by = "calendar_month") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_events_per_100_icu_days = if_else(n_icu_days > 0, 100 * n_culture_events / n_icu_days, NA_real_),
    positive_culture_events_per_100_icu_days = if_else(n_icu_days > 0, 100 * n_positive_culture_events / n_icu_days, NA_real_)
  ) %>%
  arrange(calendar_month)

first_culture_overall <- icu_culture_events %>%
  arrange(icu_admission_id, collect_dttm, specimen_type) %>%
  group_by(icu_admission_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(icu_admission_id, first_culture_dttm = collect_dttm, first_specimen_type = specimen_type, first_culture_hours = hours_since_icu_admit)

first_culture_timing <- icu_admissions %>%
  select(icu_admission_id, patient_id, hospitalization_id, icu_in_dttm, icu_out_dttm, icu_los_days) %>%
  left_join(first_culture_overall, by = "icu_admission_id") %>%
  mutate(
    any_icu_culture = !is.na(first_culture_dttm),
    first_culture_day = floor(first_culture_hours / 24) + 1L,
    timing_bin = case_when(
      is.na(first_culture_hours) ~ "No ICU culture",
      first_culture_hours < 6 ~ "0-6 hours",
      first_culture_hours < 24 ~ "6-24 hours",
      first_culture_hours < 48 ~ "ICU day 2",
      first_culture_hours < 72 ~ "ICU day 3",
      first_culture_hours < 168 ~ "ICU days 4-7",
      TRUE ~ "After ICU day 7"
    ),
    timing_bin = factor(
      timing_bin,
      levels = c("0-6 hours", "6-24 hours", "ICU day 2", "ICU day 3", "ICU days 4-7", "After ICU day 7", "No ICU culture")
    ),
    first_specimen_type = coalesce(first_specimen_type, "No ICU culture")
  )

timing_bin_summary <- first_culture_timing %>%
  count(timing_bin, first_specimen_type, name = "n_icu_admissions") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    percent_icu_admissions = 100 * n_icu_admissions / nrow(first_culture_timing)
  ) %>%
  arrange(timing_bin, desc(n_icu_admissions))

first_culture_timing_summary <- first_culture_timing %>%
  summarise(
    site_name = site_name,
    care_setting = "ICU",
    n_icu_admissions = n(),
    n_with_icu_culture = sum(any_icu_culture),
    percent_with_icu_culture = 100 * mean(any_icu_culture),
    median_hours_to_first_culture = median(first_culture_hours, na.rm = TRUE),
    p25_hours_to_first_culture = quantile(first_culture_hours, 0.25, na.rm = TRUE),
    p75_hours_to_first_culture = quantile(first_culture_hours, 0.75, na.rm = TRUE)
  )

icu_day_at_risk <- tibble(icu_day = seq_len(timing_max_day)) %>%
  mutate(
    n_icu_admissions_at_risk = map_int(
      icu_day,
      ~ sum(icu_admissions$icu_los_days > (.x - 1), na.rm = TRUE)
    )
  )

icu_day_event_rates <- icu_culture_events %>%
  filter(icu_day >= 1, icu_day <= timing_max_day) %>%
  count(icu_day, name = "n_culture_events") %>%
  complete(icu_day = seq_len(timing_max_day), fill = list(n_culture_events = 0L)) %>%
  left_join(icu_day_at_risk, by = "icu_day") %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    culture_events_per_100_icu_admissions_at_risk =
      100 * n_culture_events / n_icu_admissions_at_risk
  )

cumulative_first_culture_by_day <- first_culture_timing %>%
  crossing(icu_day = seq_len(timing_max_day)) %>%
  group_by(icu_day) %>%
  summarise(
    n_icu_admissions_with_first_culture_by_day = sum(any_icu_culture & first_culture_day <= icu_day, na.rm = TRUE),
    n_icu_admissions = n_distinct(icu_admission_id),
    .groups = "drop"
  ) %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    icu_admissions_with_first_culture_per_100_icu_admissions =
      100 * n_icu_admissions_with_first_culture_by_day / n_icu_admissions
  )

cumulative_culture_events_by_type_hour <- icu_culture_events %>%
  mutate(specimen_type_plot = if_else(specimen_type %in% top_types, specimen_type, "Other")) %>%
  group_by(specimen_type = specimen_type_plot) %>%
  summarise(n_total_culture_events = n(), .groups = "drop") %>%
  crossing(icu_hour = 0:timing_max_hour) %>%
  left_join(
    icu_culture_events %>%
      mutate(
        specimen_type = if_else(specimen_type %in% top_types, specimen_type, "Other"),
        event_icu_hour = ceiling(hours_since_icu_admit)
      ) %>%
      filter(event_icu_hour >= 0, event_icu_hour <= timing_max_hour) %>%
      count(specimen_type, event_icu_hour, name = "n_culture_events") %>%
      group_by(specimen_type) %>%
      complete(event_icu_hour = 0:timing_max_hour, fill = list(n_culture_events = 0L)) %>%
      arrange(specimen_type, event_icu_hour) %>%
      mutate(n_culture_events_by_hour = cumsum(n_culture_events)) %>%
      ungroup(),
    by = c("specimen_type", "icu_hour" = "event_icu_hour")
  ) %>%
  mutate(
    n_culture_events = coalesce(n_culture_events, 0L),
    n_culture_events_by_hour = coalesce(n_culture_events_by_hour, 0L),
    site_name = site_name,
    care_setting = "ICU",
    percent_culture_events_by_hour = 100 * n_culture_events_by_hour / n_total_culture_events,
    specimen_type = fct_relevel(factor(specimen_type), "Other", after = Inf)
  ) %>%
  arrange(specimen_type, icu_hour)

positive_organism_detections <- icu_culture_rows %>%
  mutate(
    organism_name = coalesce(na_if(organism_name, ""), organism_group, organism_category),
    organism_label = str_to_sentence(str_squish(organism_name)),
    explicit_negative_name = drop_negative_organism_name(organism_name),
    hours_since_icu_admit = as.numeric(difftime(collect_dttm, icu_in_dttm, units = "hours"))
  ) %>%
  filter(positive_culture, !explicit_negative_name, !is.na(organism_name), !is.na(hours_since_icu_admit)) %>%
  distinct(
    patient_id,
    hospitalization_id,
    icu_admission_id,
    collect_dttm,
    fluid_name,
    fluid_category,
    method_name,
    organism_name,
    organism_label,
    hours_since_icu_admit
  )

top_organisms <- positive_organism_detections %>%
  count(organism_name, organism_label, sort = TRUE) %>%
  slice_head(n = top_n_organisms)

first_organism_detection_by_admission <- positive_organism_detections %>%
  inner_join(
    top_organisms %>% select(organism_name, organism_label),
    by = c("organism_name", "organism_label")
  ) %>%
  group_by(icu_admission_id, organism_name, organism_label) %>%
  summarise(
    first_detection_hours = min(hours_since_icu_admit, na.rm = TRUE),
    .groups = "drop"
  )

cumulative_top_organism_incidence_hour <- top_organisms %>%
  select(organism_name, organism_label, n_total_detection_events = n) %>%
  crossing(icu_hour = 0:timing_max_hour) %>%
  left_join(
    first_organism_detection_by_admission %>%
      mutate(first_detection_hour = ceiling(first_detection_hours)) %>%
      filter(first_detection_hour >= 0, first_detection_hour <= timing_max_hour) %>%
      count(organism_name, organism_label, first_detection_hour, name = "n_first_detection_icu_admissions") %>%
      group_by(organism_name, organism_label) %>%
      complete(first_detection_hour = 0:timing_max_hour, fill = list(n_first_detection_icu_admissions = 0L)) %>%
      arrange(organism_name, organism_label, first_detection_hour) %>%
      mutate(n_icu_admissions_with_organism_by_hour = cumsum(n_first_detection_icu_admissions)) %>%
      ungroup(),
    by = c("organism_name", "organism_label", "icu_hour" = "first_detection_hour")
  ) %>%
  mutate(
    n_first_detection_icu_admissions = coalesce(n_first_detection_icu_admissions, 0L),
    n_icu_admissions_with_organism_by_hour = coalesce(n_icu_admissions_with_organism_by_hour, 0L),
    n_icu_admissions = nrow(icu_admissions),
    site_name = site_name,
    care_setting = "ICU",
    organism_first_collected_per_100_icu_admissions =
      100 * n_icu_admissions_with_organism_by_hour / n_icu_admissions,
    organism_label = fct_reorder(organism_label, n_total_detection_events, .desc = TRUE)
  ) %>%
  arrange(organism_label, icu_hour)

culture_type_palette <- c(
  "Blood buffy" = "#FB8072",
  "Respiratory tract" = "#BC80BD",
  "Genito urinary tract" = "#FDB462",
  "Meninges csf" = "#8DD3C7",
  "Pleural cavity fluid" = "#BEBADA",
  "Respiratory tract lower" = "#B3DE69",
  "Woundsite" = "#FCCDE5",
  "Catheter tip" = "#80B1D3",
  "Other" = "#BDBDBD",
  "No ICU culture" = "#C7C7C7"
)
specimen_levels <- levels(monthly_events_by_type_per_icu_day$specimen_type)
extra_specimen_types <- setdiff(specimen_levels, names(culture_type_palette))
extra_palette <- if (length(extra_specimen_types) > 0) {
  setNames(
    grDevices::colorRampPalette(c("#FB8072", "#BC80BD", "#FDB462", "#8DD3C7", "#80B1D3", "#B3DE69", "#BEBADA", "#FCCDE5"))(length(extra_specimen_types)),
    extra_specimen_types
  )
} else {
  character()
}
available_palette <- c(culture_type_palette, extra_palette)[specimen_levels]

plot_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = grid::unit(3, "pt"),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.caption.position = "plot"
  )

p_icu_day_rate_stacked <- ggplot(
  monthly_events_by_type_per_icu_day,
  aes(calendar_month, culture_events_per_100_icu_days, fill = specimen_type)
) +
  geom_col(width = month_bar_width, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = available_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Monthly ICU Culture Collection Rates by Specimen Type",
    x = NULL,
    y = "Culture events per 100 ICU days",
    fill = NULL
  ) +
  plot_theme

p_timing_bins <- timing_bin_summary %>%
  filter(timing_bin != "No ICU culture") %>%
  mutate(
    first_specimen_type = if_else(first_specimen_type %in% top_types, first_specimen_type, "Other"),
    first_specimen_type = fct_relevel(factor(first_specimen_type), "Other", after = Inf)
  ) %>%
  group_by(timing_bin, first_specimen_type) %>%
  summarise(
    percent_icu_admissions = sum(percent_icu_admissions),
    n_icu_admissions = sum(n_icu_admissions),
    .groups = "drop"
  ) %>%
  ggplot(aes(timing_bin, percent_icu_admissions, fill = first_specimen_type)) +
  geom_col(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = c(culture_type_palette, extra_palette)) +
  scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, NA)) +
  labs(
    title = "Timing of First ICU Culture",
    x = NULL,
    y = "ICU admissions",
    fill = NULL
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_cumulative_first_culture <- ggplot(
  cumulative_first_culture_by_day,
  aes(icu_day, icu_admissions_with_first_culture_per_100_icu_admissions)
) +
  geom_line(color = "#2F6C99", linewidth = 0.85) +
  geom_point(color = "#2F6C99", size = 1.6) +
  scale_x_continuous(breaks = seq_len(timing_max_day)) +
  scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, NA)) +
  labs(
    title = "Cumulative First ICU Culture by ICU Day",
    x = "ICU day",
    y = "ICU admissions with first culture"
  ) +
  plot_theme

p_icu_day_event_rates <- ggplot(
  icu_day_event_rates,
  aes(icu_day, culture_events_per_100_icu_admissions_at_risk)
) +
  geom_col(fill = "#4C78A8", width = 0.8) +
  scale_x_continuous(breaks = seq_len(timing_max_day)) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Culture Event Timing by ICU Day",
    x = "ICU day",
    y = "Culture events per 100 ICU admissions at risk"
  ) +
  plot_theme +
  theme(legend.position = "none")

p_cumulative_events_by_type <- ggplot(
  cumulative_culture_events_by_type_hour,
  aes(icu_hour, percent_culture_events_by_hour, color = specimen_type)
) +
  geom_line(linewidth = 0.85) +
  scale_color_manual(values = available_palette) +
  scale_x_continuous(
    breaks = seq(0, timing_max_hour, by = 24),
    labels = function(x) x / 24
  ) +
  scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, NA)) +
  labs(
    title = "Cumulative ICU Culture Events by Time Since ICU Admission and Specimen Type",
    x = "Days since ICU admission",
    y = "Culture events collected",
    color = NULL
  ) +
  plot_theme +
  guides(color = guide_legend(ncol = 4, byrow = FALSE))

organism_palette <- setNames(
  hue_pal()(nrow(top_organisms)),
  top_organisms$organism_label
)

p_cumulative_top_organism_incidence <- ggplot(
  cumulative_top_organism_incidence_hour,
  aes(icu_hour, organism_first_collected_per_100_icu_admissions, color = organism_label)
) +
  geom_line(linewidth = 0.85) +
  scale_x_continuous(
    breaks = seq(0, timing_max_hour, by = 24),
    labels = function(x) x / 24
  ) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Cumulative Incidence of Top Organism Collection",
    x = "Days since ICU admission",
    y = "ICU admissions per 100",
    color = NULL
  ) +
  plot_theme +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, label.theme = element_text(size = 9))) +
  scale_color_manual(values = organism_palette, labels = label_wrap(42))

out_dir <- file.path("output", "icu_day_timing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

paths <- c(
  monthly_icu_days = file.path(out_dir, glue("monthly_icu_days_{site_name}_{stamp}.csv")),
  monthly_overall_rates = file.path(out_dir, glue("monthly_overall_culture_rates_per_100_icu_days_{site_name}_{stamp}.csv")),
  monthly_type_rates = file.path(out_dir, glue("monthly_specimen_type_culture_rates_per_100_icu_days_{site_name}_{stamp}.csv")),
  first_culture_timing_summary = file.path(out_dir, glue("first_culture_timing_summary_{site_name}_{stamp}.csv")),
  first_culture_timing_bins = file.path(out_dir, glue("first_culture_timing_bins_{site_name}_{stamp}.csv")),
  icu_day_event_rates = file.path(out_dir, glue("icu_day_culture_event_rates_{site_name}_{stamp}.csv")),
  cumulative_first_culture_by_day = file.path(out_dir, glue("cumulative_first_culture_by_icu_day_{site_name}_{stamp}.csv")),
  cumulative_culture_events_by_type_hour = file.path(out_dir, glue("cumulative_culture_events_by_type_icu_hour_{site_name}_{stamp}.csv")),
  cumulative_top_organism_incidence_hour = file.path(out_dir, glue("cumulative_top_organism_incidence_icu_hour_{site_name}_{stamp}.csv")),
  monthly_type_rate_plot = file.path(out_dir, glue("monthly_specimen_type_culture_rates_stacked_per_100_icu_days_{site_name}_{stamp}.png")),
  first_culture_timing_plot = file.path(out_dir, glue("first_culture_timing_by_specimen_type_{site_name}_{stamp}.png")),
  cumulative_first_culture_plot = file.path(out_dir, glue("cumulative_first_culture_by_icu_day_{site_name}_{stamp}.png")),
  icu_day_event_rate_plot = file.path(out_dir, glue("icu_day_culture_event_rates_per_100_icu_admissions_at_risk_{site_name}_{stamp}.png")),
  cumulative_events_by_type_plot = file.path(out_dir, glue("cumulative_culture_events_by_type_icu_hour_{site_name}_{stamp}.png")),
  cumulative_top_organism_incidence_plot = file.path(out_dir, glue("cumulative_top_organism_incidence_icu_hour_{site_name}_{stamp}.png"))
)

write_csv(monthly_icu_days, paths[["monthly_icu_days"]])
write_csv(monthly_overall_per_icu_day, paths[["monthly_overall_rates"]])
write_csv(monthly_events_by_type_per_icu_day, paths[["monthly_type_rates"]])
write_csv(first_culture_timing_summary, paths[["first_culture_timing_summary"]])
write_csv(timing_bin_summary, paths[["first_culture_timing_bins"]])
write_csv(icu_day_event_rates, paths[["icu_day_event_rates"]])
write_csv(cumulative_first_culture_by_day, paths[["cumulative_first_culture_by_day"]])
write_csv(cumulative_culture_events_by_type_hour, paths[["cumulative_culture_events_by_type_hour"]])
write_csv(cumulative_top_organism_incidence_hour, paths[["cumulative_top_organism_incidence_hour"]])

ggsave(paths[["monthly_type_rate_plot"]], p_icu_day_rate_stacked, width = 12, height = 7, dpi = 300)
ggsave(paths[["first_culture_timing_plot"]], p_timing_bins, width = 10, height = 7, dpi = 300)
ggsave(paths[["cumulative_first_culture_plot"]], p_cumulative_first_culture, width = 10, height = 6, dpi = 300)
ggsave(paths[["icu_day_event_rate_plot"]], p_icu_day_event_rates, width = 10, height = 6, dpi = 300)
ggsave(paths[["cumulative_events_by_type_plot"]], p_cumulative_events_by_type, width = 11, height = 7, dpi = 300)
ggsave(paths[["cumulative_top_organism_incidence_plot"]], p_cumulative_top_organism_incidence, width = 11, height = 7, dpi = 300)

message("ICU-day denominator summary:")
print(monthly_icu_days %>% summarise(
  first_month = min(calendar_month),
  last_month = max(calendar_month),
  total_icu_days = sum(n_icu_days),
  median_monthly_icu_days = median(n_icu_days),
  min_monthly_icu_days = min(n_icu_days),
  max_monthly_icu_days = max(n_icu_days)
), width = Inf)

message("")
message("Overall ICU-day culture rate summary:")
print(monthly_overall_per_icu_day %>% summarise(
  total_culture_events = sum(n_culture_events),
  median_monthly_events_per_100_icu_days = median(culture_events_per_100_icu_days, na.rm = TRUE),
  min_monthly_events_per_100_icu_days = min(culture_events_per_100_icu_days, na.rm = TRUE),
  max_monthly_events_per_100_icu_days = max(culture_events_per_100_icu_days, na.rm = TRUE)
), width = Inf)

message("")
message("First culture timing summary:")
print(first_culture_timing_summary, width = Inf)

message("")
message("Specimen types displayed separately:")
print(tibble(specimen_type = top_types), n = Inf)

message("")
message("Top organisms displayed separately:")
print(top_organisms, n = Inf)

message("")
message("Wrote outputs:")
print(paths)
