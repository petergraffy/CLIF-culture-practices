# ================================================================================================
# Organism Detection Trend Screen
#
# Question:
#   Which organisms or targeted resistance-related organism labels are increasing over calendar time?
#
# Denominator:
#   All ICU admissions, defined as merged ICU ADT intervals and counted by ICU admission month.
#
# Numerator:
#   Positive ICU culture detection events. Events are collapsed by ICU culture event and organism,
#   then counted monthly.
#
# Note:
#   True MRSA/VRE/CRE phenotypes require susceptibility/resistance fields. This script detects those
#   labels only when resistance terms are present in organism text, and otherwise reports organism
#   proxies such as Staphylococcus aureus.
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

latest_file <- function(pattern, path = file.path("data", "intermediate", "cohort")) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files found in ", path, " matching pattern: ", pattern)
  }
  files[which.max(file.info(files)$mtime)]
}

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

rolling_mean_trailing <- function(x, k = 6) {
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))
}

build_monthly_icu_admissions <- function(month_seq, study_start_dttm, study_end_dttm) {
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

  adt %>%
    filter(location_category == "icu", !is.na(icu_in_dttm)) %>%
    left_join(hospitalization, by = "hospitalization_id") %>%
    mutate(
      icu_out_dttm = coalesce(icu_out_dttm_raw, discharge_dttm),
      icu_interval_missing_out = is.na(icu_out_dttm_raw)
    ) %>%
    filter(!is.na(patient_id), !is.na(icu_out_dttm), icu_out_dttm > icu_in_dttm) %>%
    arrange(patient_id, hospitalization_id, icu_in_dttm, icu_out_dttm) %>%
    group_by(patient_id, hospitalization_id) %>%
    mutate(
      prior_max_icu_out_num = lag(cummax(as.numeric(icu_out_dttm))),
      new_icu_admission = is.na(prior_max_icu_out_num) | as.numeric(icu_in_dttm) > prior_max_icu_out_num,
      icu_admission_seq = cumsum(new_icu_admission)
    ) %>%
    group_by(patient_id, hospitalization_id, icu_admission_seq) %>%
    summarise(
      icu_in_dttm = min(icu_in_dttm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(calendar_month = floor_date(icu_in_dttm, "month")) %>%
    filter(is.na(study_start_dttm) | icu_in_dttm >= study_start_dttm) %>%
    filter(is.na(study_end_dttm) | icu_in_dttm <= study_end_dttm) %>%
    count(calendar_month, name = "n_icu_admissions") %>%
    complete(calendar_month = month_seq, fill = list(n_icu_admissions = 0L)) %>%
    arrange(calendar_month)
}

target_organism_patterns <- tibble::tribble(
  ~target_label, ~pattern,
  "MRSA text label", "mrsa|methicillin[ _-]*resistant.*staph|staph.*methicillin[ _-]*resistant",
  "Staphylococcus aureus", "staphylococcus[_ ]aureus",
  "VRE text label", "\\bvre\\b|vancomycin[ _-]*resistant.*enterococcus|enterococcus.*vancomycin[ _-]*resistant",
  "Enterococcus faecium", "enterococcus[_ ]faecium",
  "ESBL text label", "\\besbl\\b|extended[ _-]*spectrum",
  "CRE text label", "\\bcre\\b|carbapenem[ _-]*resistant|\\bkpc\\b|\\bndm\\b",
  "Pseudomonas aeruginosa", "pseudomonas[_ ]aeruginosa",
  "Klebsiella pneumoniae", "klebsiella[_ ]pneumoniae",
  "Escherichia coli", "escherichia[_ ]coli|\\be[._ ]?coli\\b",
  "Candida auris", "candida[_ ]auris",
  "Clostridioides difficile", "clostridioides[_ ]difficile|clostridium[_ ]difficile"
)

fit_poisson_trend <- function(data) {
  model_data <- data %>%
    filter(!is.na(n_icu_admissions), n_icu_admissions > 0) %>%
    mutate(month_index = as.numeric(difftime(calendar_month, min(calendar_month), units = "days")) / 30.4375)

  if (nrow(model_data) < 12 || sum(model_data$n_detection_events, na.rm = TRUE) == 0) {
    return(tibble(
      monthly_irr = NA_real_,
      annual_irr = NA_real_,
      annual_percent_change = NA_real_,
      p_value = NA_real_
    ))
  }

  fit <- glm(
    n_detection_events ~ month_index,
    family = poisson(),
    offset = log(n_icu_admissions),
    data = model_data
  )
  beta <- unname(coef(fit)[["month_index"]])
  p_value <- unname(summary(fit)$coefficients["month_index", "Pr(>|z|)"])

  tibble(
    monthly_irr = exp(beta),
    annual_irr = exp(beta * 12),
    annual_percent_change = 100 * (annual_irr - 1),
    p_value = p_value
  )
}

make_monthly_rates <- function(data, label_var, month_seq, monthly_icu_admissions) {
  data %>%
    group_by(calendar_month, organism_label = .data[[label_var]]) %>%
    summarise(n_detection_events = n_distinct(detection_event_id), .groups = "drop") %>%
    complete(
      calendar_month = month_seq,
      organism_label,
      fill = list(n_detection_events = 0L)
    ) %>%
    left_join(monthly_icu_admissions, by = "calendar_month") %>%
    mutate(
      detection_events_per_100_icu_admissions = if_else(
        n_icu_admissions > 0,
        100 * n_detection_events / n_icu_admissions,
        NA_real_
      )
    ) %>%
    arrange(organism_label, calendar_month)
}

theme_trends <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = grid::unit(3, "pt"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

site_name <- clif_site_name
row_path <- Sys.getenv("ICU_CULTURE_ROWS_PATH", unset = NA_character_)
study_start_date <- Sys.getenv("STUDY_START_DATE", unset = Sys.getenv("PLOT_START_DATE", unset = NA_character_))
study_end_date <- Sys.getenv("STUDY_END_DATE", unset = Sys.getenv("PLOT_END_DATE", unset = NA_character_))
top_n_organisms <- as.integer(Sys.getenv("TOP_N_TREND_ORGANISMS", unset = "25"))
plot_n_increasing <- as.integer(Sys.getenv("PLOT_N_INCREASING_ORGANISMS", unset = "12"))

if (is.na(row_path) || !nzchar(row_path)) {
  row_path <- latest_file("^icu_culture_rows_.*\\.csv$")
}

study_start_dttm <- if (!is.na(study_start_date) && nzchar(study_start_date)) safe_ts(study_start_date) else as.POSIXct(NA)
study_end_dttm <- if (!is.na(study_end_date) && nzchar(study_end_date)) safe_ts(study_end_date) + days(1) - seconds(1) else as.POSIXct(NA)

message("Reading ICU culture rows: ", row_path)
if (!is.na(study_start_dttm)) message("Study start: ", study_start_dttm)
if (!is.na(study_end_dttm)) message("Study end: ", study_end_dttm)

rows <- readr::read_csv(row_path, show_col_types = FALSE) %>%
  mutate(
    collect_dttm = safe_ts(collect_dttm),
    calendar_month = floor_date(collect_dttm, "month"),
    fluid_name = coalesce(na_if(fluid_name, ""), "missing"),
    method_name = coalesce(na_if(method_name, ""), "missing"),
    organism_group = coalesce(na_if(str_to_lower(str_trim(as.character(organism_group))), ""), "missing"),
    organism_category = coalesce(na_if(str_to_lower(str_trim(as.character(organism_category))), ""), organism_group),
    organism_name = coalesce(na_if(str_to_lower(str_trim(as.character(organism_name))), ""), organism_category),
    organism_category_label = clean_label(organism_category),
    organism_text = str_squish(str_c(organism_name, organism_category, organism_group, sep = " ")),
    positive_culture = as.logical(positive_culture),
    explicit_negative_name = drop_negative_organism_name(organism_name),
    detection_event_id = str_c(
      patient_id,
      hospitalization_id,
      icu_interval_id,
      collect_dttm,
      fluid_name,
      method_name,
      organism_category,
      sep = "|"
    )
  ) %>%
  filter(positive_culture, !explicit_negative_name, !is.na(calendar_month)) %>%
  filter(is.na(study_start_dttm) | collect_dttm >= study_start_dttm) %>%
  filter(is.na(study_end_dttm) | collect_dttm <= study_end_dttm)

if (nrow(rows) == 0) {
  stop("No positive organism rows after filters.")
}

month_seq <- seq(min(rows$calendar_month), max(rows$calendar_month), by = "month")
monthly_icu_admissions <- build_monthly_icu_admissions(month_seq, study_start_dttm, study_end_dttm)

top_organism_labels <- rows %>%
  distinct(detection_event_id, organism_category_label) %>%
  count(organism_category_label, name = "total_detection_events", sort = TRUE) %>%
  slice_head(n = top_n_organisms)

top_monthly_rates <- rows %>%
  semi_join(top_organism_labels, by = "organism_category_label") %>%
  make_monthly_rates("organism_category_label", month_seq, monthly_icu_admissions) %>%
  left_join(top_organism_labels, by = c("organism_label" = "organism_category_label"))

target_detections <- target_organism_patterns %>%
  tidyr::crossing(row_id = seq_len(nrow(rows))) %>%
  mutate(row_match = str_detect(rows$organism_text[row_id], regex(pattern, ignore_case = TRUE))) %>%
  filter(row_match) %>%
  transmute(
    target_label,
    detection_event_id = rows$detection_event_id[row_id],
    calendar_month = rows$calendar_month[row_id]
  ) %>%
  distinct()

target_monthly_rates <- if (nrow(target_detections) > 0) {
  make_monthly_rates(target_detections, "target_label", month_seq, monthly_icu_admissions) %>%
    mutate(total_detection_events = sum(n_detection_events), .by = organism_label)
} else {
  tidyr::expand_grid(calendar_month = month_seq, organism_label = target_organism_patterns$target_label) %>%
    left_join(monthly_icu_admissions, by = "calendar_month") %>%
    mutate(
      n_detection_events = 0L,
      detection_events_per_100_icu_admissions = if_else(n_icu_admissions > 0, 0, NA_real_),
      total_detection_events = 0L
    )
}

trend_summary_top <- top_monthly_rates %>%
  group_by(organism_label) %>%
  group_modify(~ fit_poisson_trend(.x)) %>%
  ungroup() %>%
  left_join(
    top_monthly_rates %>%
      group_by(organism_label) %>%
      summarise(
        total_detection_events = max(total_detection_events, na.rm = TRUE),
        first_nonzero_month = suppressWarnings(min(calendar_month[n_detection_events > 0], na.rm = TRUE)),
        last_nonzero_month = suppressWarnings(max(calendar_month[n_detection_events > 0], na.rm = TRUE)),
        mean_monthly_rate_per_100_icu_admissions = mean(detection_events_per_100_icu_admissions, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "organism_label"
  ) %>%
  mutate(
    trend_direction = case_when(
      is.na(annual_percent_change) | is.na(p_value) ~ "Not estimated",
      p_value < 0.05 & annual_percent_change > 0 ~ "Increasing",
      p_value < 0.05 & annual_percent_change < 0 ~ "Decreasing",
      TRUE ~ "No clear trend"
    )
  ) %>%
  arrange(desc(annual_percent_change))

trend_summary_targets <- target_monthly_rates %>%
  group_by(organism_label) %>%
  group_modify(~ fit_poisson_trend(.x)) %>%
  ungroup() %>%
  left_join(
    target_monthly_rates %>%
      group_by(organism_label) %>%
      summarise(
        total_detection_events = max(total_detection_events, na.rm = TRUE),
        first_nonzero_month = suppressWarnings(min(calendar_month[n_detection_events > 0], na.rm = TRUE)),
        last_nonzero_month = suppressWarnings(max(calendar_month[n_detection_events > 0], na.rm = TRUE)),
        mean_monthly_rate_per_100_icu_admissions = mean(detection_events_per_100_icu_admissions, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "organism_label"
  ) %>%
  mutate(
    first_nonzero_month = if_else(is.infinite(first_nonzero_month), as.POSIXct(NA), first_nonzero_month),
    last_nonzero_month = if_else(is.infinite(last_nonzero_month), as.POSIXct(NA), last_nonzero_month),
    trend_direction = case_when(
      total_detection_events == 0 ~ "Not observed",
      is.na(annual_percent_change) | is.na(p_value) ~ "Not estimated",
      p_value < 0.05 & annual_percent_change > 0 ~ "Increasing",
      p_value < 0.05 & annual_percent_change < 0 ~ "Decreasing",
      TRUE ~ "No clear trend"
    )
  ) %>%
  arrange(desc(annual_percent_change))

increasing_labels <- trend_summary_top %>%
  filter(total_detection_events >= 25, !is.na(annual_percent_change)) %>%
  arrange(desc(annual_percent_change)) %>%
  slice_head(n = plot_n_increasing) %>%
  pull(organism_label)

plot_increasing_data <- top_monthly_rates %>%
  filter(organism_label %in% increasing_labels) %>%
  arrange(organism_label, calendar_month) %>%
  group_by(organism_label) %>%
  mutate(rate_rolling_6mo = rolling_mean_trailing(detection_events_per_100_icu_admissions, 6)) %>%
  ungroup() %>%
  mutate(organism_label = factor(organism_label, levels = increasing_labels))

target_plot_labels <- target_monthly_rates %>%
  group_by(organism_label) %>%
  summarise(total_detection_events = max(total_detection_events, na.rm = TRUE), .groups = "drop") %>%
  filter(total_detection_events > 0 | organism_label %in% c("MRSA text label", "Staphylococcus aureus")) %>%
  arrange(desc(total_detection_events), organism_label) %>%
  pull(organism_label)

plot_target_data <- target_monthly_rates %>%
  filter(organism_label %in% target_plot_labels) %>%
  arrange(organism_label, calendar_month) %>%
  group_by(organism_label) %>%
  mutate(rate_rolling_6mo = rolling_mean_trailing(detection_events_per_100_icu_admissions, 6)) %>%
  ungroup() %>%
  mutate(organism_label = factor(organism_label, levels = target_plot_labels))

out_dir <- file.path("output", "organism_trends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

paths <- c(
  trend_summary_top = file.path(out_dir, glue("organism_detection_trend_screen_top_{site_name}_{stamp}.csv")),
  trend_summary_targets = file.path(out_dir, glue("organism_detection_trend_screen_targets_{site_name}_{stamp}.csv")),
  monthly_top_rates = file.path(out_dir, glue("monthly_top_organism_detection_rates_per_100_icu_admissions_{site_name}_{stamp}.csv")),
  monthly_target_rates = file.path(out_dir, glue("monthly_target_organism_detection_rates_per_100_icu_admissions_{site_name}_{stamp}.csv")),
  monthly_icu_admissions = file.path(out_dir, glue("monthly_icu_admissions_for_organism_trends_{site_name}_{stamp}.csv")),
  increasing_plot = file.path(out_dir, glue("monthly_fastest_increasing_organism_detection_rates_{site_name}_{stamp}.png")),
  target_plot = file.path(out_dir, glue("monthly_target_organism_detection_rates_{site_name}_{stamp}.png"))
)

write_csv(trend_summary_top, paths[["trend_summary_top"]])
write_csv(trend_summary_targets, paths[["trend_summary_targets"]])
write_csv(top_monthly_rates, paths[["monthly_top_rates"]])
write_csv(target_monthly_rates, paths[["monthly_target_rates"]])
write_csv(monthly_icu_admissions, paths[["monthly_icu_admissions"]])

p_increasing <- ggplot(plot_increasing_data, aes(calendar_month, detection_events_per_100_icu_admissions)) +
  geom_col(fill = "#4C78A8", width = 25 * 24 * 60 * 60) +
  geom_line(aes(y = rate_rolling_6mo), color = "black", linewidth = 0.65, na.rm = TRUE) +
  facet_wrap(vars(organism_label), scales = "free_y", ncol = 3) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Fastest Increasing Organism Detection Rates",
    x = NULL,
    y = "Detection events per 100 ICU admissions"
  ) +
  theme_trends

p_targets <- ggplot(plot_target_data, aes(calendar_month, detection_events_per_100_icu_admissions)) +
  geom_col(fill = "#D55E00", width = 25 * 24 * 60 * 60) +
  geom_line(aes(y = rate_rolling_6mo), color = "black", linewidth = 0.65, na.rm = TRUE) +
  facet_wrap(vars(organism_label), scales = "free_y", ncol = 2) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(
    title = "Target Organism Detection Rates",
    x = NULL,
    y = "Detection events per 100 ICU admissions"
  ) +
  theme_trends

ggsave(paths[["increasing_plot"]], p_increasing, width = 14, height = 12, dpi = 300)
ggsave(paths[["target_plot"]], p_targets, width = 12, height = 10, dpi = 300)

message("Wrote trend summaries:")
print(paths[c("trend_summary_top", "trend_summary_targets", "monthly_top_rates", "monthly_target_rates")])
message("Wrote plots:")
print(paths[c("increasing_plot", "target_plot")])

message("Top increasing organisms by annual percent change:")
print(
  trend_summary_top %>%
    select(organism_label, total_detection_events, annual_percent_change, p_value, trend_direction) %>%
    head(15),
  n = 15,
  width = Inf
)

message("Target organism trend summary:")
print(
  trend_summary_targets %>%
    select(organism_label, total_detection_events, annual_percent_change, p_value, trend_direction) %>%
    arrange(desc(total_detection_events)),
  n = Inf,
  width = Inf
)
