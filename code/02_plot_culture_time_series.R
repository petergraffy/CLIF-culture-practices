# ================================================================================================
# Plot ICU Culture Time Series
#
# Input:
#   ICU culture event output from code/01_identify_icu_culture_cohort.R.
#
# Outputs:
#   Monthly culture volume and positivity summaries plus PNG figures.
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

latest_file <- function(pattern, path = file.path("data", "intermediate", "cohort")) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files found in ", path, " matching pattern: ", pattern)
  }
  files[which.max(file.info(files)$mtime)]
}

safe_ts <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = tz))
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

site_name <- Sys.getenv("CLIF_SITE_NAME", unset = "SITE")
event_path <- Sys.getenv("ICU_CULTURE_EVENTS_PATH", unset = NA_character_)
top_n_types <- as.integer(Sys.getenv("TOP_N_CULTURE_TYPES", unset = "8"))
plot_start_date <- Sys.getenv("PLOT_START_DATE", unset = NA_character_)
plot_end_date <- Sys.getenv("PLOT_END_DATE", unset = NA_character_)

if (is.na(event_path) || !nzchar(event_path)) {
  event_path <- latest_file("^icu_culture_events_.*\\.csv$")
}

message("Reading ICU culture events: ", event_path)

events <- readr::read_csv(event_path, show_col_types = FALSE) %>%
  mutate(
    collect_dttm = safe_ts(collect_dttm),
    culture_month = floor_date(collect_dttm, "month"),
    fluid_category = coalesce(na_if(fluid_category, ""), "missing"),
    culture_type = clean_label(fluid_category),
    any_positive_culture = as.logical(any_positive_culture)
  ) %>%
  filter(!is.na(culture_month))

plot_start_dttm <- if (!is.na(plot_start_date) && nzchar(plot_start_date)) safe_ts(plot_start_date) else as.POSIXct(NA)
plot_end_dttm <- if (!is.na(plot_end_date) && nzchar(plot_end_date)) safe_ts(plot_end_date) + days(1) - seconds(1) else as.POSIXct(NA)

events <- events %>%
  filter(is.na(plot_start_dttm) | collect_dttm >= plot_start_dttm) %>%
  filter(is.na(plot_end_dttm) | collect_dttm <= plot_end_dttm)

if (nrow(events) == 0) {
  stop("No events with non-missing collection month.")
}

if (!is.na(plot_start_dttm)) message("Plot start: ", plot_start_dttm)
if (!is.na(plot_end_dttm)) message("Plot end: ", plot_end_dttm)

top_types <- events %>%
  count(culture_type, sort = TRUE) %>%
  slice_head(n = top_n_types) %>%
  pull(culture_type)

events_plot <- events %>%
  mutate(culture_type_plot = if_else(culture_type %in% top_types, culture_type, "Other")) %>%
  mutate(culture_type_plot = fct_reorder(culture_type_plot, culture_type_plot == "Other", .desc = TRUE))

monthly_overall <- events %>%
  group_by(culture_month) %>%
  summarise(
    n_events = n(),
    n_positive_events = sum(any_positive_culture, na.rm = TRUE),
    positive_event_rate = n_positive_events / n_events,
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  )

monthly_by_type <- events_plot %>%
  group_by(culture_month, culture_type = culture_type_plot) %>%
  summarise(
    n_events = n(),
    n_positive_events = sum(any_positive_culture, na.rm = TRUE),
    positive_event_rate = n_positive_events / n_events,
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  complete(
    culture_month = seq(min(events$culture_month), max(events$culture_month), by = "month"),
    culture_type,
    fill = list(
      n_events = 0L,
      n_positive_events = 0L,
      n_hospitalizations = 0L,
      n_patients = 0L
    )
  ) %>%
  mutate(
    positive_event_rate = if_else(n_events > 0, n_positive_events / n_events, NA_real_),
    culture_type = fct_relevel(factor(culture_type), "Other", after = Inf)
  )

out_dir <- file.path("output", "time_series")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

overall_path <- file.path(out_dir, glue("monthly_overall_culture_events_{site_name}_{stamp}.csv"))
type_path <- file.path(out_dir, glue("monthly_culture_events_by_type_{site_name}_{stamp}.csv"))
readr::write_csv(monthly_overall, overall_path)
readr::write_csv(monthly_by_type, type_path)

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
month_bar_width <- 25 * 24 * 60 * 60
culture_type_palette <- c(
  "Blood buffy" = "#4C78A8",
  "Respiratory tract" = "#F58518",
  "Genito urinary tract" = "#54A24B",
  "Meninges csf" = "#B279A2",
  "Other unspecified" = "#72B7B2",
  "Pleural cavity fluid" = "#E45756",
  "Respiratory tract lower" = "#EECA3B",
  "Woundsite" = "#FF9DA6",
  "Other" = "#9D9D9D"
)
available_palette <- culture_type_palette[names(culture_type_palette) %in% levels(monthly_by_type$culture_type)]

monthly_by_type_stacked <- monthly_by_type %>%
  mutate(culture_type = fct_collapse(culture_type, Other = c("Other", "Other unspecified"))) %>%
  group_by(culture_month, culture_type) %>%
  summarise(
    n_events = sum(n_events),
    n_positive_events = sum(n_positive_events),
    positive_event_rate = if_else(n_events > 0, n_positive_events / n_events, NA_real_),
    n_hospitalizations = sum(n_hospitalizations),
    n_patients = sum(n_patients),
    .groups = "drop"
  ) %>%
  mutate(culture_type = fct_relevel(factor(culture_type), "Other", after = Inf))

stacked_palette <- culture_type_palette[names(culture_type_palette) %in% levels(monthly_by_type_stacked$culture_type)]

p_overall_volume <- ggplot(monthly_overall, aes(culture_month, n_events)) +
  geom_col(fill = "#2f6f73", width = month_bar_width) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Monthly ICU Culture Events",
    x = NULL,
    y = "Culture events"
  ) +
  plot_theme

p_overall_positivity <- ggplot(monthly_overall, aes(culture_month, positive_event_rate)) +
  geom_line(color = "#8f3d56", linewidth = 0.8) +
  geom_point(color = "#8f3d56", size = 1.2) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Monthly ICU Culture Event Positivity",
    x = NULL,
    y = "Positive events"
  ) +
  plot_theme

p_type_volume <- ggplot(monthly_by_type, aes(culture_month, n_events, color = culture_type)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = available_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Monthly ICU Culture Events by Culture Type",
    x = NULL,
    y = "Culture events",
    color = NULL
  ) +
  plot_theme

p_type_stacked_volume <- ggplot(monthly_by_type_stacked, aes(culture_month, n_events, fill = culture_type)) +
  geom_col(width = month_bar_width, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = stacked_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Monthly ICU Culture Events by Culture Type",
    subtitle = "Stacked monthly event counts; Other includes Other unspecified and less common culture types",
    x = NULL,
    y = "Culture events",
    fill = NULL
  ) +
  plot_theme

p_type_positivity <- monthly_by_type %>%
  filter(n_events >= 10) %>%
  ggplot(aes(culture_month, positive_event_rate, color = culture_type)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = available_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Monthly ICU Culture Event Positivity by Culture Type",
    subtitle = "Months with at least 10 events for the displayed culture type",
    x = NULL,
    y = "Positive events",
    color = NULL
  ) +
  plot_theme

p_type_facets <- monthly_by_type %>%
  filter(culture_type != "Other") %>%
  ggplot(aes(culture_month, n_events)) +
  geom_col(fill = "#557a46", width = month_bar_width) +
  facet_wrap(vars(culture_type), scales = "free_y", ncol = 2) +
  scale_x_datetime(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Monthly ICU Culture Events for Major Culture Types",
    x = NULL,
    y = "Culture events"
  ) +
  plot_theme +
  theme(legend.position = "none")

plot_paths <- c(
  overall_volume = file.path(out_dir, glue("monthly_overall_culture_volume_{site_name}_{stamp}.png")),
  overall_positivity = file.path(out_dir, glue("monthly_overall_culture_positivity_{site_name}_{stamp}.png")),
  type_volume = file.path(out_dir, glue("monthly_culture_volume_by_type_{site_name}_{stamp}.png")),
  type_stacked_volume = file.path(out_dir, glue("monthly_culture_volume_stacked_by_type_{site_name}_{stamp}.png")),
  type_positivity = file.path(out_dir, glue("monthly_culture_positivity_by_type_{site_name}_{stamp}.png")),
  type_facets = file.path(out_dir, glue("monthly_culture_volume_major_type_facets_{site_name}_{stamp}.png"))
)

ggsave(plot_paths[["overall_volume"]], p_overall_volume, width = 9, height = 5, dpi = 300)
ggsave(plot_paths[["overall_positivity"]], p_overall_positivity, width = 9, height = 5, dpi = 300)
ggsave(plot_paths[["type_volume"]], p_type_volume, width = 11, height = 6, dpi = 300)
ggsave(plot_paths[["type_stacked_volume"]], p_type_stacked_volume, width = 11, height = 6.5, dpi = 300)
ggsave(plot_paths[["type_positivity"]], p_type_positivity, width = 11, height = 6, dpi = 300)
ggsave(plot_paths[["type_facets"]], p_type_facets, width = 11, height = 9, dpi = 300)

message("Monthly overall summary:")
print(monthly_overall %>% summarise(
  first_month = min(culture_month),
  last_month = max(culture_month),
  total_events = sum(n_events),
  total_positive_events = sum(n_positive_events),
  event_positivity = total_positive_events / total_events,
  median_monthly_events = median(n_events),
  min_monthly_events = min(n_events),
  max_monthly_events = max(n_events)
), width = Inf)

message("")
message("Culture types displayed separately:")
print(tibble(culture_type = top_types), n = Inf)
message("")
message("Wrote monthly overall summary: ", overall_path)
message("Wrote monthly type summary: ", type_path)
message("Wrote plots:")
print(plot_paths)
