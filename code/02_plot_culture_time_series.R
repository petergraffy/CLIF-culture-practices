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

collapse_culture_panel <- function(x) {
  case_when(
    x == "Blood buffy" ~ "Blood buffy",
    x == "Genito urinary tract" ~ "Genito urinary tract",
    x == "Respiratory tract" ~ "Respiratory tract",
    TRUE ~ "Other"
  )
}

five_panel_levels <- c(
  "Overall",
  "Blood buffy",
  "Genito urinary tract",
  "Respiratory tract",
  "Other"
)

monthly_event_density <- function(data, panel_var, count_var = "n_events") {
  data %>%
    group_by(.data[[panel_var]]) %>%
    mutate(
      total_events = sum(.data[[count_var]], na.rm = TRUE),
      monthly_total_event_density = if_else(total_events > 0, .data[[count_var]] / total_events, NA_real_)
    ) %>%
    ungroup() %>%
    transmute(
      culture_panel = .data[[panel_var]],
      culture_month,
      n_events = .data[[count_var]],
      total_events,
      monthly_total_event_density
    )
}

site_name <- Sys.getenv("CLIF_SITE_NAME", unset = "SITE")
event_path <- Sys.getenv("ICU_CULTURE_EVENTS_PATH", unset = NA_character_)
top_n_types <- as.integer(Sys.getenv("TOP_N_CULTURE_TYPES", unset = "8"))
top_n_positivity_types <- as.integer(Sys.getenv("TOP_N_POSITIVITY_FLUID_CATEGORIES", unset = as.character(top_n_types)))
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

month_seq <- seq(min(events$culture_month), max(events$culture_month), by = "month")

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

positivity_fluid_categories <- events %>%
  count(fluid_category, culture_type, sort = TRUE) %>%
  slice_head(n = top_n_positivity_types)

monthly_positivity_by_fluid_category <- events %>%
  inner_join(
    positivity_fluid_categories %>% select(fluid_category, culture_type),
    by = c("fluid_category", "culture_type")
  ) %>%
  group_by(culture_month, fluid_category, fluid_category_label = culture_type) %>%
  summarise(
    n_events = n(),
    n_positive_events = sum(any_positive_culture, na.rm = TRUE),
    positive_event_rate = n_positive_events / n_events,
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  complete(
    culture_month = month_seq,
    nesting(fluid_category, fluid_category_label),
    fill = list(
      n_events = 0L,
      n_positive_events = 0L,
      n_hospitalizations = 0L,
      n_patients = 0L
    )
  ) %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    positive_event_rate = if_else(n_events > 0, n_positive_events / n_events, NA_real_),
    fluid_category_label = fct_reorder(fluid_category_label, n_events, .fun = sum, .desc = TRUE)
  )

monthly_positivity_five_panel <- bind_rows(
  events %>% mutate(culture_panel = "Overall"),
  events %>% mutate(culture_panel = collapse_culture_panel(culture_type))
) %>%
  mutate(culture_panel = factor(culture_panel, levels = five_panel_levels)) %>%
  group_by(culture_month, culture_panel) %>%
  summarise(
    n_events = n(),
    n_positive_events = sum(any_positive_culture, na.rm = TRUE),
    positive_event_rate = n_positive_events / n_events,
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  complete(
    culture_month = month_seq,
    culture_panel = factor(five_panel_levels, levels = five_panel_levels),
    fill = list(
      n_events = 0L,
      n_positive_events = 0L,
      n_hospitalizations = 0L,
      n_patients = 0L
    )
  ) %>%
  mutate(
    site_name = site_name,
    care_setting = "ICU",
    positive_event_rate = if_else(n_events > 0, n_positive_events / n_events, NA_real_)
  )

five_panel_total_event_density <- monthly_event_density(monthly_positivity_five_panel, "culture_panel")
density_plot_scale <- 1 / max(five_panel_total_event_density$monthly_total_event_density, na.rm = TRUE)

out_dir <- file.path("output", "time_series")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

overall_path <- file.path(out_dir, glue("monthly_overall_culture_events_{site_name}_{stamp}.csv"))
type_path <- file.path(out_dir, glue("monthly_culture_events_by_type_{site_name}_{stamp}.csv"))
fluid_positivity_path <- file.path(out_dir, glue("monthly_positive_culture_event_rate_by_fluid_category_{site_name}_{stamp}.csv"))
five_panel_positivity_path <- file.path(out_dir, glue("monthly_positive_culture_event_rate_five_panel_{site_name}_{stamp}.csv"))
five_panel_density_path <- file.path(out_dir, glue("monthly_total_culture_event_density_five_panel_{site_name}_{stamp}.csv"))
readr::write_csv(monthly_overall, overall_path)
readr::write_csv(monthly_by_type, type_path)
readr::write_csv(monthly_positivity_by_fluid_category, fluid_positivity_path)
readr::write_csv(monthly_positivity_five_panel, five_panel_positivity_path)
readr::write_csv(five_panel_total_event_density, five_panel_density_path)

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
  "Blood buffy" = "#BC80BD",
  "Respiratory tract" = "#8DD3C7",
  "Genito urinary tract" = "#BEBADA",
  "Meninges csf" = "#FB8072",
  "Other unspecified" = "#80B1D3",
  "Pleural cavity fluid" = "#FDB462",
  "Respiratory tract lower" = "#B3DE69",
  "Woundsite" = "#FCCDE5",
  "Other" = "#BDBDBD"
)
five_panel_palette <- c(
  "Overall" = "#333333",
  "Blood buffy" = culture_type_palette[["Blood buffy"]],
  "Genito urinary tract" = culture_type_palette[["Genito urinary tract"]],
  "Respiratory tract" = culture_type_palette[["Respiratory tract"]],
  "Other" = culture_type_palette[["Other"]]
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
    x = NULL,
    y = "Positive events",
    color = NULL
  ) +
  plot_theme

p_fluid_category_positivity_facets <- monthly_positivity_by_fluid_category %>%
  filter(n_events > 0) %>%
  ggplot(aes(culture_month)) +
  geom_col(aes(y = n_events), fill = "#D9D9D9", width = month_bar_width, color = NA) +
  geom_col(
    aes(y = n_positive_events, fill = fluid_category_label),
    width = month_bar_width * 0.82,
    color = "white",
    linewidth = 0.05
  ) +
  facet_wrap(vars(fluid_category_label), ncol = 2) +
  scale_x_datetime(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  scale_fill_manual(values = culture_type_palette) +
  labs(
    title = "Monthly Total (Gray) and Positive Culture Events by Fluid Category",
    x = NULL,
    y = "Culture events",
    fill = NULL
  ) +
  plot_theme +
  guides(fill = guide_legend(ncol = 4, byrow = TRUE))

fluid_category_positivity_levels <- levels(monthly_positivity_by_fluid_category$fluid_category_label)
fluid_category_positivity_palette <- culture_type_palette[
  names(culture_type_palette) %in% fluid_category_positivity_levels
]

p_fluid_category_positivity_overlay <- monthly_positivity_by_fluid_category %>%
  filter(n_events > 0) %>%
  ggplot(aes(culture_month, positive_event_rate, color = fluid_category_label)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = fluid_category_positivity_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Monthly Positive Culture Event Rate by Fluid Category",
    x = NULL,
    y = "Positive culture events",
    color = NULL
  ) +
  plot_theme +
  guides(color = guide_legend(ncol = 4, byrow = TRUE))

p_five_panel_positive_rate_density <- ggplot() +
  geom_col(
    data = monthly_positivity_five_panel %>% filter(n_events > 0),
    aes(culture_month, positive_event_rate, fill = culture_panel),
    width = month_bar_width,
    color = "white",
    linewidth = 0.05
  ) +
  geom_area(
    data = five_panel_total_event_density,
    aes(culture_month, monthly_total_event_density * density_plot_scale),
    fill = "#7A7A7A",
    alpha = 0.10
  ) +
  facet_grid(rows = vars(culture_panel)) +
  scale_fill_manual(values = five_panel_palette) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    sec.axis = sec_axis(
      ~ . / density_plot_scale,
      labels = percent_format(accuracy = 0.1),
      name = "Monthly total event density"
    )
  ) +
  labs(
    title = "Monthly Positive Culture Event Rate and Total Event Density",
    x = NULL,
    y = "Positive culture event rate",
    fill = NULL
  ) +
  plot_theme +
  guides(fill = guide_legend(nrow = 1))

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
  fluid_category_positivity_facets = file.path(out_dir, glue("monthly_positive_culture_event_rate_by_fluid_category_{site_name}_{stamp}.png")),
  fluid_category_positivity_overlay = file.path(out_dir, glue("monthly_positive_culture_event_rate_by_fluid_category_overlay_{site_name}_{stamp}.png")),
  five_panel_positive_rate_density = file.path(out_dir, glue("monthly_positive_culture_event_rate_density_five_panel_{site_name}_{stamp}.png")),
  type_facets = file.path(out_dir, glue("monthly_culture_volume_major_type_facets_{site_name}_{stamp}.png"))
)

ggsave(plot_paths[["overall_volume"]], p_overall_volume, width = 9, height = 5, dpi = 300)
ggsave(plot_paths[["overall_positivity"]], p_overall_positivity, width = 9, height = 5, dpi = 300)
ggsave(plot_paths[["type_volume"]], p_type_volume, width = 11, height = 6, dpi = 300)
ggsave(plot_paths[["type_stacked_volume"]], p_type_stacked_volume, width = 11, height = 6.5, dpi = 300)
ggsave(plot_paths[["type_positivity"]], p_type_positivity, width = 11, height = 6, dpi = 300)
ggsave(plot_paths[["fluid_category_positivity_facets"]], p_fluid_category_positivity_facets, width = 12, height = 12, dpi = 300)
ggsave(plot_paths[["fluid_category_positivity_overlay"]], p_fluid_category_positivity_overlay, width = 11, height = 6.5, dpi = 300)
ggsave(plot_paths[["five_panel_positive_rate_density"]], p_five_panel_positive_rate_density, width = 12, height = 14, dpi = 300)
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
message("Wrote monthly fluid-category positivity summary: ", fluid_positivity_path)
message("Wrote monthly five-panel positivity summary: ", five_panel_positivity_path)
message("Wrote monthly five-panel total event density summary: ", five_panel_density_path)
message("Wrote plots:")
print(plot_paths)
