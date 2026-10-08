# ================================================================================================
# Organism Detection Trend Screen
#
# Question:
#   Which organisms or targeted resistance-related organism labels are increasing or decreasing over calendar time?
#
# Denominators:
#   1. All ICU admissions, defined as merged ICU ADT intervals and counted by ICU admission month.
#   2. All ICU days, allocated to calendar months from merged ICU ADT interval overlap time.
#
# Numerator:
#   Positive ICU culture detection events. Events are collapsed by ICU culture event and organism,
#   then counted monthly.
#
# Note: text-reported resistance is a separate screen. Standardized AST analysis is in script 09.
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
source("utils/trends.R")


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

classify_organism_type <- classify_microbe_taxonomy

organism_type_levels <- c(
  "Gram positive bacteria",
  "Gram negative bacteria",
  "Fungi/yeast",
  "Mycobacteria/AFB",
  "Anaerobes",
  "Atypical/other bacteria",
  "Other bacteria",
  "Virus",
  "Parasite/protozoa",
  "Other/unspecified"
)

organism_type_base_colors <- c(
  "Gram positive bacteria" = "#2F6C99",
  "Gram negative bacteria" = "#D55E00",
  "Fungi/yeast" = "#7B4AB8",
  "Mycobacteria/AFB" = "#008B8B",
  "Anaerobes" = "#8C6D31",
  "Atypical/other bacteria" = "#6A994E",
  "Other bacteria" = "#4E8F4A",
  "Virus" = "#C44E52",
  "Parasite/protozoa" = "#B56576",
  "Other/unspecified" = "#6B7280"
)

build_monthly_icu_denominators <- function(month_seq, study_start_dttm, study_end_dttm) {
  data <- read_culture_data(study_start_dttm, study_end_dttm)
  monthly_icu_denominators(data$icu_admissions, month_seq) %>%
    left_join(data$events %>% count(calendar_month = floor_date(collect_dttm, "month"), name = "n_observed_culture_events"), by = "calendar_month") %>%
    mutate(n_observed_culture_events = coalesce(n_observed_culture_events, 0L))

}

target_organism_patterns <- tibble::tribble(
  ~target_label, ~pattern,
  "MRSA (organism text)", "mrsa|methicillin[ _-]*resistant.*staph|staph.*methicillin[ _-]*resistant",
  "Staphylococcus aureus", "staphylococcus[_ ]aureus",
  "VRE (organism text)", "\\bvre\\b|vancomycin[ _-]*resistant.*enterococcus|enterococcus.*vancomycin[ _-]*resistant",
  "Enterococcus faecium", "enterococcus[_ ]faecium",
  "ESBL (organism text)", "\\besbl\\b|extended[ _-]*spectrum",
  "CRE (organism text)", "\\bcre\\b|carbapenem[ _-]*resistant|\\bkpc\\b|\\bndm\\b",
  "Pseudomonas aeruginosa", "pseudomonas[_ ]aeruginosa",
  "Klebsiella pneumoniae", "klebsiella[_ ]pneumoniae",
  "Escherichia coli", "escherichia[_ ]coli|\\be[._ ]?coli\\b",
  "Candida auris", "candida[_ ]auris",
  "Clostridioides difficile", "clostridioides[_ ]difficile|clostridium[_ ]difficile"
) %>%
  mutate(
    organism_type = factor(classify_organism_type(target_label), levels = organism_type_levels)
)

fit_detection_trend <- function(data, denominator_var) {
  fit_temporal_model(data, denominator_var)$summary
}

make_monthly_rates <- function(data, label_var, month_seq, monthly_icu_denominators) {
  data %>%
    group_by(calendar_month, organism_label = .data[[label_var]]) %>%
    summarise(n_detection_events = n_distinct(detection_event_id), .groups = "drop") %>%
    complete(
      calendar_month = month_seq,
      organism_label,
      fill = list(n_detection_events = 0L)
    ) %>%
    left_join(monthly_icu_denominators, by = "calendar_month") %>%
    mutate(
      detection_events_per_100_icu_admissions = if_else(
        n_icu_admissions > 0,
        100 * n_detection_events / n_icu_admissions,
        NA_real_
      ),
      detection_events_per_100_icu_days = if_else(
        n_icu_days > 0,
        100 * n_detection_events / n_icu_days,
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
study_start_date <- study_settings$study_start_date
study_end_date <- study_settings$study_end_date
top_n_organisms <- as.integer(Sys.getenv("TOP_N_TREND_ORGANISMS", unset = "25"))
plot_n_increasing <- as.integer(Sys.getenv("PLOT_N_INCREASING_ORGANISMS", unset = "12"))
plot_n_decreasing <- as.integer(Sys.getenv("PLOT_N_DECREASING_ORGANISMS", unset = as.character(plot_n_increasing)))

if (is.na(row_path) || !nzchar(row_path)) {
  row_path <- latest_project_intermediate_file("^icu_culture_rows_.*\\.csv$", "cohort")
}

study_start_dttm <- if (!is.na(study_start_date) && nzchar(study_start_date)) safe_ts(study_start_date) else as.POSIXct(NA)
study_end_dttm <- if (!is.na(study_end_date) && nzchar(study_end_date)) safe_ts(study_end_date) + days(1) else as.POSIXct(NA)

message("Reading ICU culture rows: ", row_path)
if (!is.na(study_start_dttm)) message("Study start: ", study_start_dttm)
if (!is.na(study_end_dttm)) message("Study end: ", study_end_dttm)

rows <- read_clif_csv(row_path, show_col_types = FALSE) %>%
  mutate(
    collect_dttm = safe_ts(collect_dttm),
    calendar_month = floor_date(collect_dttm, "month"),
    fluid_name = coalesce(na_if(fluid_name, ""), "missing"),
    method_name = coalesce(na_if(method_name, ""), "missing"),
    organism_id = as.character(organism_id),
    organism_group = coalesce(na_if(str_to_lower(str_trim(as.character(organism_group))), ""), "missing"),
    organism_category = coalesce(na_if(str_to_lower(str_trim(as.character(organism_category))), ""), organism_group),
    organism_name = coalesce(na_if(str_to_lower(str_trim(as.character(organism_name))), ""), organism_category),
    organism_category_label = clean_label(organism_category),
    organism_text = str_squish(str_c(organism_name, organism_category, organism_group, sep = " ")),
    positive_culture = as.logical(positive_culture),
    explicit_negative_name = drop_negative_organism_name(organism_name),
    detection_event_id = str_c(culture_event_id, organism_category, sep = "|")
  ) %>%
  filter(positive_culture, !explicit_negative_name, !is.na(calendar_month)) %>%
  filter(is.na(study_start_dttm) | collect_dttm >= study_start_dttm) %>%
  filter(is.na(study_end_dttm) | collect_dttm < study_end_dttm)

observation_stays <- read_culture_data(study_start_dttm, study_end_dttm)$icu_admissions
month_seq <- seq(floor_date(min(observation_stays$icu_in_dttm_clipped), "month"), floor_date(max(observation_stays$icu_out_dttm_clipped - seconds(1)), "month"), by = "month")
monthly_icu_denominators <- build_monthly_icu_denominators(month_seq, study_start_dttm, study_end_dttm)

if (nrow(rows) == 0) {
  out_dir <- project_output_dir("organism_trends"); stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write_csv(tibble(calendar_month=as.POSIXct(character(),tz="UTC"),organism_category=character(),n_detection_events=integer(),site_name=character()),file.path(out_dir,glue("monthly_all_organism_detection_counts_{site_name}_{stamp}.csv")))
  write_csv(monthly_icu_denominators,file.path(out_dir,glue("monthly_icu_denominators_for_organism_trends_{site_name}_{stamp}.csv")))
  write_csv(tibble(site_name=site_name,analysis_status="skipped_no_positive_organisms"),file.path(out_dir,glue("organism_trend_availability_{site_name}_{stamp}.csv")))
  message("No positive organisms; trend models skipped.")
  quit(save="no",status=0)
}

top_organism_labels <- rows %>%
  distinct(detection_event_id, organism_category_label) %>%
  count(organism_category_label, name = "total_detection_events", sort = TRUE) %>%
  slice_head(n = top_n_organisms) %>%
  mutate(
    organism_type = factor(classify_organism_type(organism_category_label), levels = organism_type_levels)
  )

top_monthly_rates <- rows %>%
  semi_join(top_organism_labels, by = "organism_category_label") %>%
  make_monthly_rates("organism_category_label", month_seq, monthly_icu_denominators) %>%
  left_join(top_organism_labels, by = c("organism_label" = "organism_category_label"))

target_text_detections <- target_organism_patterns %>%
  tidyr::crossing(row_id = seq_len(nrow(rows))) %>%
  mutate(row_match = str_detect(rows$organism_text[row_id], regex(pattern, ignore_case = TRUE))) %>%
  filter(row_match) %>%
  transmute(
    target_label,
    detection_event_id = rows$detection_event_id[row_id],
    calendar_month = rows$calendar_month[row_id],
    phenotype_source = "organism_text"
  ) %>%
  distinct()

target_detections <- target_text_detections

target_detection_source_summary <- target_detections %>%
  count(target_label, phenotype_source, name = "n_detection_events") %>%
  complete(
    target_label = target_organism_patterns$target_label,
    phenotype_source = "organism_text",
    fill = list(n_detection_events = 0L)
  ) %>%
  arrange(target_label, phenotype_source)

target_monthly_rates <- if (nrow(target_detections) > 0) {
  make_monthly_rates(target_detections, "target_label", month_seq, monthly_icu_denominators) %>%
    mutate(total_detection_events = sum(n_detection_events), .by = organism_label) %>%
    left_join(
      target_organism_patterns %>% select(organism_label = target_label, organism_type),
      by = "organism_label"
    )
} else {
  tidyr::expand_grid(calendar_month = month_seq, organism_label = target_organism_patterns$target_label) %>%
    left_join(monthly_icu_denominators, by = "calendar_month") %>%
    left_join(
      target_organism_patterns %>% select(organism_label = target_label, organism_type),
      by = "organism_label"
    ) %>%
    mutate(
      n_detection_events = 0L,
      detection_events_per_100_icu_admissions = if_else(n_icu_admissions > 0, 0, NA_real_),
      detection_events_per_100_icu_days = if_else(n_icu_days > 0, 0, NA_real_),
      total_detection_events = 0L
    )
}

summarise_trends <- function(data, denominator_var, rate_var, zero_label = "Not estimated") {
  data %>%
    group_by(organism_label) %>%
    group_modify(~ fit_detection_trend(.x, denominator_var)) %>%
    ungroup() %>%
    left_join(
      data %>%
        group_by(organism_label, organism_type) %>%
        summarise(
          total_detection_events = max(total_detection_events, na.rm = TRUE),
          first_nonzero_month = suppressWarnings(min(calendar_month[n_detection_events > 0], na.rm = TRUE)),
          last_nonzero_month = suppressWarnings(max(calendar_month[n_detection_events > 0], na.rm = TRUE)),
          mean_monthly_rate = mean(.data[[rate_var]], na.rm = TRUE),
          .groups = "drop"
        ),
      by = "organism_label"
    ) %>%
    mutate(
      first_nonzero_month = if_else(is.infinite(first_nonzero_month), as.POSIXct(NA), first_nonzero_month),
      last_nonzero_month = if_else(is.infinite(last_nonzero_month), as.POSIXct(NA), last_nonzero_month),
      fdr_p_value = p.adjust(p_value, method = "BH"),
      trend_direction = case_when(
        total_detection_events == 0 ~ zero_label,
        is.na(annual_percent_change) | is.na(p_value) ~ "Not estimated",
        fdr_p_value < 0.05 & !coalesce(residual_dependence_flag, TRUE) & annual_percent_change > 0 ~ "Increasing",
        fdr_p_value < 0.05 & !coalesce(residual_dependence_flag, TRUE) & annual_percent_change < 0 ~ "Decreasing",
        TRUE ~ "No clear trend"
      )
    ) %>%
    arrange(desc(annual_percent_change))
}

trend_summary_top_admissions <- summarise_trends(
  top_monthly_rates,
  "n_icu_admissions",
  "detection_events_per_100_icu_admissions"
)
trend_summary_top_icu_days <- summarise_trends(
  top_monthly_rates,
  "n_icu_days",
  "detection_events_per_100_icu_days"
)
trend_summary_targets_admissions <- summarise_trends(
  target_monthly_rates,
  "n_icu_admissions",
  "detection_events_per_100_icu_admissions",
  zero_label = "Not observed"
)
trend_summary_targets_icu_days <- summarise_trends(
  target_monthly_rates,
  "n_icu_days",
  "detection_events_per_100_icu_days",
  zero_label = "Not observed"
)

increasing_labels <- trend_summary_top_icu_days %>%
  filter(total_detection_events >= 25, !is.na(annual_percent_change), annual_percent_change > 1e-8) %>%
  arrange(desc(annual_percent_change)) %>%
  slice_head(n = plot_n_increasing) %>%
  pull(organism_label)

decreasing_labels <- trend_summary_top_icu_days %>%
  filter(total_detection_events >= 25, !is.na(annual_percent_change), annual_percent_change < -1e-8) %>%
  arrange(annual_percent_change) %>%
  slice_head(n = plot_n_decreasing) %>%
  pull(organism_label)

target_plot_labels <- target_monthly_rates %>%
  group_by(organism_label) %>%
  summarise(total_detection_events = max(total_detection_events, na.rm = TRUE), .groups = "drop") %>%
  filter(total_detection_events > 0 | organism_label %in% c("MRSA (organism text)", "Staphylococcus aureus")) %>%
  arrange(desc(total_detection_events), organism_label) %>%
  pull(organism_label)

prepare_plot_data <- function(data, plot_labels, rate_var) {
  selected <- data %>% filter(organism_label %in% plot_labels)
  denominator <- if (rate_var == "detection_events_per_100_icu_days") "n_icu_days" else "n_icu_admissions"
  fitted <- selected %>% group_by(organism_label) %>% group_modify(~ fit_temporal_model(.x, denominator)$predictions) %>% ungroup()
  selected %>% left_join(fitted, by = c("organism_label", "calendar_month")) %>%
    mutate(plot_rate = .data[[rate_var]], organism_type = factor(organism_type, levels = organism_type_levels), organism_label = factor(organism_label, levels = plot_labels)) %>% arrange(organism_label, calendar_month)
}

plot_increasing_admissions_data <- prepare_plot_data(
  top_monthly_rates,
  increasing_labels,
  "detection_events_per_100_icu_admissions"
)
plot_increasing_icu_days_data <- prepare_plot_data(
  top_monthly_rates,
  increasing_labels,
  "detection_events_per_100_icu_days"
)
plot_decreasing_admissions_data <- prepare_plot_data(
  top_monthly_rates,
  decreasing_labels,
  "detection_events_per_100_icu_admissions"
)
plot_decreasing_icu_days_data <- prepare_plot_data(
  top_monthly_rates,
  decreasing_labels,
  "detection_events_per_100_icu_days"
)
plot_target_admissions_data <- prepare_plot_data(
  target_monthly_rates,
  target_plot_labels,
  "detection_events_per_100_icu_admissions"
)
plot_target_icu_days_data <- prepare_plot_data(
  target_monthly_rates,
  target_plot_labels,
  "detection_events_per_100_icu_days"
)

taxonomy_palette <- organism_type_base_colors[organism_type_levels]

plot_trend_facets <- function(data, title, y_label, ncol = 3) {
  if (!nrow(data)) return(ggplot() + annotate("text", x = 0, y = 0, label = "No estimable net changes in this direction") + labs(title = title) + theme_void())
  available_taxonomy_palette <- taxonomy_palette[names(taxonomy_palette) %in% unique(as.character(data$organism_type))]

  ggplot(data, aes(calendar_month, plot_rate)) +
    geom_col(aes(fill = organism_type), width = 25 * 24 * 60 * 60) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey60", alpha = 0.2, na.rm = TRUE) +
    geom_line(aes(y = fitted_value), color = "black", linewidth = 0.65, na.rm = TRUE) +
    facet_wrap(vars(organism_label), scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = available_taxonomy_palette, drop = FALSE) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = comma, limits = c(0, NA)) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      fill = "Taxonomy",
      caption = NULL) +
    theme_trends
}

out_dir <- project_output_dir("organism_trends")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

paths <- c(
  trend_summary_top_admissions = file.path(out_dir, glue("organism_detection_trend_screen_top_per_100_icu_admissions_{site_name}_{stamp}.csv")),
  trend_summary_top_icu_days = file.path(out_dir, glue("organism_detection_trend_screen_top_per_100_icu_days_{site_name}_{stamp}.csv")),
  trend_summary_targets_admissions = file.path(out_dir, glue("organism_detection_trend_screen_targets_per_100_icu_admissions_{site_name}_{stamp}.csv")),
  trend_summary_targets_icu_days = file.path(out_dir, glue("organism_detection_trend_screen_targets_per_100_icu_days_{site_name}_{stamp}.csv")),
  target_detection_source_summary = file.path(out_dir, glue("target_resistance_detection_source_summary_{site_name}_{stamp}.csv")),
  monthly_top_rates = file.path(out_dir, glue("monthly_top_organism_detection_rates_{site_name}_{stamp}.csv")),
  monthly_target_rates = file.path(out_dir, glue("monthly_target_organism_detection_rates_{site_name}_{stamp}.csv")),
  monthly_icu_denominators = file.path(out_dir, glue("monthly_icu_denominators_for_organism_trends_{site_name}_{stamp}.csv")),
  plotted_increasing_organisms = file.path(out_dir, glue("fastest_increasing_organisms_plotted_{site_name}_{stamp}.csv")),
  plotted_decreasing_organisms = file.path(out_dir, glue("fastest_decreasing_organisms_plotted_{site_name}_{stamp}.csv")),
  increasing_plot_admissions = file.path(out_dir, glue("monthly_fastest_increasing_organism_detection_rates_per_100_icu_admissions_{site_name}_{stamp}.png")),
  increasing_plot_icu_days = file.path(out_dir, glue("monthly_fastest_increasing_organism_detection_rates_per_100_icu_days_{site_name}_{stamp}.png")),
  decreasing_plot_admissions = file.path(out_dir, glue("monthly_fastest_decreasing_organism_detection_rates_per_100_icu_admissions_{site_name}_{stamp}.png")),
  decreasing_plot_icu_days = file.path(out_dir, glue("monthly_fastest_decreasing_organism_detection_rates_per_100_icu_days_{site_name}_{stamp}.png")),
  target_plot_admissions = file.path(out_dir, glue("monthly_target_organism_detection_rates_per_100_icu_admissions_{site_name}_{stamp}.png")),
  target_plot_icu_days = file.path(out_dir, glue("monthly_target_organism_detection_rates_per_100_icu_days_{site_name}_{stamp}.png"))
)

fitted_trends <- top_monthly_rates %>% group_by(organism_label) %>% group_modify(~ fit_temporal_model(.x, "n_icu_days")$predictions) %>% ungroup()
write_csv(fitted_trends, file.path(out_dir, glue("monthly_organism_gam_fitted_rates_per_100_icu_days_{site_name}_{stamp}.csv")))

write_csv(trend_summary_top_admissions, paths[["trend_summary_top_admissions"]])
write_csv(trend_summary_top_icu_days, paths[["trend_summary_top_icu_days"]])
write_csv(trend_summary_targets_admissions, paths[["trend_summary_targets_admissions"]])
write_csv(trend_summary_targets_icu_days, paths[["trend_summary_targets_icu_days"]])
write_csv(target_detection_source_summary, paths[["target_detection_source_summary"]])
# Complete aggregate counts for central pooling, independent of site-specific top-N ranks.
all_monthly_counts <- rows %>% group_by(calendar_month, organism_category) %>%
  summarise(n_detection_events = n_distinct(detection_event_id), .groups = "drop") %>%
  complete(calendar_month = month_seq, organism_category, fill = list(n_detection_events = 0L)) %>%
  mutate(site_name = site_name)
write_csv(all_monthly_counts, file.path(out_dir, glue("monthly_all_organism_detection_counts_{site_name}_{stamp}.csv")))
write_csv(top_monthly_rates, paths[["monthly_top_rates"]])
write_csv(target_monthly_rates, paths[["monthly_target_rates"]])
write_csv(monthly_icu_denominators, paths[["monthly_icu_denominators"]])
write_csv(
  trend_summary_top_icu_days %>% filter(organism_label %in% increasing_labels),
  paths[["plotted_increasing_organisms"]]
)
write_csv(
  trend_summary_top_icu_days %>% filter(organism_label %in% decreasing_labels),
  paths[["plotted_decreasing_organisms"]]
)

p_increasing_admissions <- plot_trend_facets(
  plot_increasing_admissions_data,
  "Largest Estimated Increases in Organism Detection Rates per 100 ICU Admissions",
  "Detection events per 100 ICU admissions",
  ncol = 3
)
p_increasing_icu_days <- plot_trend_facets(
  plot_increasing_icu_days_data,
  "Largest Estimated Increases in Organism Detection Rates per 100 ICU Days",
  "Detection events per 100 ICU days",
  ncol = 3
)
p_decreasing_admissions <- plot_trend_facets(
  plot_decreasing_admissions_data,
  "Largest Estimated Decreases in Organism Detection Rates per 100 ICU Admissions",
  "Detection events per 100 ICU admissions",
  ncol = 3
)
p_decreasing_icu_days <- plot_trend_facets(
  plot_decreasing_icu_days_data,
  "Largest Estimated Decreases in Organism Detection Rates per 100 ICU Days",
  "Detection events per 100 ICU days",
  ncol = 3
)
p_targets_admissions <- plot_trend_facets(
  plot_target_admissions_data,
  "Target Organism Detection Rates per 100 ICU Admissions",
  "Detection events per 100 ICU admissions",
  ncol = 2
)
p_targets_icu_days <- plot_trend_facets(
  plot_target_icu_days_data,
  "Target Organism Detection Rates per 100 ICU Days",
  "Detection events per 100 ICU days",
  ncol = 2
)

ggsave(paths[["increasing_plot_admissions"]], p_increasing_admissions, width = 14, height = 12, dpi = 300)
ggsave(paths[["increasing_plot_icu_days"]], p_increasing_icu_days, width = 14, height = 12, dpi = 300)
ggsave(paths[["decreasing_plot_admissions"]], p_decreasing_admissions, width = 14, height = 12, dpi = 300)
ggsave(paths[["decreasing_plot_icu_days"]], p_decreasing_icu_days, width = 14, height = 12, dpi = 300)
ggsave(paths[["target_plot_admissions"]], p_targets_admissions, width = 12, height = 10, dpi = 300)
ggsave(paths[["target_plot_icu_days"]], p_targets_icu_days, width = 12, height = 10, dpi = 300)

message("Wrote trend summaries:")
print(paths[c(
  "trend_summary_top_admissions",
  "trend_summary_top_icu_days",
  "trend_summary_targets_admissions",
  "trend_summary_targets_icu_days",
  "target_detection_source_summary",
  "monthly_top_rates",
  "monthly_target_rates",
  "plotted_increasing_organisms",
  "plotted_decreasing_organisms"
)])
message("Wrote plots:")
print(paths[c(
  "increasing_plot_admissions",
  "increasing_plot_icu_days",
  "decreasing_plot_admissions",
  "decreasing_plot_icu_days",
  "target_plot_admissions",
  "target_plot_icu_days"
)])

message("Top increasing organisms by annual percent change:")
print(
  trend_summary_top_admissions %>%
    select(organism_label, total_detection_events, annual_percent_change, p_value, trend_direction) %>%
    head(15),
  n = 15,
  width = Inf
)

message("Top decreasing organisms by annual percent change:")
print(
  trend_summary_top_admissions %>%
    arrange(annual_percent_change) %>%
    select(organism_label, total_detection_events, annual_percent_change, p_value, trend_direction) %>%
    head(15),
  n = 15,
  width = Inf
)

message("Target organism trend summary:")
print(
  trend_summary_targets_admissions %>%
    select(organism_label, total_detection_events, annual_percent_change, p_value, trend_direction) %>%
    arrange(desc(total_detection_events)),
  n = Inf,
  width = Inf
)
