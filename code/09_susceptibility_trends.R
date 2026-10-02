# Optional analysis for sites with BOTH microbiology tables.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(lubridate); library(ggplot2); library(glue)})
source("utils/clif_io.R")
source("utils/susceptibility.R")
source("utils/trends.R")
out_dir <- project_output_dir("susceptibility")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
path <- find_table_path("microbiology_susceptibility", required = FALSE)
availability <- tibble(site_name = clif_site_name, susceptibility_table_available = !is.na(path), analysis_status = if (is.na(path)) "skipped_table_unavailable" else "available")
write_csv(availability, file.path(out_dir, glue("susceptibility_availability_{clif_site_name}_{stamp}.csv")))
if (!is.na(path)) {
  start <- safe_ts(config_value(config, "study_start_date", env = "STUDY_START_DATE", default = NA_character_))
  end <- safe_ts(config_value(config, "study_end_date", env = "STUDY_END_DATE", default = NA_character_)) + days(1)
  rows <- read_clif_csv(latest_project_intermediate_file("^icu_culture_rows_.*\\.csv$", "cohort"), show_col_types = FALSE) %>% mutate(collect_dttm = safe_ts(collect_dttm), organism_id = as.character(organism_id))
  stays <- read_culture_data(start, end)$icu_admissions
  month_seq <- seq(floor_date(min(stays$icu_in_dttm_clipped), "month"), floor_date(max(stays$icu_out_dttm_clipped - seconds(1)), "month"), by = "month")
  result <- build_susceptibility_analysis(rows, read_tbl("microbiology_susceptibility"), monthly_icu_denominators(stays, month_seq), project_path("config", "mcide"))
  write_csv(mutate(result$qc, site_name = clif_site_name), file.path(out_dir, glue("susceptibility_qc_{clif_site_name}_{stamp}.csv")))
  if (nrow(result$monthly)) {
    monthly <- result$monthly %>% mutate(site_name = clif_site_name)
    write_csv(monthly, file.path(out_dir, glue("monthly_organism_antimicrobial_susceptibility_{clif_site_name}_{stamp}.csv")))
    coverage_min <- as.numeric(Sys.getenv("AST_MIN_TESTING_FRACTION", "0.5"))
    if (!is.finite(coverage_min) || coverage_min < 0 || coverage_min > 1) stop("AST_MIN_TESTING_FRACTION must be between 0 and 1.")
    keys <- c("organism_category", "antimicrobial_category", "specimen_stratum")
    fits <- monthly %>% group_by(across(all_of(keys))) %>% group_modify(function(.x, .y) {
      screened <- .x %>% filter(!is.na(testing_fraction), testing_fraction >= coverage_min)
      bind_rows(
        fit_temporal_model(screened, "n_icu_days", "n_susceptible")$summary %>% mutate(outcome = "susceptible_detection_rate"),
        fit_temporal_model(screened, "n_icu_days", "n_non_susceptible")$summary %>% mutate(outcome = "non_susceptible_detection_rate"),
        fit_temporal_model(.x, "n_interpretable", "n_non_susceptible", proportion = TRUE)$summary %>% mutate(outcome = "non_susceptible_fraction")
      ) %>% mutate(minimum_testing_fraction_for_rate = coverage_min)
    }) %>% ungroup() %>% group_by(outcome, specimen_stratum) %>% mutate(fdr_p_value = p.adjust(p_value, "BH"), direction = classify_temporal_direction(model_status, fdr_p_value, annual_percent_change, residual_dependence_flag)) %>% ungroup() %>% mutate(site_name = clif_site_name, effect_scale = if_else(outcome == "non_susceptible_fraction", "annualized endpoint odds ratio", "annualized endpoint rate ratio"))
    write_csv(fits, file.path(out_dir, glue("susceptibility_temporal_models_{clif_site_name}_{stamp}.csv")))
    # Curves, including intervals, for every estimable pair; keep specimens separate.
    curves <- monthly %>% group_by(across(all_of(keys))) %>% group_modify(function(.x, .y) {
      screened <- filter(.x, !is.na(testing_fraction), testing_fraction >= coverage_min)
      bind_rows(fit_temporal_model(screened, "n_icu_days", "n_susceptible")$predictions %>% mutate(outcome = "susceptible_detection_rate"), fit_temporal_model(screened, "n_icu_days", "n_non_susceptible")$predictions %>% mutate(outcome = "non_susceptible_detection_rate"), fit_temporal_model(.x, "n_interpretable", "n_non_susceptible", proportion = TRUE)$predictions %>% mutate(outcome = "non_susceptible_fraction"))
    }) %>% ungroup()
    write_csv(curves, file.path(out_dir, glue("susceptibility_gam_fitted_curves_{clif_site_name}_{stamp}.csv")))
    top_pairs <- monthly %>% filter(specimen_stratum == "Overall") %>% group_by(organism_category, antimicrobial_category) %>% summarise(n_tested = sum(n_interpretable), .groups = "drop") %>% slice_max(n_tested, n = 12, with_ties = FALSE)
    plot_data <- monthly %>% filter(specimen_stratum == "Overall") %>% semi_join(top_pairs, by = c("organism_category", "antimicrobial_category")) %>% mutate(pair = paste(organism_category, antimicrobial_category, sep = " / "))
    if (nrow(plot_data)) {
      p <- ggplot(plot_data, aes(calendar_month, non_susceptible_fraction)) + geom_point(aes(size = n_interpretable), alpha = 0.5, na.rm = TRUE) + facet_wrap(~pair, ncol = 3) + scale_y_continuous(limits = c(0, 1), labels = scales::percent) + labs(x = NULL, y = "Non-susceptible fraction among interpretable tests", size = "Tested isolates", caption = "Points: observed fractions. Curve and 95% CI: season-adjusted long-term mean.
Unavailable and indeterminate results excluded; interpret alongside monthly testing coverage.") + theme_bw()
      if (nrow(curves)) p <- p + geom_ribbon(data = curves %>% filter(specimen_stratum == "Overall", outcome == "non_susceptible_fraction") %>% semi_join(top_pairs, by = c("organism_category", "antimicrobial_category")) %>% mutate(pair = paste(organism_category, antimicrobial_category, sep = " / ")), aes(y = fitted_value, ymin = ci_low, ymax = ci_high), alpha = 0.15) + geom_line(data = curves %>% filter(specimen_stratum == "Overall", outcome == "non_susceptible_fraction") %>% semi_join(top_pairs, by = c("organism_category", "antimicrobial_category")) %>% mutate(pair = paste(organism_category, antimicrobial_category, sep = " / ")), aes(y = fitted_value), color = "#9C3333")
      ggsave(file.path(out_dir, glue("non_susceptible_tested_fraction_{clif_site_name}_{stamp}.png")), p, width = 14, height = 10, dpi = 200)
      rates <- plot_data %>% pivot_longer(c(susceptible_per_100_icu_days, non_susceptible_per_100_icu_days), names_to = "status", values_to = "rate")
      ggsave(file.path(out_dir, glue("susceptibility_detection_rates_{clif_site_name}_{stamp}.png")), ggplot(rates, aes(calendar_month, rate, color = status)) + geom_line(na.rm = TRUE) + facet_wrap(~pair, scales = "free_y", ncol = 3) + labs(x = NULL, y = "Tested isolate detections per 100 ICU days", color = NULL) + theme_bw() + theme(legend.position = "bottom"), width = 14, height = 10, dpi = 200)
    }
  } else message("No linked, mapped ICU susceptibility tests; see QC output.")
} else message("Susceptibility analysis skipped: table unavailable. Text resistance labels are not substituted for AST.")
