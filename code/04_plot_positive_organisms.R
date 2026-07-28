# ================================================================================================
# Plot Organisms Detected in Positive ICU Cultures
#
# Input:
#   ICU culture row output from code/01_identify_icu_culture_cohort.R.
#
# Outputs:
#   Positive organism detection summaries and PNG figures overall and by culture category.
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

latest_file <- function(pattern, path = file.path("output", "cohort")) {
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

drop_negative_organism_name <- function(x) {
  str_detect(
    coalesce(x, ""),
    regex("^(no|none|not) .*isolated|no growth|no .* detected|negative for", ignore_case = TRUE)
  )
}

make_within_label <- function(label, facet) paste(label, facet, sep = "___")
strip_within_label <- function(x) str_replace(x, "___.*$", "")

site_name <- Sys.getenv("CLIF_SITE_NAME", unset = "SITE")
row_path <- Sys.getenv("ICU_CULTURE_ROWS_PATH", unset = NA_character_)
top_n_culture_types <- as.integer(Sys.getenv("TOP_N_CULTURE_TYPES", unset = "8"))
top_n_overall <- as.integer(Sys.getenv("TOP_N_ORGANISMS_OVERALL", unset = "20"))
top_n_per_type <- as.integer(Sys.getenv("TOP_N_ORGANISMS_PER_TYPE", unset = "10"))
plot_start_date <- Sys.getenv("PLOT_START_DATE", unset = NA_character_)
plot_end_date <- Sys.getenv("PLOT_END_DATE", unset = NA_character_)

if (is.na(row_path) || !nzchar(row_path)) {
  row_path <- latest_file("^icu_culture_rows_.*\\.csv$")
}

message("Reading ICU culture rows: ", row_path)

rows <- readr::read_csv(row_path, show_col_types = FALSE) %>%
  mutate(
    collect_dttm = safe_ts(collect_dttm),
    fluid_category = coalesce(na_if(fluid_category, ""), "missing"),
    culture_type = clean_label(fluid_category),
    organism_group = coalesce(na_if(organism_group, ""), "missing"),
    organism_category = coalesce(na_if(organism_category, ""), organism_group, "missing"),
    organism_name = coalesce(na_if(organism_name, ""), organism_category),
    organism_group_label = clean_label(organism_group),
    organism_category_label = clean_label(organism_category),
    organism_name_label = str_to_sentence(str_squish(organism_name)),
    positive_culture = as.logical(positive_culture),
    explicit_negative_name = drop_negative_organism_name(organism_name)
  )

plot_start_dttm <- if (!is.na(plot_start_date) && nzchar(plot_start_date)) safe_ts(plot_start_date) else as.POSIXct(NA)
plot_end_dttm <- if (!is.na(plot_end_date) && nzchar(plot_end_date)) safe_ts(plot_end_date) + days(1) - seconds(1) else as.POSIXct(NA)

positive_rows <- rows %>%
  filter(positive_culture, !explicit_negative_name) %>%
  filter(is.na(plot_start_dttm) | collect_dttm >= plot_start_dttm) %>%
  filter(is.na(plot_end_dttm) | collect_dttm <= plot_end_dttm)

if (nrow(positive_rows) == 0) {
  stop("No positive organism rows after filters.")
}

if (!is.na(plot_start_dttm)) message("Plot start: ", plot_start_dttm)
if (!is.na(plot_end_dttm)) message("Plot end: ", plot_end_dttm)

top_culture_types <- positive_rows %>%
  count(culture_type, sort = TRUE) %>%
  slice_head(n = top_n_culture_types) %>%
  pull(culture_type)

positive_rows <- positive_rows %>%
  mutate(culture_type_plot = if_else(culture_type %in% top_culture_types, culture_type, "Other"))

summarise_detection <- function(data, organism_var, label_var) {
  data %>%
    group_by(.data[[organism_var]], .data[[label_var]]) %>%
    summarise(
      n_detection_rows = n(),
      n_culture_events = n_distinct(paste(patient_id, hospitalization_id, icu_interval_id, collect_dttm, fluid_name, method_name)),
      n_hospitalizations = n_distinct(hospitalization_id),
      n_patients = n_distinct(patient_id),
      .groups = "drop"
    ) %>%
    rename(organism = all_of(organism_var), organism_label = all_of(label_var)) %>%
    arrange(desc(n_detection_rows), organism_label)
}

summarise_detection_by_type <- function(data, organism_var, label_var) {
  data %>%
    group_by(culture_type = culture_type_plot, .data[[organism_var]], .data[[label_var]]) %>%
    summarise(
      n_detection_rows = n(),
      n_culture_events = n_distinct(paste(patient_id, hospitalization_id, icu_interval_id, collect_dttm, fluid_name, method_name)),
      n_hospitalizations = n_distinct(hospitalization_id),
      n_patients = n_distinct(patient_id),
      .groups = "drop"
    ) %>%
    rename(organism = all_of(organism_var), organism_label = all_of(label_var)) %>%
    arrange(culture_type, desc(n_detection_rows), organism_label)
}

organism_group_overall <- summarise_detection(positive_rows, "organism_group", "organism_group_label")
organism_category_overall <- summarise_detection(positive_rows, "organism_category", "organism_category_label")
organism_name_overall <- summarise_detection(positive_rows, "organism_name", "organism_name_label")

organism_group_by_type <- summarise_detection_by_type(positive_rows, "organism_group", "organism_group_label")
organism_category_by_type <- summarise_detection_by_type(positive_rows, "organism_category", "organism_category_label")
organism_name_by_type <- summarise_detection_by_type(positive_rows, "organism_name", "organism_name_label")

top_by_type <- function(data, n = top_n_per_type) {
  data %>%
    group_by(culture_type) %>%
    slice_max(n_detection_rows, n = n, with_ties = FALSE) %>%
    ungroup()
}

out_dir <- file.path("output", "organisms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

paths <- c(
  organism_group_overall = file.path(out_dir, glue("positive_organism_group_overall_{site_name}_{stamp}.csv")),
  organism_group_by_type = file.path(out_dir, glue("positive_organism_group_by_culture_type_{site_name}_{stamp}.csv")),
  organism_category_overall = file.path(out_dir, glue("positive_organism_category_overall_{site_name}_{stamp}.csv")),
  organism_category_by_type = file.path(out_dir, glue("positive_organism_category_by_culture_type_{site_name}_{stamp}.csv")),
  organism_name_overall = file.path(out_dir, glue("positive_organism_name_overall_{site_name}_{stamp}.csv")),
  organism_name_by_type = file.path(out_dir, glue("positive_organism_name_by_culture_type_{site_name}_{stamp}.csv"))
)

write_csv(organism_group_overall, paths[["organism_group_overall"]])
write_csv(organism_group_by_type, paths[["organism_group_by_type"]])
write_csv(organism_category_overall, paths[["organism_category_overall"]])
write_csv(organism_category_by_type, paths[["organism_category_by_type"]])
write_csv(organism_name_overall, paths[["organism_name_overall"]])
write_csv(organism_name_by_type, paths[["organism_name_by_type"]])

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    legend.position = "none"
  )

plot_overall_bar <- function(data, title, fill) {
  data %>%
    slice_head(n = top_n_overall) %>%
    mutate(organism_label = fct_reorder(organism_label, n_detection_rows)) %>%
    ggplot(aes(n_detection_rows, organism_label)) +
    geom_col(fill = fill) +
    scale_x_continuous(labels = comma) +
    labs(title = title, x = "Positive organism detection rows", y = NULL) +
    plot_theme
}

plot_faceted_bar <- function(data, title, fill) {
  top_by_type(data) %>%
    mutate(
      organism_label_within = make_within_label(organism_label, culture_type),
      organism_label_within = fct_reorder(organism_label_within, n_detection_rows)
    ) %>%
    ggplot(aes(n_detection_rows, organism_label_within)) +
    geom_col(fill = fill) +
    facet_wrap(vars(culture_type), scales = "free", ncol = 2) +
    scale_x_continuous(labels = comma) +
    scale_y_discrete(labels = strip_within_label) +
    labs(title = title, x = "Positive organism detection rows", y = NULL) +
    plot_theme
}

p_group_overall <- plot_overall_bar(
  organism_group_overall,
  "Predominant Organism Groups in Positive ICU Cultures",
  "#4C78A8"
)

p_category_overall <- plot_overall_bar(
  organism_category_overall,
  "Predominant Organisms in Positive ICU Cultures",
  "#54A24B"
)

p_group_by_type <- plot_faceted_bar(
  organism_group_by_type,
  "Predominant Organism Groups by Culture Type",
  "#4C78A8"
)

p_category_by_type <- plot_faceted_bar(
  organism_category_by_type,
  "Predominant Organisms by Culture Type",
  "#54A24B"
)

plot_paths <- c(
  organism_group_overall = file.path(out_dir, glue("positive_organism_group_overall_{site_name}_{stamp}.png")),
  organism_category_overall = file.path(out_dir, glue("positive_organism_category_overall_{site_name}_{stamp}.png")),
  organism_group_by_type = file.path(out_dir, glue("positive_organism_group_by_culture_type_{site_name}_{stamp}.png")),
  organism_category_by_type = file.path(out_dir, glue("positive_organism_category_by_culture_type_{site_name}_{stamp}.png"))
)

ggsave(plot_paths[["organism_group_overall"]], p_group_overall, width = 9, height = 7, dpi = 300)
ggsave(plot_paths[["organism_category_overall"]], p_category_overall, width = 9, height = 7, dpi = 300)
ggsave(plot_paths[["organism_group_by_type"]], p_group_by_type, width = 12, height = 12, dpi = 300)
ggsave(plot_paths[["organism_category_by_type"]], p_category_by_type, width = 12, height = 12, dpi = 300)

message("Positive organism rows after excluding explicit negative organism names: ", nrow(positive_rows))
message("Excluded explicit negative organism-name rows from detected-organism summaries: ", sum(rows$positive_culture & rows$explicit_negative_name, na.rm = TRUE))
message("")
message("Top organism groups:")
print(head(organism_group_overall, 20), n = 20, width = Inf)
message("")
message("Top organisms:")
print(head(organism_category_overall, 20), n = 20, width = Inf)
message("")
message("Culture types displayed separately:")
print(tibble(culture_type = top_culture_types), n = Inf)
message("")
message("Wrote summaries:")
print(paths)
message("")
message("Wrote plots:")
print(plot_paths)
