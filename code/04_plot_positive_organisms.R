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

classify_organism_type <- function(x) {
  x_clean <- str_to_lower(coalesce(x, ""))
  case_when(
    str_detect(x_clean, "candida|yeast|fung|aspergillus|cryptococcus|mold|mould|saccharomyces") ~ "Fungi/yeast",
    str_detect(x_clean, "mycobacter|\\bafb\\b|acid fast") ~ "Mycobacteria/AFB",
    str_detect(x_clean, "anaerob|bacteroides|clostrid|prevotella|fusobacter|cutibacter|propionibacter") ~ "Anaerobes",
    str_detect(x_clean, "staphylococcus|streptococcus|enterococcus|bacillus|corynebacter|lactobacillus|listeria|gram_positive|gram positive|gpc|coag_pos|coag_neg|coagneg") ~ "Gram positive bacteria",
    str_detect(x_clean, "pseudomonas|stenotrophomonas|xanthomonas|klebsiella|enterobacter|escherichia|serratia|haemophilus|citrobacter|acinetobacter|proteus|moraxella|neisseria|salmonella|shigella|campylobacter|burkholderia|cepacia|legionella|gram_negative|gram negative|gnr") ~ "Gram negative bacteria",
    str_detect(x_clean, "virus|viral|covid|sars|influenza|rsv|herpes|cmv|adenovirus") ~ "Virus",
    str_detect(x_clean, "bacteria|gram|cocci|bacilli|rods") ~ "Other bacteria",
    TRUE ~ "Other/unspecified"
  )
}

organism_type_levels <- c(
  "Gram positive bacteria",
  "Gram negative bacteria",
  "Fungi/yeast",
  "Mycobacteria/AFB",
  "Anaerobes",
  "Other bacteria",
  "Virus",
  "Other/unspecified"
)

organism_type_base_colors <- c(
  "Gram positive bacteria" = "#2F6C99",
  "Gram negative bacteria" = "#D55E00",
  "Fungi/yeast" = "#7B4AB8",
  "Mycobacteria/AFB" = "#008B8B",
  "Anaerobes" = "#8C6D31",
  "Other bacteria" = "#4E8F4A",
  "Virus" = "#C44E52",
  "Other/unspecified" = "#6B7280"
)

organism_palette_caption <- paste(
  "Color:",
  "blue Gram+;",
  "orange/brown Gram-;",
  "purple fungi/yeast;",
  "green other bacteria;",
  "gray other/unspecified."
)

make_shade_ramp <- function(base_color, n) {
  if (n <= 1) return(base_color)
  grDevices::colorRampPalette(c("#F2F4F8", base_color, "#263238"))(n + 2)[2:(n + 1)]
}

make_organism_palette <- function(data) {
  palette_data <- data %>%
    mutate(
      organism_type = factor(coalesce(organism_type, "Other/unspecified"), levels = organism_type_levels),
      n_detection_rows = coalesce(n_detection_rows, 0L)
    ) %>%
    group_by(organism_label, organism_type) %>%
    summarise(total_detection_rows = sum(n_detection_rows), .groups = "drop") %>%
    arrange(organism_type, desc(total_detection_rows), organism_label) %>%
    group_by(organism_type) %>%
    mutate(
      type_index = row_number(),
      type_n = n(),
      color = make_shade_ramp(organism_type_base_colors[as.character(first(organism_type))], first(type_n))[type_index]
    ) %>%
    ungroup()

  setNames(palette_data$color, palette_data$organism_label)
}

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
    culture_month = floor_date(collect_dttm, "month"),
    fluid_category = coalesce(na_if(fluid_category, ""), "missing"),
    culture_type = clean_label(fluid_category),
    culture_type = if_else(culture_type %in% c("Other", "Other unspecified"), "Other", culture_type),
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
    mutate(organism_type = classify_organism_type(organism)) %>%
    relocate(organism_type, .after = organism_label) %>%
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
    mutate(organism_type = classify_organism_type(organism)) %>%
    relocate(organism_type, .after = organism_label) %>%
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

month_seq <- seq(min(positive_rows$culture_month), max(positive_rows$culture_month), by = "month")

monthly_detection <- function(data, organism_var, label_var, top_data) {
  organism_sym <- rlang::sym(organism_var)
  data %>%
    inner_join(
      top_data %>% select(organism, organism_label) %>% distinct(),
      by = setNames("organism", organism_var)
    ) %>%
    group_by(culture_month, !!organism_sym, organism_label) %>%
    summarise(
      n_detection_rows = n(),
      n_culture_events = n_distinct(paste(patient_id, hospitalization_id, icu_interval_id, collect_dttm, fluid_name, method_name)),
      .groups = "drop"
    ) %>%
    complete(
      culture_month = month_seq,
      nesting(!!organism_sym, organism_label),
      fill = list(n_detection_rows = 0L, n_culture_events = 0L)
    ) %>%
    rename(organism = all_of(organism_var)) %>%
    mutate(organism_type = classify_organism_type(organism)) %>%
    relocate(organism_type, .after = organism_label)
}

monthly_detection_by_type <- function(data, organism_var, label_var, top_data) {
  organism_sym <- rlang::sym(organism_var)
  data %>%
    inner_join(
      top_data %>%
        select(culture_type, organism, organism_label) %>%
        distinct(),
      by = c("culture_type_plot" = "culture_type", setNames("organism", organism_var))
    ) %>%
    group_by(culture_month, culture_type = culture_type_plot, !!organism_sym, organism_label) %>%
    summarise(
      n_detection_rows = n(),
      n_culture_events = n_distinct(paste(patient_id, hospitalization_id, icu_interval_id, collect_dttm, fluid_name, method_name)),
      .groups = "drop"
    ) %>%
    complete(
      culture_month = month_seq,
      nesting(culture_type, !!organism_sym, organism_label),
      fill = list(n_detection_rows = 0L, n_culture_events = 0L)
    ) %>%
    rename(organism = all_of(organism_var)) %>%
    mutate(organism_type = classify_organism_type(organism)) %>%
    relocate(organism_type, .after = organism_label)
}

top_group_overall <- organism_group_overall %>% slice_head(n = top_n_overall)
top_category_overall <- organism_category_overall %>% slice_head(n = top_n_overall)
top_group_by_type <- top_by_type(organism_group_by_type)
top_category_by_type <- top_by_type(organism_category_by_type)

monthly_organism_group_overall <- monthly_detection(
  positive_rows,
  "organism_group",
  "organism_group_label",
  top_group_overall
)
monthly_organism_category_overall <- monthly_detection(
  positive_rows,
  "organism_category",
  "organism_category_label",
  top_category_overall
)
monthly_organism_group_by_type <- monthly_detection_by_type(
  positive_rows,
  "organism_group",
  "organism_group_label",
  top_group_by_type
)
monthly_organism_category_by_type <- monthly_detection_by_type(
  positive_rows,
  "organism_category",
  "organism_category_label",
  top_category_by_type
)

out_dir <- file.path("output", "organisms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

paths <- c(
  organism_group_overall = file.path(out_dir, glue("positive_organism_group_overall_{site_name}_{stamp}.csv")),
  organism_group_by_type = file.path(out_dir, glue("positive_organism_group_by_culture_type_{site_name}_{stamp}.csv")),
  organism_category_overall = file.path(out_dir, glue("positive_organism_category_overall_{site_name}_{stamp}.csv")),
  organism_category_by_type = file.path(out_dir, glue("positive_organism_category_by_culture_type_{site_name}_{stamp}.csv")),
  organism_name_overall = file.path(out_dir, glue("positive_organism_name_overall_{site_name}_{stamp}.csv")),
  organism_name_by_type = file.path(out_dir, glue("positive_organism_name_by_culture_type_{site_name}_{stamp}.csv")),
  monthly_organism_group_overall = file.path(out_dir, glue("monthly_positive_organism_group_overall_{site_name}_{stamp}.csv")),
  monthly_organism_group_by_type = file.path(out_dir, glue("monthly_positive_organism_group_by_culture_type_{site_name}_{stamp}.csv")),
  monthly_organism_category_overall = file.path(out_dir, glue("monthly_positive_organism_category_overall_{site_name}_{stamp}.csv")),
  monthly_organism_category_by_type = file.path(out_dir, glue("monthly_positive_organism_category_by_culture_type_{site_name}_{stamp}.csv"))
)

write_csv(organism_group_overall, paths[["organism_group_overall"]])
write_csv(organism_group_by_type, paths[["organism_group_by_type"]])
write_csv(organism_category_overall, paths[["organism_category_overall"]])
write_csv(organism_category_by_type, paths[["organism_category_by_type"]])
write_csv(organism_name_overall, paths[["organism_name_overall"]])
write_csv(organism_name_by_type, paths[["organism_name_by_type"]])
write_csv(monthly_organism_group_overall, paths[["monthly_organism_group_overall"]])
write_csv(monthly_organism_group_by_type, paths[["monthly_organism_group_by_type"]])
write_csv(monthly_organism_category_overall, paths[["monthly_organism_category_overall"]])
write_csv(monthly_organism_category_by_type, paths[["monthly_organism_category_by_type"]])

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    legend.position = "none"
  )

stacked_time_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    legend.position = "bottom"
  )

plot_overall_bar <- function(data, title) {
  plot_data <- data %>%
    slice_head(n = top_n_overall) %>%
    mutate(organism_label = fct_reorder(organism_label, n_detection_rows)) %>%
    arrange(factor(organism_type, levels = organism_type_levels), desc(n_detection_rows), organism_label)

  ggplot(plot_data, aes(n_detection_rows, organism_label, fill = organism_label)) +
    geom_col(color = "white", linewidth = 0.15) +
    scale_fill_manual(values = make_organism_palette(plot_data)) +
    scale_x_continuous(labels = comma) +
    labs(title = title, x = "Positive organism detection rows", y = NULL, caption = organism_palette_caption) +
    plot_theme
}

plot_faceted_bar <- function(data, title) {
  plot_data <- top_by_type(data) %>%
    mutate(
      organism_label_within = make_within_label(organism_label, culture_type),
      organism_label_within = fct_reorder(organism_label_within, n_detection_rows)
    )

  ggplot(plot_data, aes(n_detection_rows, organism_label_within, fill = organism_label)) +
    geom_col(color = "white", linewidth = 0.15) +
    facet_wrap(vars(culture_type), scales = "free", ncol = 2) +
    scale_fill_manual(values = make_organism_palette(plot_data)) +
    scale_x_continuous(labels = comma) +
    scale_y_discrete(labels = strip_within_label) +
    labs(title = title, x = "Positive organism detection rows", y = NULL, caption = organism_palette_caption) +
    plot_theme
}

plot_overall_stacked_bar <- function(data, title) {
  plot_data <- data %>%
    mutate(organism_label = fct_reorder(organism_label, n_detection_rows, .fun = sum, .desc = TRUE)) %>%
    arrange(factor(organism_type, levels = organism_type_levels), organism_label)

  ggplot(plot_data, aes(culture_month, n_detection_rows, fill = organism_label)) +
    geom_col(width = 25 * 24 * 60 * 60, color = "white", linewidth = 0.05) +
    scale_fill_manual(values = make_organism_palette(plot_data)) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    labs(title = title, x = NULL, y = "Positive organism detection rows", fill = NULL, caption = organism_palette_caption) +
    stacked_time_theme +
    guides(fill = guide_legend(ncol = 2))
}

plot_faceted_stacked_bar <- function(data, title) {
  plot_data <- data %>%
    group_by(culture_type, organism_label) %>%
    mutate(total_detection_rows = sum(n_detection_rows)) %>%
    ungroup() %>%
    mutate(organism_label = fct_reorder(organism_label, total_detection_rows, .desc = TRUE)) %>%
    arrange(factor(organism_type, levels = organism_type_levels), organism_label)

  ggplot(plot_data, aes(culture_month, n_detection_rows, fill = organism_label)) +
    geom_col(width = 25 * 24 * 60 * 60, color = "white", linewidth = 0.05) +
    facet_wrap(vars(culture_type), scales = "free_y", ncol = 2) +
    scale_fill_manual(values = make_organism_palette(plot_data)) +
    scale_x_datetime(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    labs(title = title, x = NULL, y = "Positive organism detection rows", fill = NULL, caption = organism_palette_caption) +
    stacked_time_theme +
    guides(fill = guide_legend(ncol = 2))
}

p_group_overall <- plot_overall_bar(
  organism_group_overall,
  "Predominant Organism Groups in Positive ICU Cultures"
)

p_category_overall <- plot_overall_bar(
  organism_category_overall,
  "Predominant Organisms in Positive ICU Cultures"
)

p_group_by_type <- plot_faceted_bar(
  organism_group_by_type,
  "Predominant Organism Groups by Culture Type"
)

p_category_by_type <- plot_faceted_bar(
  organism_category_by_type,
  "Predominant Organisms by Culture Type"
)

p_group_overall_time <- plot_overall_stacked_bar(
  monthly_organism_group_overall,
  "Monthly Predominant Organism Groups in Positive ICU Cultures"
)

p_category_overall_time <- plot_overall_stacked_bar(
  monthly_organism_category_overall,
  "Monthly Predominant Organisms in Positive ICU Cultures"
)

p_group_by_type_time <- plot_faceted_stacked_bar(
  monthly_organism_group_by_type,
  "Monthly Predominant Organism Groups by Culture Type"
)

p_category_by_type_time <- plot_faceted_stacked_bar(
  monthly_organism_category_by_type,
  "Monthly Predominant Organisms by Culture Type"
)

plot_paths <- c(
  organism_group_overall = file.path(out_dir, glue("positive_organism_group_overall_{site_name}_{stamp}.png")),
  organism_category_overall = file.path(out_dir, glue("positive_organism_category_overall_{site_name}_{stamp}.png")),
  organism_group_by_type = file.path(out_dir, glue("positive_organism_group_by_culture_type_{site_name}_{stamp}.png")),
  organism_category_by_type = file.path(out_dir, glue("positive_organism_category_by_culture_type_{site_name}_{stamp}.png")),
  monthly_organism_group_overall = file.path(out_dir, glue("monthly_positive_organism_group_overall_{site_name}_{stamp}.png")),
  monthly_organism_category_overall = file.path(out_dir, glue("monthly_positive_organism_category_overall_{site_name}_{stamp}.png")),
  monthly_organism_group_by_type = file.path(out_dir, glue("monthly_positive_organism_group_by_culture_type_{site_name}_{stamp}.png")),
  monthly_organism_category_by_type = file.path(out_dir, glue("monthly_positive_organism_category_by_culture_type_{site_name}_{stamp}.png"))
)

ggsave(plot_paths[["organism_group_overall"]], p_group_overall, width = 9, height = 7, dpi = 300)
ggsave(plot_paths[["organism_category_overall"]], p_category_overall, width = 9, height = 7, dpi = 300)
ggsave(plot_paths[["organism_group_by_type"]], p_group_by_type, width = 12, height = 12, dpi = 300)
ggsave(plot_paths[["organism_category_by_type"]], p_category_by_type, width = 12, height = 12, dpi = 300)
ggsave(plot_paths[["monthly_organism_group_overall"]], p_group_overall_time, width = 12, height = 8, dpi = 300)
ggsave(plot_paths[["monthly_organism_category_overall"]], p_category_overall_time, width = 12, height = 8, dpi = 300)
ggsave(plot_paths[["monthly_organism_group_by_type"]], p_group_by_type_time, width = 14, height = 12, dpi = 300)
ggsave(plot_paths[["monthly_organism_category_by_type"]], p_category_by_type_time, width = 14, height = 12, dpi = 300)

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
