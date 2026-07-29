# ================================================================================================
# Organism Detection Trend Screen
#
# Question:
#   Which organisms or targeted resistance-related organism labels are increasing over calendar time?
#
# Denominators:
#   1. All ICU admissions, defined as merged ICU ADT intervals and counted by ICU admission month.
#   2. All ICU days, allocated to calendar months from merged ICU ADT interval overlap time.
#
# Numerator:
#   Positive ICU culture detection events. Events are collapsed by ICU culture event and organism,
#   then counted monthly.
#
# Note:
#   When CLIF microbiology_susceptibility is available, MRSA/VRE/CRE phenotypes are derived from
#   resistant antimicrobial susceptibility results. When it is unavailable, target labels fall back
#   to explicit organism text only and the aggregate source summary records that limitation.
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

first_existing_col <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

normalise_sir <- function(x) {
  x_clean <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x_clean %in% c("r", "resistant", "res") ~ "resistant",
    x_clean %in% c("i", "intermediate", "intermediate/resistant", "non-susceptible", "nonsusceptible") ~ "intermediate",
    x_clean %in% c("s", "susceptible", "sus") ~ "susceptible",
    TRUE ~ x_clean
  )
}

classify_organism_type <- function(x) {
  x_clean <- str_to_lower(coalesce(x, ""))
  case_when(
    str_detect(x_clean, "candida|yeast|fung|aspergillus|cryptococcus|mold|mould|saccharomyces|fusarium|mucor|rhizopus|pneumocystis|torulopsis") ~ "Fungi/yeast",
    str_detect(x_clean, "mycobacter|tuberculosis|\\bafb\\b|acid fast") ~ "Mycobacteria/AFB",
    str_detect(x_clean, "anaerob|bacteroides|clostrid|prevotella|fusobacter|cutibacter|propionibacter|propionbacterium|leptotrichia|clostridioides") ~ "Anaerobes",
    str_detect(x_clean, "amebiasis|cryptosporidium|echinoco|giardia|protozo|toxoplasma|trichomonas|parasite") ~ "Parasite/protozoa",
    str_detect(x_clean, "adenovirus|cytomegalovirus|enterovirus|epstein|hepatitis|herpes|hhv|hiv|influenza|measles|mumps|papovavirus|parainfluenza|polyomavirus|respiratory_syncytial|rsv|rhinovirus|rotavirus|rubella|virus|viral|covid|sars|cmv") ~ "Virus",
    str_detect(x_clean, "borrelia|chlamydia|coxiella|leptospira|mycoplasma|rickettsia|treponema") ~ "Atypical/other bacteria",
    str_detect(x_clean, "staphylococcus|streptococcus|enterococcus|bacillus|corynebacter|lactobacillus|listeria|leuconostoc|micrococcus|nocardia|rhodococcus|stomatococcus|mrsa|vre|gram_positive|gram positive|gpc|coag_pos|coag_neg|coagneg") ~ "Gram positive bacteria",
    str_detect(x_clean, "acinetobacter|agrobacterium|alcaligenes|branhamelia|moraxella|pseudomonas|stenotrophomonas|xanthomonas|klebsiella|enterobacter|escherichia|serratia|haemophilus|citrobacter|proteus|neisseria|salmonella|shigella|campylobacter|burkholderia|cepacia|legionella|flavimonas|flavobacterium|helicobacter|methylobacterium|vibrio|esbl|cre|carbapenem|gram_negative|gram negative|gnr") ~ "Gram negative bacteria",
    str_detect(x_clean, "bacteria|gram|cocci|bacilli|rods|flora") ~ "Other bacteria",
    TRUE ~ "Other/unspecified"
  )
}

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

  icu_admissions <- adt %>%
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
      icu_out_dttm = max(icu_out_dttm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      icu_in_dttm = if_else(!is.na(study_start_dttm) & icu_in_dttm < study_start_dttm, study_start_dttm, icu_in_dttm),
      icu_out_dttm = if_else(!is.na(study_end_dttm) & icu_out_dttm > study_end_dttm, study_end_dttm, icu_out_dttm)
    ) %>%
    filter(icu_out_dttm > icu_in_dttm)

  monthly_admissions <- icu_admissions %>%
    mutate(calendar_month = floor_date(icu_in_dttm, "month")) %>%
    count(calendar_month, name = "n_icu_admissions") %>%
    complete(calendar_month = month_seq, fill = list(n_icu_admissions = 0L)) %>%
    arrange(calendar_month)

  month_windows <- tibble(
    calendar_month = month_seq,
    month_start = month_seq,
    month_end = month_seq %m+% months(1)
  )

  monthly_days <- icu_admissions %>%
    tidyr::crossing(month_windows) %>%
    mutate(
      overlap_start = pmax(icu_in_dttm, month_start),
      overlap_end = pmin(icu_out_dttm, month_end),
      overlap_days = as.numeric(difftime(overlap_end, overlap_start, units = "hours")) / 24
    ) %>%
    filter(overlap_days > 0) %>%
    group_by(calendar_month) %>%
    summarise(n_icu_days = sum(overlap_days), .groups = "drop") %>%
    complete(calendar_month = month_seq, fill = list(n_icu_days = 0)) %>%
    arrange(calendar_month)

  monthly_admissions %>%
    left_join(monthly_days, by = "calendar_month") %>%
    mutate(n_icu_days = coalesce(n_icu_days, 0))
}

target_organism_patterns <- tibble::tribble(
  ~target_label, ~pattern,
  "MRSA", "mrsa|methicillin[ _-]*resistant.*staph|staph.*methicillin[ _-]*resistant",
  "Staphylococcus aureus", "staphylococcus[_ ]aureus",
  "VRE", "\\bvre\\b|vancomycin[ _-]*resistant.*enterococcus|enterococcus.*vancomycin[ _-]*resistant",
  "Enterococcus faecium", "enterococcus[_ ]faecium",
  "ESBL", "\\besbl\\b|extended[ _-]*spectrum",
  "CRE", "\\bcre\\b|carbapenem[ _-]*resistant|\\bkpc\\b|\\bndm\\b",
  "Pseudomonas aeruginosa", "pseudomonas[_ ]aeruginosa",
  "Klebsiella pneumoniae", "klebsiella[_ ]pneumoniae",
  "Escherichia coli", "escherichia[_ ]coli|\\be[._ ]?coli\\b",
  "Candida auris", "candida[_ ]auris",
  "Clostridioides difficile", "clostridioides[_ ]difficile|clostridium[_ ]difficile"
) %>%
  mutate(
    organism_type = factor(classify_organism_type(target_label), levels = organism_type_levels)
)

read_susceptibility_table <- function() {
  path <- find_table_path("microbiology_susceptibility", required = FALSE)
  if (is.na(path)) {
    message("CLIF microbiology_susceptibility table not found; target resistance phenotypes will use organism text fallback only.")
    return(NULL)
  }

  message("Reading CLIF microbiology_susceptibility: ", path)
  read_tbl("microbiology_susceptibility")
}

build_susceptibility_phenotypes <- function(rows) {
  susceptibility <- read_susceptibility_table()
  if (is.null(susceptibility) || nrow(susceptibility) == 0) {
    return(tibble(
      target_label = character(),
      detection_event_id = character(),
      calendar_month = as.POSIXct(character()),
      phenotype_source = character()
    ))
  }

  organism_id_col <- first_existing_col(susceptibility, c("organism_id", "microbiology_organism_id"))
  antimicrobial_col <- first_existing_col(
    susceptibility,
    c("antimicrobial_name", "antimicrobial_category", "antibiotic_name", "antibiotic", "drug_name", "susceptibility_name")
  )
  interpretation_col <- first_existing_col(
    susceptibility,
    c("susceptibility_interpretation", "interpretation", "susceptibility_result", "result", "result_category")
  )

  if (is.na(organism_id_col) || is.na(antimicrobial_col) || is.na(interpretation_col)) {
    warning(
      "microbiology_susceptibility is present but missing expected columns. Found columns: ",
      paste(names(susceptibility), collapse = ", "),
      ". Expected organism_id plus antimicrobial and interpretation/result fields."
    )
    return(tibble(
      target_label = character(),
      detection_event_id = character(),
      calendar_month = as.POSIXct(character()),
      phenotype_source = character()
    ))
  }

  susceptibility_resistant <- susceptibility %>%
    transmute(
      organism_id = as.character(.data[[organism_id_col]]),
      antimicrobial = str_to_lower(str_trim(as.character(.data[[antimicrobial_col]]))),
      interpretation = normalise_sir(.data[[interpretation_col]])
    ) %>%
    filter(!is.na(organism_id), nzchar(organism_id), interpretation == "resistant")

  if (nrow(susceptibility_resistant) == 0) {
    return(tibble(
      target_label = character(),
      detection_event_id = character(),
      calendar_month = as.POSIXct(character()),
      phenotype_source = character()
    ))
  }

  rows_for_ast <- rows %>%
    mutate(organism_id = as.character(organism_id)) %>%
    filter(!is.na(organism_id), nzchar(organism_id)) %>%
    select(detection_event_id, calendar_month, organism_id, organism_text, organism_name, organism_category, organism_group)

  rows_for_ast %>%
    inner_join(susceptibility_resistant, by = "organism_id", relationship = "many-to-many") %>%
    mutate(
      target_label = case_when(
        str_detect(organism_text, "staphylococcus[_ ]aureus") &
          str_detect(antimicrobial, "oxacillin|cefoxitin|methicillin") ~ "MRSA",
        str_detect(organism_text, "enterococcus") &
          str_detect(antimicrobial, "vancomycin") ~ "VRE",
        str_detect(organism_text, "escherichia|klebsiella|enterobacter|citrobacter|serratia|proteus") &
          str_detect(antimicrobial, "ertapenem|imipenem|meropenem|doripenem") ~ "CRE",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(target_label)) %>%
    transmute(
      target_label,
      detection_event_id,
      calendar_month,
      phenotype_source = "susceptibility"
    ) %>%
    distinct()
}

fit_poisson_trend <- function(data, denominator_var) {
  model_data <- data %>%
    filter(!is.na(.data[[denominator_var]]), .data[[denominator_var]] > 0) %>%
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
    offset = log(.subset2(model_data, denominator_var)),
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
    organism_id = as.character(organism_id),
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
monthly_icu_denominators <- build_monthly_icu_denominators(month_seq, study_start_dttm, study_end_dttm)

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

target_susceptibility_detections <- build_susceptibility_phenotypes(rows)

target_detections <- bind_rows(
  target_susceptibility_detections,
  target_text_detections
) %>%
  arrange(
    target_label,
    detection_event_id,
    desc(phenotype_source == "susceptibility")
  ) %>%
  distinct(target_label, detection_event_id, calendar_month, .keep_all = TRUE)

target_detection_source_summary <- target_detections %>%
  count(target_label, phenotype_source, name = "n_detection_events") %>%
  complete(
    target_label = target_organism_patterns$target_label,
    phenotype_source = c("susceptibility", "organism_text"),
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
    group_modify(~ fit_poisson_trend(.x, denominator_var)) %>%
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
      trend_direction = case_when(
        total_detection_events == 0 ~ zero_label,
        is.na(annual_percent_change) | is.na(p_value) ~ "Not estimated",
        p_value < 0.05 & annual_percent_change > 0 ~ "Increasing",
        p_value < 0.05 & annual_percent_change < 0 ~ "Decreasing",
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

increasing_labels <- trend_summary_top_admissions %>%
  filter(total_detection_events >= 25, !is.na(annual_percent_change)) %>%
  arrange(desc(annual_percent_change)) %>%
  slice_head(n = plot_n_increasing) %>%
  pull(organism_label)

target_plot_labels <- target_monthly_rates %>%
  group_by(organism_label) %>%
  summarise(total_detection_events = max(total_detection_events, na.rm = TRUE), .groups = "drop") %>%
  filter(total_detection_events > 0 | organism_label %in% c("MRSA", "Staphylococcus aureus")) %>%
  arrange(desc(total_detection_events), organism_label) %>%
  pull(organism_label)

prepare_plot_data <- function(data, plot_labels, rate_var) {
  data %>%
    filter(organism_label %in% plot_labels) %>%
    mutate(
      plot_rate = .data[[rate_var]],
      organism_type = factor(organism_type, levels = organism_type_levels)
    ) %>%
    arrange(organism_label, calendar_month) %>%
    group_by(organism_label) %>%
    mutate(rate_rolling_6mo = rolling_mean_trailing(plot_rate, 6)) %>%
    ungroup() %>%
    mutate(organism_label = factor(organism_label, levels = plot_labels))
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
  available_taxonomy_palette <- taxonomy_palette[names(taxonomy_palette) %in% unique(as.character(data$organism_type))]

  ggplot(data, aes(calendar_month, plot_rate)) +
    geom_col(aes(fill = organism_type), width = 25 * 24 * 60 * 60) +
    geom_line(aes(y = rate_rolling_6mo), color = "black", linewidth = 0.65, na.rm = TRUE) +
    facet_wrap(vars(organism_label), scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = available_taxonomy_palette, drop = FALSE) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = comma, limits = c(0, NA)) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      fill = "Taxonomy"
    ) +
    theme_trends
}

out_dir <- file.path("output", "organism_trends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
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
  increasing_plot_admissions = file.path(out_dir, glue("monthly_fastest_increasing_organism_detection_rates_per_100_icu_admissions_{site_name}_{stamp}.png")),
  increasing_plot_icu_days = file.path(out_dir, glue("monthly_fastest_increasing_organism_detection_rates_per_100_icu_days_{site_name}_{stamp}.png")),
  target_plot_admissions = file.path(out_dir, glue("monthly_target_organism_detection_rates_per_100_icu_admissions_{site_name}_{stamp}.png")),
  target_plot_icu_days = file.path(out_dir, glue("monthly_target_organism_detection_rates_per_100_icu_days_{site_name}_{stamp}.png"))
)

write_csv(trend_summary_top_admissions, paths[["trend_summary_top_admissions"]])
write_csv(trend_summary_top_icu_days, paths[["trend_summary_top_icu_days"]])
write_csv(trend_summary_targets_admissions, paths[["trend_summary_targets_admissions"]])
write_csv(trend_summary_targets_icu_days, paths[["trend_summary_targets_icu_days"]])
write_csv(target_detection_source_summary, paths[["target_detection_source_summary"]])
write_csv(top_monthly_rates, paths[["monthly_top_rates"]])
write_csv(target_monthly_rates, paths[["monthly_target_rates"]])
write_csv(monthly_icu_denominators, paths[["monthly_icu_denominators"]])

p_increasing_admissions <- plot_trend_facets(
  plot_increasing_admissions_data,
  "Fastest Increasing Organism Detection Rates per 100 ICU Admissions",
  "Detection events per 100 ICU admissions",
  ncol = 3
)
p_increasing_icu_days <- plot_trend_facets(
  plot_increasing_icu_days_data,
  "Fastest Increasing Organism Detection Rates per 100 ICU Days",
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
  "monthly_target_rates"
)])
message("Wrote plots:")
print(paths[c(
  "increasing_plot_admissions",
  "increasing_plot_icu_days",
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

message("Target organism trend summary:")
print(
  trend_summary_targets_admissions %>%
    select(organism_label, total_detection_events, annual_percent_change, p_value, trend_direction) %>%
    arrange(desc(total_detection_events)),
  n = Inf,
  width = Inf
)
