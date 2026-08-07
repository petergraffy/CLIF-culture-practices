# Shared site configuration and table readers for CLIF project scripts.

source("utils/config.R")

clif_site_name <- config_value(config, "site_name", env = "CLIF_SITE_NAME", default = "SITE")
clif_repo_path <- normalizePath(
  config_value(config, "repo", env = "CLIF_REPO", default = getwd()),
  winslash = "/",
  mustWork = FALSE
)
clif_tables_path <- config_value(config, "tables_path", env = "CLIF_TABLES_PATH", required = TRUE)
clif_file_type <- tolower(config_value(config, "file_type", env = "CLIF_FILE_TYPE", default = "parquet"))

project_path <- function(...) {
  file.path(clif_repo_path, ...)
}

project_output_path <- function(...) {
  project_path("output", ...)
}

project_output_dir <- function(...) {
  path <- project_output_path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

project_intermediate_path <- function(...) {
  project_path("data", "intermediate", ...)
}

project_intermediate_dir <- function(...) {
  path <- project_intermediate_path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

latest_project_intermediate_file <- function(pattern, ...) {
  path <- project_intermediate_path(...)
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files found in ", path, " matching pattern: ", pattern)
  }
  files[which.max(file.info(files)$mtime)]
}

clif_specimen_type_palette <- c(
  "Blood buffy" = "#FB8072",
  "Nasopharynx upperairway" = "#80B1D3",
  "Respiratory tract" = "#BC80BD",
  "Genito urinary tract" = "#FDB462",
  "Respiratory tract lower" = "#B3DE69",
  "Meninges csf" = "#8DD3C7",
  "Pleural cavity fluid" = "#BEBADA",
  "Peritoneum" = "#FCCDE5",
  "Woundsite" = "#CCEBC5",
  "Other unspecified" = "#FFFFB3",
  "Catheter tip" = "#FFED6F",
  "Joints" = "#A6CEE3",
  "Feces stool" = "#B2DF8A",
  "Cardiac" = "#FDBF6F",
  "Gallbladder billary pancreas" = "#CAB2D6",
  "Gastrointestinal tract" = "#FFFF99",
  "Oropharynx tongue oralcavity" = "#1F78B4",
  "Skin" = "#33A02C",
  "Ears" = "#E31A1C",
  "Bone cortex" = "#FF7F00",
  "Eyes" = "#6A3D9A",
  "Bone marrow" = "#B15928",
  "Kidneys renal pelvis ureters bladder" = "#A1D99B",
  "Genital area" = "#9E9AC8",
  "Vagina" = "#FDD0A2",
  "Stomach" = "#9ECAE1",
  "Lymph nodes" = "#F2B6C6",
  "Fallopians uterus cervix" = "#C49C94",
  "Other" = "#BDBDBD",
  "No ICU culture" = "#C7C7C7"
)

clif_complete_specimen_palette <- function(specimen_levels, include_no_icu = FALSE) {
  specimen_levels <- as.character(specimen_levels)
  if (!include_no_icu) {
    base_palette <- clif_specimen_type_palette[names(clif_specimen_type_palette) != "No ICU culture"]
  } else {
    base_palette <- clif_specimen_type_palette
  }

  missing_levels <- setdiff(specimen_levels, names(base_palette))
  if (length(missing_levels) > 0) {
    fallback_colors <- grDevices::hcl.colors(
      n = length(missing_levels),
      palette = "Dark 3",
      alpha = 1
    )
    base_palette <- c(base_palette, stats::setNames(fallback_colors, missing_levels))
  }

  base_palette[names(base_palette) %in% specimen_levels]
}

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  out <- switch(
    ext,
    "csv" = readr::read_csv(path, show_col_types = FALSE),
    "parquet" = arrow::read_parquet(path),
    "fst" = {
      if (!requireNamespace("fst", quietly = TRUE)) stop("Package 'fst' is required to read fst files.")
      fst::read_fst(path)
    },
    stop("Unsupported extension: ", ext, " for path: ", path)
  )

  janitor::clean_names(out)
}

find_table_path <- function(tbl_base, tables_path = clif_tables_path, file_type = clif_file_type, required = TRUE) {
  wanted <- tolower(tbl_base)
  if (!startsWith(wanted, "clif_")) wanted <- paste0("clif_", wanted)

  allowed_ext <- if (!is.null(file_type) && nzchar(file_type) && file_type != "auto") {
    tolower(file_type)
  } else {
    c("csv", "parquet", "fst")
  }

  files <- list.files(tables_path, full.names = TRUE, recursive = TRUE)
  files <- files[tolower(tools::file_ext(files)) %in% allowed_ext]
  base <- tolower(tools::file_path_sans_ext(basename(files)))
  base_norm <- ifelse(startsWith(base, "clif_"), base, paste0("clif_", base))
  hit <- files[base_norm == wanted]

  if (length(hit) == 1) return(hit)
  if (required) {
    stop("Could not uniquely locate ", wanted, " in ", tables_path, ". Matches: ", length(hit))
  }

  NA_character_
}

read_tbl <- function(tbl_base, required = TRUE) {
  path <- find_table_path(tbl_base, required = required)
  if (is.na(path)) return(NULL)
  read_any(path)
}
