# Shared site configuration and table readers for CLIF project scripts.

source("utils/config.R")

clif_site_name <- config_value(config, "site_name", env = "CLIF_SITE_NAME", default = "SITE")
clif_repo_path <- config_value(config, "repo", env = "CLIF_REPO", default = getwd())
clif_tables_path <- config_value(config, "tables_path", env = "CLIF_TABLES_PATH", required = TRUE)
clif_file_type <- tolower(config_value(config, "file_type", env = "CLIF_FILE_TYPE", default = "parquet"))

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
