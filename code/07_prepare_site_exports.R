# ================================================================================================
# Prepare Site Exports
#
# Purpose:
#   Create a manifest of shareable site-level CSV outputs and fail fast if output/ contains
#   obvious row-level identifiers or exact timestamp columns.
#
# Governance:
#   output/ is the shareable export tree. PHI-bearing intermediates belong in ignored
#   data/intermediate/ and should not be pooled across sites.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(readr)
  library(stringr)
  library(tibble)
})

source("utils/clif_io.R")

site_name <- clif_site_name
output_dir <- Sys.getenv("SITE_EXPORT_OUTPUT_DIR", unset = "output")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
manifest_dir <- file.path(output_dir, "site_exports")
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

disallowed_columns <- c(
  "patient_id",
  "hospitalization_id",
  "icu_admission_id",
  "icu_interval_id",
  "microbiology_row_id",
  "organism_id",
  "admission_dttm",
  "discharge_dttm",
  "icu_in_dttm",
  "icu_out_dttm",
  "icu_in_dttm_clipped",
  "icu_out_dttm_clipped",
  "order_dttm",
  "collect_dttm",
  "result_dttm",
  "first_culture_dttm"
)

csv_files <- list.files(output_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
csv_files <- csv_files[!str_detect(csv_files, fixed(file.path("site_exports", "")))]

inspect_csv <- function(path) {
  header <- names(read_csv(path, n_max = 0, show_col_types = FALSE))
  output_abs <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  path_abs <- normalizePath(path, winslash = "/", mustWork = FALSE)
  tibble(
    site_name = site_name,
    export_relative_path = str_remove(path_abs, paste0("^", output_abs, "/?")),
    file_size_bytes = file.info(path)$size,
    n_rows = length(count.fields(path, sep = ",")) - 1L,
    n_columns = length(header),
    columns = paste(header, collapse = "; "),
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}

manifest <- if (length(csv_files) > 0) {
  bind_rows(lapply(csv_files, inspect_csv))
} else {
  tibble(
    site_name = site_name,
    export_relative_path = character(),
    file_size_bytes = numeric(),
    n_rows = integer(),
    n_columns = integer(),
    columns = character(),
    created_at = character()
  )
}

privacy_audit <- manifest %>%
  tidyr::separate_longer_delim(columns, delim = "; ") %>%
  mutate(
    lower_column = str_to_lower(columns),
    disallowed = lower_column %in% disallowed_columns | str_detect(lower_column, "(^|_)datetime$")
  ) %>%
  filter(disallowed) %>%
  select(site_name, export_relative_path, disallowed_column = columns)

manifest_path <- file.path(manifest_dir, glue("site_export_manifest_{site_name}_{stamp}.csv"))
audit_path <- file.path(manifest_dir, glue("site_export_privacy_audit_{site_name}_{stamp}.csv"))

write_csv(manifest, manifest_path)
write_csv(privacy_audit, audit_path)

message("Wrote site export manifest: ", manifest_path)
message("Wrote site export privacy audit: ", audit_path)

if (nrow(privacy_audit) > 0) {
  print(privacy_audit, n = Inf)
  stop("Potential PHI-bearing columns found in output/. Move row-level exports to data/intermediate/ before pooling.")
}

message("No disallowed row-level identifier or exact timestamp columns found in output CSV exports.")
