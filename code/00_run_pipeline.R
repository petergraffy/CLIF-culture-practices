# Single command creates an isolated, reproducible run; downstream scripts never select another run.
source("utils/preflight.R")
preflight_packages()
suppressPackageStartupMessages(library(jsonlite))
source("utils/config.R")
if (nzchar(Sys.getenv("ICU_CULTURE_ROWS_PATH", "")) || nzchar(Sys.getenv("ICU_CULTURE_EVENTS_PATH", ""))) stop("Remove explicit intermediate path overrides for an isolated pipeline run.")
if (tolower(Sys.getenv("WRITE_ROW_LEVEL_INTERMEDIATES", "true")) %in% c("false", "0", "no", "n")) stop("The full pipeline requires private row-level intermediates.")
run_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
Sys.setenv(CLIF_RUN_ID = run_id)
source("utils/clif_io.R")
preflight <- run_site_preflight()
manifest_dir <- project_output_dir("provenance")
readr::write_csv(preflight$summary,file.path(manifest_dir,"preflight_source_availability.csv"))
scripts <- c("01_identify_icu_culture_cohort.R", "10_quality_checks.R", "02_plot_culture_time_series.R", "04_plot_positive_organisms.R", "05_culture_rates_per_icu_admission.R", "06_icu_day_denominators_and_timing.R", "08_organism_trends.R", "09_susceptibility_trends.R", "11_ase_stratified_analysis.R", "07_prepare_site_exports.R")
files <- c(".Rprofile", "renv/activate.R", "renv/settings.json", list.files("code", full.names = TRUE, pattern = "[.]R$"), list.files("utils", full.names = TRUE, recursive = TRUE, pattern = "[.](R|sql|json|txt)$"), list.files("config", full.names = TRUE, pattern = "[.]csv$"), list.files("config/mcide", full.names = TRUE), if (file.exists("renv.lock")) "renv.lock" else character())
config_path <- Sys.getenv("CLIF_CONFIG_PATH", "config/config.json")
manifest <- list(run_id = run_id, site_name = clif_site_name, started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config_md5 = unname(tools::md5sum(config_path)), study_start_date = study_settings$study_start_date, study_end_date = study_settings$study_end_date, code_and_mapping_md5 = as.list(tools::md5sum(files)), R_version = R.version.string, packages = as.list(setNames(vapply(c("dplyr", "tidyr", "readr", "lubridate", "mgcv", "ggplot2"), function(p) as.character(packageVersion(p)), character(1)), c("dplyr", "tidyr", "readr", "lubridate", "mgcv", "ggplot2"))), environment_overrides = as.list(Sys.getenv()[grepl("^(CLIF_|STUDY_|PLOT_|TOP_N_|AST_|TIMING_|WRITE_ROW_)", names(Sys.getenv()))]), analysis_status = "running")
# Exact paths and environment values stay private.
private_manifest <- file.path(project_intermediate_dir("provenance"), "run_manifest.json")
write_json(manifest, private_manifest, pretty = TRUE, auto_unbox = TRUE)
# Preserve the exact uncommitted analysis source behind the manifest hashes.
for (file in files) {
  snapshot_path <- file.path(dirname(private_manifest), "source_snapshot", file)
  dir.create(dirname(snapshot_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file, snapshot_path, overwrite = FALSE)) stop("Could not snapshot analysis source: ", file)
}
file.copy(config_path, file.path(dirname(private_manifest), "config_used.json"))
input_files <- vapply(c("hospitalization", "adt", "microbiology_culture", "patient", "labs", "medication_admin_intermittent", "medication_admin_continuous", "respiratory_support", "hospital_diagnosis"), find_table_path, character(1), required = FALSE)
input_files <- input_files[!is.na(input_files)]
optional <- find_table_path("microbiology_susceptibility", required = FALSE)
if (!is.na(optional)) input_files <- c(input_files, microbiology_susceptibility = optional)
readr::write_csv(data.frame(table = names(input_files), file_basename = basename(input_files), size_bytes = file.info(input_files)$size, modified_utc = format(file.info(input_files)$mtime, tz = "UTC", usetz = TRUE)), file.path(manifest_dir, "input_file_metadata.csv"))
expected_hashes <- tools::md5sum(files)
expected_config_hash <- tools::md5sum(config_path)
validate_run_inputs <- function() {
  if (!identical(tools::md5sum(files), expected_hashes) || !identical(tools::md5sum(config_path), expected_config_hash)) stop("Code, mappings, dependencies, or configuration changed during the run; start a new isolated run.")
}
status <- tryCatch({
  for (script in scripts) {
    validate_run_inputs()
    message("Running ", script)
    exit_status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path("code", script)))
    if (exit_status != 0) stop("Pipeline failed at ", script)
  }
  validate_run_inputs()
  "completed"
}, error = function(e) {manifest$error <<- conditionMessage(e); "failed"})
manifest$analysis_status <- status; manifest$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
write_json(manifest, private_manifest, pretty = TRUE, auto_unbox = TRUE)
public <- manifest; public$environment_overrides <- NULL; public$error <- NULL
write_json(public, file.path(manifest_dir, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
if (status != "completed") stop(manifest$error)
message("Completed run ", run_id, "; aggregate exports: ", project_output_path())
