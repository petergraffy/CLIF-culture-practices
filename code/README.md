# Running the culture-practices pipeline

Copy `config/config_template.json` to ignored `config/config.json` and set site, repository/table paths, and file type. The study team maintains the shared inclusive dates in `utils/study_settings.R`, currently January 1, 2018–December 31, 2024. Culture coverage is assumed complete within that window. From the repository root, install the pinned packages before the first pipeline run, then run the pipeline:

```sh
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript code/00_run_pipeline.R
```

The restore command is a separate first-time setup step. The pipeline does not install missing packages; `.Rprofile` activates the project environment automatically on subsequent runs. Restore again when `renv.lock` changes or preflight reports a dependency mismatch.

This is the recommended site workflow. Every run has an isolated `output/runs/<run_id>/` aggregate tree and `data/intermediate/runs/<run_id>/` private tree. Source rows and exact identifiers/timestamps remain private. The manifest records code/config/mapping hashes, input metadata, packages, and completion status. Review QC and the export audit before sharing.

## Manual runs

Use one unique `CLIF_RUN_ID` for all scripts in a manual run:

```sh
export CLIF_RUN_ID=SITE_20261002_review
Rscript code/01_identify_icu_culture_cohort.R
Rscript code/10_quality_checks.R
Rscript code/02_plot_culture_time_series.R
Rscript code/04_plot_positive_organisms.R
Rscript code/05_culture_rates_per_icu_admission.R
Rscript code/06_icu_day_denominators_and_timing.R
Rscript code/08_organism_trends.R
Rscript code/09_susceptibility_trends.R
Rscript code/07_prepare_site_exports.R
```

Without `CLIF_RUN_ID`, standalone scripts use the original output/intermediate folders. Legacy intermediates lacking `culture_event_id` must be rebuilt. Explicit `ICU_CULTURE_ROWS_PATH`/`ICU_CULTURE_EVENTS_PATH` overrides are for debugging and bypass normal run selection; the automated runner rejects them.

Susceptibility analysis skips cleanly when the optional table is absent. Sites with that table must provide the standard mCIDE fields. `AST_MIN_TESTING_FRACTION` controls the minimum interpretable testing coverage for detection-rate models (default 0.5). Pair-specific tested-fraction models use interpretable tests and export coverage separately.

Other overrides include `CLIF_CONFIG_PATH`, `CLIF_SITE_NAME`, `CLIF_TABLES_PATH`, `CLIF_FILE_TYPE`, `CLIF_REPO`, and plot-specific top-N controls. Verify plots cover the intended full analysis window. `WRITE_ROW_LEVEL_INTERMEDIATES=false` is supported for cohort-only exports; the full pipeline requires private intermediates.

Definitions, model interpretation, mapping overrides, QC limitations, and AST testing rules are in [analysis_definitions.md](../docs/analysis_definitions.md).

## Verification

```sh
Rscript tests/test_preflight.R
Rscript tests/test_config.R
Rscript tests/test_core.R
Rscript tests/run_integration.R
```

The synthetic tests exercise stay/event boundaries, carry-in denominators, missing/mixed result categories, date-specific mappings, duplicate/conflicting AST reports, mCIDE field requirements, and count/proportion GAMs. Real susceptibility ETL still needs a buddy test at a site with both tables.

`11_pool_site_trends.R` is a separate central workflow using completed aggregate runs and `config/pooling_sites.csv`; it is not invoked by the site runner. See `docs/pooling.md`.

Dependency setup uses the committed `.Rprofile` and `renv/activate.R` bootstrap (renv 1.1.5). Restore with `Rscript -e 'renv::restore(prompt = FALSE)'` from the repository root. The site runner checks dependencies and table schemas before analysis, and writes aggregate date availability to `provenance/preflight_source_availability.csv` in the completed run.
