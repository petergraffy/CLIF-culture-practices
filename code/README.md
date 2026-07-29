# Running the CLIF Culture Practices Pipeline

These scripts are intended to be run by each CLIF site from the repository root after creating a local `config/config.json`.

## One-Time Setup

1. Copy `config/config_template.json` to `config/config.json`.
2. Edit `config/config.json`:
   - `site_name`: short site label for output filenames.
   - `repo`: absolute path to this repository on the local machine.
   - `tables_path`: absolute path to the local CLIF table directory.
   - `file_type`: `parquet`, `csv`, `fst`, or `auto`.
   - `study_start_date` and `study_end_date`: inclusive calendar window.

Do not commit `config/config.json`; it contains local site paths.

## Run Order

Run these commands from the repository root:

```sh
Rscript code/01_identify_icu_culture_cohort.R
Rscript code/02_plot_culture_time_series.R
Rscript code/04_plot_positive_organisms.R
Rscript code/05_culture_rates_per_icu_admission.R
Rscript code/06_icu_day_denominators_and_timing.R
Rscript code/08_organism_trends.R
Rscript code/07_prepare_site_exports.R
```

`code/01_identify_icu_culture_cohort.R` must run first. It writes private row-level intermediates under:

```text
<repo>/data/intermediate/cohort/
```

Downstream scripts automatically read the latest matching intermediate file from that folder. Sites should not provide intermediate paths manually for the standard pipeline.

## Outputs

Shareable aggregate CSVs, figures, and the export manifest are written under:

```text
<repo>/output/
```

The final script, `code/07_prepare_site_exports.R`, creates a manifest and privacy audit under `output/site_exports/`. Confirm this script completes without errors before sharing or pooling site outputs.

## Private Data

Do not share or upload `data/intermediate/`. It can include row-level identifiers and exact timestamps. Only files under `output/` are intended for cross-site comparison after the privacy audit passes.

## Optional Overrides

The config file should be enough for ordinary site runs. Environment variables are still available for troubleshooting or one-off reruns:

- `CLIF_CONFIG_PATH`: alternate config JSON path.
- `CLIF_SITE_NAME`, `CLIF_TABLES_PATH`, `CLIF_FILE_TYPE`, `CLIF_REPO`: override matching config fields.
- `STUDY_START_DATE`, `STUDY_END_DATE`: override the study window.
- `PLOT_START_DATE`, `PLOT_END_DATE`: override only plot windows.
- `ICU_CULTURE_EVENTS_PATH`, `ICU_CULTURE_ROWS_PATH`: advanced/debug-only explicit intermediate inputs.
