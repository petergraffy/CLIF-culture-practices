# Configuration

Copy `config_template.json` to `config.json` and update it for the local environment.

Required fields:

1. `site_name`: short site label used in output filenames.
2. `repo`: absolute path to this repository.
3. `tables_path`: absolute path to the CLIF table directory.
4. `file_type`: CLIF table file type, usually `parquet`, `csv`, or `fst`.
5. `study_start_date`: first culture collection date to include, formatted as `YYYY-MM-DD`.
6. `study_end_date`: last culture collection date to include, formatted as `YYYY-MM-DD`.

The code locates CLIF tables recursively under `tables_path` and accepts filenames with or without the `clif_` prefix, as long as the base table name is unique. For example, `clif_hospitalization.parquet` and `hospitalization.parquet` are both valid.

Common environment variable overrides:

1. `CLIF_CONFIG_PATH`: alternate config JSON path.
2. `CLIF_SITE_NAME`: override `site_name`.
3. `CLIF_REPO`: override `repo`.
4. `CLIF_TABLES_PATH`: override `tables_path`.
5. `CLIF_FILE_TYPE`: override `file_type`; use `auto` to scan `csv`, `parquet`, and `fst`.
6. `STUDY_START_DATE` and `STUDY_END_DATE`: override the configured study window.

All standard script outputs are written under `<repo>/output/`. Private row-level intermediates are written under `<repo>/data/intermediate/` and are read automatically by downstream scripts.

The `.gitignore` file prevents `config/config.json` from being pushed to GitHub. Keep site-specific paths and credentials local.
