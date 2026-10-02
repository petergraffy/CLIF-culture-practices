# Configuration

From the repository root, copy `config/config_template.json` to `config/config.json` and edit the four local fields below. Then install the pinned packages before the first pipeline run:

```sh
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript code/00_run_pipeline.R
```

Package installation is a separate first-time step; the pipeline only checks dependencies. The project environment activates automatically for subsequent runs from the repository root. Restore again after `renv.lock` changes or if preflight reports missing packages/version mismatches.

Required fields:

1. `site_name`: short site label used in output filenames.
2. `repo`: absolute path to this repository.
3. `tables_path`: absolute path to the CLIF table directory.
4. `file_type`: CLIF table file type, usually `parquet`, `csv`, or `fst`.

The code locates CLIF tables recursively under `tables_path` and accepts filenames with or without the `clif_` prefix, as long as the base table name is unique. For example, `clif_hospitalization.parquet` and `hospitalization.parquet` are both valid.

Common environment variable overrides:

1. `CLIF_CONFIG_PATH`: alternate config JSON path.
2. `CLIF_SITE_NAME`: override `site_name`.
3. `CLIF_REPO`: override `repo`.
4. `CLIF_TABLES_PATH`: override `tables_path`.
5. `CLIF_FILE_TYPE`: override `file_type`; use `auto` to scan `csv`, `parquet`, and `fst`.

All standard script outputs are written under `<repo>/output/`. Private row-level intermediates are written under `<repo>/data/intermediate/` and are read automatically by downstream scripts.

The `.gitignore` file prevents `config/config.json` from being pushed to GitHub. Keep site-specific paths and credentials local.

## Optional analyses and mapping corrections

Susceptibility uses canonical CLIF mCIDE fields and the pinned reference CSVs in `mcide/`. Sites without the table skip this analysis. Set `AST_MIN_TESTING_FRACTION` (default 0.5) to control interpretable testing / all observed organism isolates for detection-rate models. `AST_MIN_LINKAGE_FRACTION` (default 0.9) controls linkable IDs / all observed organism isolates for both rate and fraction models.

The study team maintains the shared inclusive analysis dates in `utils/study_settings.R` (currently January 1, 2018–December 31, 2024). All analyses and plots use those dates; sites configure only their label, local paths, and file type. Culture coverage is assumed complete within the shared window. Months without observed specimen-source activity and months with incomplete organism identification remain excluded from zero-organism AST inference. See the [buddy-testing guide](../docs/buddy_testing.md).

`specimen_category_overrides.csv` is initially empty. Add only source-verified exact fluid-name/original-category/date-range corrections, with a documented reason. No mapping is inferred from an order name alone.

Use `Rscript code/00_run_pipeline.R` for isolated run directories and manifests. See [analysis definitions](../docs/analysis_definitions.md).

Copy `pooling_sites_template.csv` to the ignored local `pooling_sites.csv` to register reviewed aggregate runs for central pooling. Specify validated observation dates and culture/AST QC eligibility per site. See `docs/pooling.md`; do not declare unresolved mappings validated.

Dependency setup uses the committed `.Rprofile` and `renv/activate.R` bootstrap (renv 1.1.5). Restore with `Rscript -e 'renv::restore(prompt = FALSE)'` from the repository root. The site runner checks dependencies and table schemas before analysis, and writes aggregate date availability to `provenance/preflight_source_availability.csv` in the completed run.
