# CLIF Culture Practices

This repository supports a CLIF-wide study of microbiology culture practices and culture results over time.

## Study Objective

Describe variation in microbiology culture acquisition across CLIF sites and characterize positive culture rates and organism distributions among patients with any culture collected.

## Cohort

Use all ICU stays as the acquisition denominator, including stays without cultures. Describe results among cultures collected during those stays. Retain all specimen types and summarize yield by specimen, site, and time. See [analysis definitions](docs/analysis_definitions.md).

## Core Questions

1. How often are cultures collected across sites, care settings, specimen types, and calendar time?
2. What proportion of cultures are positive overall and within specimen type, site, care setting, and time strata?
3. Which organisms and organism groups predominate overall and over time?
4. How much variation in culture positivity and organism mix reflects specimen type, site practice, patient case mix, and secular trends?

## Analysis Plan

- Identify all ICU admissions from CLIF ADT rows, merge overlapping or back-to-back ICU intervals within hospitalization, and identify ICU culture events collected during ICU time.
- Use all ICU admissions and ICU days as denominators for practice-rate analyses.
- Summarize cultures per patient, cultures per encounter, specimen type mix, and timing relative to admission or ICU time where available.
- Classify culture result status as positive, negative/no growth, contaminated/mixed flora when distinguishable, and indeterminate/missing.
- Map organisms to clinically meaningful groups, preserving organism-level detail for common isolates.
- Estimate positive culture rates by specimen type, site, calendar month or quarter, and care setting.
- Describe temporal trends in organism groups and common organisms.
- Evaluate between-site variation using stratified summaries first, then multivariable or hierarchical models if needed.

## Repository Layout

- `code/`: analysis scripts
- `config/`: study configuration and organism/specimen mapping files
- `data/`: non-sensitive data documentation and derived public metadata only
- `docs/`: protocol notes, data dictionaries, and manuscript materials
- `output/`: generated aggregate tables, figures, and site export manifests only

## Site Export Rule

Everything intended for pooling or cross-site comparison is written under `output/`.

Do not place row-level CLIF extracts in `output/`. Scripts that need local row-level intermediates write them under ignored `data/intermediate/` by default. Those files can include patient, hospitalization, ICU interval, timestamp, and microbiology row identifiers and should not be shared.

## First Analysis Step

Configure local CLIF table paths with `config/config.json`, then run:

```sh
Rscript code/01_identify_icu_culture_cohort.R
```

This writes aggregate cohort summaries under `output/cohort/` and private intermediates under `data/intermediate/cohort/`. For an isolated, fully documented run, use the recommended pipeline command below.

By default, this also writes private row-level intermediates under `data/intermediate/cohort/` for scripts that still use a local cohort extract. Disable those private intermediates with:

```sh
WRITE_ROW_LEVEL_INTERMEDIATES=false Rscript code/01_identify_icu_culture_cohort.R
```

The study team maintains the shared inclusive study window in `utils/study_settings.R`, currently January 1, 2018–December 31, 2024. All site analyses and plots use that window. Sites configure only their label, repository/table paths, and file type.

## Time-Series Plots

After cohort identification, run:

```sh
Rscript code/02_plot_culture_time_series.R
```

Optional plot controls:

```sh
TOP_N_CULTURE_TYPES=8 Rscript code/02_plot_culture_time_series.R
```

The script reads the latest private event file from `<repo>/data/intermediate/cohort/`. It writes monthly aggregate summaries and PNG figures under `<repo>/output/time_series/`.

## Positive Organism Plots

After cohort identification, run:

```sh
Rscript code/04_plot_positive_organisms.R
```

Optional controls:

```sh
TOP_N_CULTURE_TYPES=8 TOP_N_ORGANISMS_PER_TYPE=10 Rscript code/04_plot_positive_organisms.R
```

The script reads the latest private culture row file from `<repo>/data/intermediate/cohort/`. It writes aggregate positive organism summaries and PNG figures under `<repo>/output/organisms/`.

## Organism Trend Screen

After cohort identification, run:

```sh
Rscript code/08_organism_trends.R
```

This screens organism detection with negative-binomial GAMs, smooth calendar time, annual seasonality, offsets, endpoint contrasts, FDR correction, and residual diagnostics. Rank plots describe net fitted change, not a constant annual slope. Text-reported resistance labels remain explicitly separate from susceptibility-derived analyses.

## Optional Susceptibility Analysis

`code/09_susceptibility_trends.R` runs at sites with both microbiology tables. It uses the CLIF mCIDE `organism_id`, `antimicrobial_category`, and `susceptibility_category` fields. For each organism–antimicrobial pair it reports susceptible/non-susceptible detections per ICU day, the non-susceptible fraction among interpretable tests, testing coverage, and flexible temporal models. Missing tables produce a clean skip; missing tests never become susceptible. See [definitions and limitations](docs/analysis_definitions.md).

## Recommended Multi-Site Run

For each site, create `config/config.json` with `site_name`, `repo`, `tables_path`, `file_type`, and the study window. Install pinned dependencies with `renv::restore()`, then run from the repository root:

```sh
Rscript code/00_run_pipeline.R
```

The runner executes cohort construction, continuity QC, descriptive plots, rates, ICU-day/timing analyses, organism trends, optional susceptibility trends, and the export audit. It writes aggregate outputs to `output/runs/<run_id>/` and private intermediates to `data/intermediate/runs/<run_id>/`. Config/code/mapping hashes and package versions are recorded in run manifests. Share a single completed run after reviewing continuity flags and site release rules.

For manual run order and overrides, see [code/README.md](code/README.md). Run synthetic unit checks with `Rscript tests/test_core.R`.

## Data Governance

Do not commit PHI, row-level CLIF extracts, credentials, or institution-specific restricted files. Use local paths, environment variables, or ignored private directories for sensitive inputs. Share only aggregate files from `output/` after the site export privacy audit passes.

## Cross-site pooling

After sites generate completed aggregate runs, use `Rscript code/11_pool_site_trends.R config/pooling_sites.csv`. This refits site models over common validated months, pools endpoint changes with random-effects meta-analysis, and fits secondary joint season-adjusted curves. See [pooling definitions and registry setup](docs/pooling.md). No cross-site estimates are produced by the single-site runner.

Before multisite testing, follow the [buddy-testing guide](docs/buddy_testing.md), including source linkage/coverage reconciliation at a site with real AST.
