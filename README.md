# CLIF Culture Practices

This repository supports a CLIF-wide study of microbiology culture practices and culture results over time.

## Study Objective

Describe variation in microbiology culture acquisition across CLIF sites and characterize positive culture rates and organism distributions among patients with any culture collected.

## Cohort

Include all patients in CLIF with at least one microbiology culture collected during the study window. All culture specimen types should be retained initially, with specimen-type-specific summaries used to describe practice variation and interpret organism yield.

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

```r
source("code/01_identify_icu_culture_cohort.R")
```

This writes timestamped ICU culture cohort exports under `output/cohort/`.

By default, this also writes private row-level intermediates under `data/intermediate/cohort/` for scripts that still use a local cohort extract. Disable those private intermediates with:

```sh
WRITE_ROW_LEVEL_INTERMEDIATES=false Rscript code/01_identify_icu_culture_cohort.R
```

Optional date-window environment variables restrict by culture collection time:

```sh
STUDY_START_DATE=2018-01-01 STUDY_END_DATE=2025-12-31 Rscript code/01_identify_icu_culture_cohort.R
```

## Time-Series Plots

After cohort identification, run:

```sh
ICU_CULTURE_EVENTS_PATH=data/intermediate/cohort/icu_culture_events_UCMC_YYYYMMDD_HHMMSS.csv Rscript code/02_plot_culture_time_series.R
```

Optional plot controls:

```sh
TOP_N_CULTURE_TYPES=8 PLOT_END_DATE=2024-12-31 Rscript code/02_plot_culture_time_series.R
```

If `ICU_CULTURE_EVENTS_PATH` is not set, the script reads the latest private event file from `data/intermediate/cohort/`. It writes monthly aggregate summaries and PNG figures under `output/time_series/`.

## Positive Organism Plots

After cohort identification, run:

```sh
ICU_CULTURE_ROWS_PATH=data/intermediate/cohort/icu_culture_rows_UCMC_YYYYMMDD_HHMMSS.csv Rscript code/04_plot_positive_organisms.R
```

Optional controls:

```sh
TOP_N_CULTURE_TYPES=8 TOP_N_ORGANISMS_PER_TYPE=10 PLOT_END_DATE=2024-12-31 Rscript code/04_plot_positive_organisms.R
```

If `ICU_CULTURE_ROWS_PATH` is not set, the script reads the latest private culture row file from `data/intermediate/cohort/`. It writes aggregate positive organism summaries and PNG figures under `output/organisms/`.

## Organism Trend Screen

After cohort identification, run:

```sh
ICU_CULTURE_ROWS_PATH=data/intermediate/cohort/icu_culture_rows_UCMC_YYYYMMDD_HHMMSS.csv Rscript code/08_organism_trends.R
```

This screens top organisms and targeted organism/resistance text labels for monthly detection-rate trends per 100 ICU admissions and per 100 ICU days. Trend plots are color-coded by organism taxonomy. True MRSA, VRE, ESBL, and CRE phenotypes require susceptibility or resistance fields; this script only identifies those labels when resistance terms appear in organism text.

## Recommended Multi-Site Run

For each site, set `CLIF_SITE_NAME`, `CLIF_TABLES_PATH`, and the study window, then run the aggregate-producing scripts:

```sh
STUDY_START_DATE=2018-01-01 STUDY_END_DATE=2024-12-31 Rscript code/01_identify_icu_culture_cohort.R
STUDY_START_DATE=2018-01-01 STUDY_END_DATE=2024-12-31 Rscript code/05_culture_rates_per_icu_admission.R
STUDY_START_DATE=2018-01-01 STUDY_END_DATE=2024-12-31 Rscript code/06_icu_day_denominators_and_timing.R
PLOT_START_DATE=2018-01-01 PLOT_END_DATE=2024-12-31 Rscript code/04_plot_positive_organisms.R
STUDY_START_DATE=2018-01-01 STUDY_END_DATE=2024-12-31 Rscript code/08_organism_trends.R
Rscript code/07_prepare_site_exports.R
```

Before pooling, confirm that `code/07_prepare_site_exports.R` completes without finding disallowed identifier or exact timestamp columns in `output/`.

## Data Governance

Do not commit PHI, row-level CLIF extracts, credentials, or institution-specific restricted files. Use local paths, environment variables, or ignored private directories for sensitive inputs. Share only aggregate files from `output/` after the site export privacy audit passes.
