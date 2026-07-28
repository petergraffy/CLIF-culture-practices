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

## Initial Analysis Plan

- Identify the ICU culture cohort as patients/hospitalizations with at least one culture collected during an ICU interval, then retain every culture row collected during ICU time across all specimen types.
- Define the denominator as patients with at least one culture collected.
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
- `output/`: generated tables, figures, and intermediate non-sensitive outputs

## First Analysis Step

Configure local CLIF table paths with `config/config.json`, then run:

```r
source("code/01_identify_icu_culture_cohort.R")
```

This writes timestamped ICU culture cohort exports under `output/cohort/`.

Optional date-window environment variables restrict by culture collection time:

```sh
STUDY_START_DATE=2018-01-01 STUDY_END_DATE=2025-12-31 Rscript code/01_identify_icu_culture_cohort.R
```

## Time-Series Plots

After cohort identification, run:

```sh
ICU_CULTURE_EVENTS_PATH=output/cohort/icu_culture_events_UCMC_YYYYMMDD_HHMMSS.csv Rscript code/02_plot_culture_time_series.R
```

Optional plot controls:

```sh
TOP_N_CULTURE_TYPES=8 PLOT_END_DATE=2024-12-31 Rscript code/02_plot_culture_time_series.R
```

This writes monthly summaries and PNG figures under `output/time_series/`.

## Positive Organism Plots

After cohort identification, run:

```sh
ICU_CULTURE_ROWS_PATH=output/cohort/icu_culture_rows_UCMC_YYYYMMDD_HHMMSS.csv Rscript code/04_plot_positive_organisms.R
```

Optional controls:

```sh
TOP_N_CULTURE_TYPES=8 TOP_N_ORGANISMS_PER_TYPE=10 PLOT_END_DATE=2024-12-31 Rscript code/04_plot_positive_organisms.R
```

This writes positive organism summaries and PNG figures under `output/organisms/`.

## Data Governance

Do not commit PHI, row-level CLIF extracts, credentials, or institution-specific restricted files. Use local paths, environment variables, or ignored private directories for sensitive inputs.
