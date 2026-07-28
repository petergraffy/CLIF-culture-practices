# ICU Culture Cohort Identification

## Cohort Definition

The initial analytic cohort includes all CLIF patients/hospitalizations with at least one microbiology culture collected during an ICU interval.

ICU intervals are identified from `adt` rows where `location_category == "icu"`. Culture timing is based on `microbiology_culture.collect_dttm`.

## Culture Export Definition

For qualifying hospitalizations, export every `microbiology_culture` row satisfying:

1. `method_category == "culture"`
2. non-missing `collect_dttm`
3. `collect_dttm` falls between ICU `in_dttm` and ICU `out_dttm`

No specimen/fluid category restrictions are applied. Blood, urine, respiratory, wound, CSF, stool/GI, and all other culture types are retained.

## Study Window

By default, the script includes all available culture collection dates. To restrict the run, set:

- `STUDY_START_DATE`, inclusive, formatted as `YYYY-MM-DD`
- `STUDY_END_DATE`, inclusive, formatted as `YYYY-MM-DD`

For example:

```sh
STUDY_START_DATE=2018-01-01 STUDY_END_DATE=2025-12-31 Rscript code/01_identify_icu_culture_cohort.R
```

## Outputs

`code/01_identify_icu_culture_cohort.R` writes timestamped CSVs under `output/cohort/`:

- ICU culture rows
- collapsed culture-event summaries
- cohort hospitalizations
- ICU intervals with culture
- cohort summary
- fluid/specimen summary

## Notes

- ICU `out_dttm` is filled with hospital discharge time when ADT has a missing ICU out time.
- The row-level export preserves organism-level rows. A separate event-level summary collapses rows sharing patient, hospitalization, ICU interval, order/collection time, fluid, and method.
- Positive culture is initially defined as a non-missing organism group/category other than `no_growth` or `no growth`.
