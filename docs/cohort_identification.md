# ICU culture cohort identification

All valid ICU stays form the practices denominator. `code/01_identify_icu_culture_cohort.R` exports aggregate summaries of the cultured subset and private ICU culture/event intermediates.

ICU ADT intervals are merged when overlapping or contiguous within patient/hospitalization. Culture collection belongs to half-open `[ICU entry, ICU exit)` time. All specimen categories are included. Study start/end dates are inclusive calendar dates.

Culture events share patient/hospitalization/merged ICU stay/order time/collection time/source fluid name/source method name. Standardized categories are descriptive attributes. Latest available result per event/isolate is retained, and polymicrobial isolates remain distinct.

`output/cohort/` contains only aggregate cohort and fluid summaries for standalone runs. Private rows, events, cultured hospitalizations, and ICU stays are under ignored `data/intermediate/cohort/`. Automated runs place these within corresponding run-specific folders. Do not share private intermediates.

Use `code/00_run_pipeline.R` for normal site execution. See [analysis_definitions.md](analysis_definitions.md) for result precedence, positivity denominators, carry-in handling, temporal models, susceptibility analyses, and mapping QC.
