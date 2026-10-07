# Culture practices by ASE hospitalization status

The additional analysis divides adult ICU hospitalizations into **ASE** and **Non-ASE** groups. It runs after the existing whole-cohort analyses and writes aggregate results under `output/runs/<run_id>/ase/`. Sites keep the same four configuration fields and run `Rscript code/00_run_pipeline.R` after restoring the updated lockfile.

## Membership and interpretation

Primary membership is retrospective: a hospitalization is ASE if **any** blood-culture episode during the full hospitalization meets presumed serious infection plus qualifying acute organ dysfunction. Otherwise an evaluable adult hospitalization is Non-ASE. All ICU stays, ICU-days, and cultures from that hospitalization inherit the same label, including observations before the event. This is not a time-varying classification and is not limited to ASE events that begin inside the ICU or the study window. Only ICU exposure and culture outcomes inside the shared study dates are analyzed.

Non-ASE means surveillance criteria were not met; it does not establish absence of infection. Missing source tables/columns or unverified component activity skip the additional stage with an explicit status. Children, missing ages, and invalid hospitalization boundaries are counted in exclusion QC rather than assigned Non-ASE. The whole-cohort analysis is unchanged.

ASE includes blood-culture ordering in its definition. Differences in blood-culture collection between these groups are partly built into group membership. Changes in culturing, treatment duration, or organ-function measurement can also change who qualifies. Results describe surveillance-defined populations and should not be interpreted as causal effects of sepsis on culturing or as risk-adjusted site comparisons.

## Definition and reference

Reference: [clifpy ASE implementation](https://github.com/Common-Longitudinal-ICU-data-Format/clifpy/blob/12e16ba9891ea6512fd464859565fdc3ed5ed8da/clifpy/utils/ase.py), pinned at commit `12e16ba9891ea6512fd464859565fdc3ed5ed8da`; [CDC surveillance toolkit](https://www.cdc.gov/sepsis/media/pdfs/sepsis-surveillance-toolkit-aug-2018-508.pdf). Attributed SQL and the Apache-2.0 license are stored in `utils/ase/`. R executes the queries using pinned DuckDB/DBI packages; no Python setup is required.

- Presumed infection requires a blood culture (any result) and a new qualifying antimicrobial course starting within ±2 calendar days, including new IV/IM treatment and normally at least four qualifying antimicrobial days. The reference's single-day-gap, drug-switching, oral/IV vancomycin distinction, and short-course death/transfer/hospice censoring logic are retained.
- Primary organ dysfunction excludes lactate, following the linked code's default. It includes new qualifying vasopressor infusion outside procedural care, new invasive ventilation, doubling creatinine excluding ESRD, bilirubin at least 2 mg/dL and doubled, or platelets below 100 × 10³/µL with at least a 50% fall from an eligible baseline.
- Community and hospital laboratory baselines follow the reference queries. Values above the reference outlier caps are excluded. Sites should verify CLIF lab units and clinical mapping before interpreting results. The numeric thresholds expect creatinine/bilirubin in mg/dL, platelets in 10³/µL, and lactate in mmol/L; UCMC’s creatinine, bilirubin, and platelet units were checked against these conventions.
- Lactate-inclusive flags remain in the private classification and the aggregate QC reports additional hospitalizations that would qualify. Separate lactate-inclusive trend analyses are not produced.
- A repeat-infection timeframe changes episode counts, but not whether a hospitalization has any qualifying episode. This analysis does not count ASE episodes and does not apply an episode filter.
- Calendar-day windows and month allocation use the project's UTC timestamp convention, consistently with existing analyses. Sites must ensure their timestamps are correctly standardized.

Documented adaptations from the linked code:

1. Blood-culture anchors require `method_category == culture`, with duplicate hospitalization/collection timestamps collapsed.
2. QAD hospitalization limits compare dates with dates, preserving admission-day treatment and up to two pre-admission calendar days of ED care.
3. Vasopressor and ventilation windows use ±2 **calendar days**, instead of the reference's timestamp-based ±48 hours.
4. Only administered, positive-dose qualifying antimicrobial records and the five CDC vasopressors are used; held/ordered drugs are excluded.
5. ESRD exclusion uses explicit `N18.6`; the broader upstream list includes `I27.2`, which is not an ESRD code.
6. The censoring query additionally recognizes an explicit `comfort_measures_only` discharge category. Unrecorded comfort transitions cannot be inferred.
7. Raw lab values may supply the reference SQL's numeric fallback rather than being discarded by the upstream loading wrapper.

## Required source data

The stage needs `patient`, `hospitalization`, `adt`, `microbiology_culture`, `medication_admin_intermittent`, `medication_admin_continuous`, `labs`, `respiratory_support`, and `hospital_diagnosis`. Full hospitalization data are needed, including baseline labs, antibiotic history, follow-up treatment, and discharge disposition outside the culture study dates. Merely supplying ICU-only extracts is insufficient.

A missing susceptibility table skips only the ASE-stratified susceptibility analysis. Missing ASE inputs never cause every patient to be labeled Non-ASE. Availability statuses, component activity, classification/exclusion QC, and the pinned definition are recorded in the aggregate export.

## Outputs and denominators

- Cohort counts and monthly ICU admissions/ICU-days, including uncultured stays. Carry-in stays contribute ICU-days but not new admissions.
- Monthly culture collection, positive-event rates, positivity, and exact specimen-category summaries. Categories are not selected independently by subgroup or collapsed into different Other bins.
- Admission-level proportion cultured and first-culture timing among new ICU admissions. Each admission's first culture is measured from its own ICU entry; admissions without a culture remain in the proportion denominator.
- Organism distribution, event/category detection counts, and monthly detection rates. Each organism category is counted once per culture event. Rates reflect detection events, including repeat cultures, not incident infections.
- Separate season-adjusted negative-binomial GAMs for culture and organism counts, using each group's ICU-day and ICU-admission offsets. Annualized endpoint contrasts, pointwise curves, residual diagnostics, and BH-adjusted screens use the existing shared modeling code. These are within-group time contrasts, not ASE-versus-Non-ASE effect estimates.
- If AST is available: subgroup organism–drug–specimen counts, tested fractions, QC, susceptible/non-susceptible detection models, and quasi-binomial non-susceptibility models using the existing testing/linkage rules.

Both groups share the full study calendar. Months without subgroup exposure have unavailable rates. Organism zero counts can enter models when the subgroup has exposure and the site still has observed culture-source activity. AST testing remains unknown when interpretable results are unavailable; it is never inferred susceptible.

Identifier-bearing hospitalization labels and criteria timestamps stay under ignored `data/intermediate/runs/<run_id>/ase/`. The existing export audit runs after this stage. Current central pooling does not automatically include these subgroup exports; subgroup pooling requires a separate extension.

## Validation

Run `Rscript tests/test_ase.R` for clinical boundary cases and subgroup denominator/count reconciliation. Before buddy release, reconcile samples of ASE and Non-ASE hospitalizations against source records, review lab units, antimicrobial/route/MAR mapping and observation coverage, and confirm exclusion counts. Synthetic fixtures validate implementation behavior, not clinical accuracy at a site.

## Subgroup demographics and characteristics

`ase_characteristics_table1_*.csv` provides a side-by-side ASE and Non-ASE Table 1. `ase_characteristics_long_*.csv` contains numeric counts, percentages, medians and quartiles for reuse across sites. Both are automatically included in the site export manifest. `ase_characteristics_qc_*.csv` reports hospitalization totals, distinct patients overall and patients represented in both groups.

The analysis unit is one classified adult hospitalization with ICU time overlapping the study window, including hospitalizations without cultures. Demographics are age at hospitalization admission and standardized patient sex, race and ethnicity categories. Characteristics include admission type, discharge disposition, full hospital length of stay, ICU days and number of ICU stays overlapping the study window, ICU culture counts and whether any ICU culture or positive culture was collected. Hospital length of stay and discharge disposition describe the completed hospitalization; they are not baseline severity measures. ICU and culture measures cover only study-window ICU time.

Categorical percentages use all hospitalizations in each subgroup as the denominator, including missing/unknown values shown as a separate row. Continuous summaries are median [25th, 75th percentile] among observed values. Every characteristic includes observed and missing counts; absent optional demographic columns are explicitly marked `source_column_unavailable`. Missing, unknown or empty source labels are treated as missing. Patients with multiple hospitalizations contribute once per hospitalization and may appear in both groups; subgroup patient counts should not be added to obtain a unique cohort count. No row-level identifiers or dates are exported. These are descriptive comparisons without hypothesis tests or causal interpretation.

Additional clinical descriptors report in-hospital death (`discharge_category = expired`), recorded IMV (`device_category = imv`) and recorded nonprocedural vasopressor use (administered positive-dose norepinephrine, dopamine, epinephrine, phenylephrine or vasopressin). Support observations must fall within admission-inclusive/discharge-exclusive hospitalization boundaries. These are whole-hospitalization recorded treatments, not baseline variables or ASE-specific organ dysfunction: no blood-culture proximity or new-initiation requirement is applied. A “no” indicates no qualifying recorded treatment, conditional on the available source tables; it does not prove complete chart capture. Missing discharge disposition leaves mortality unknown. These measures can overlap the criteria used to classify ASE and should be interpreted descriptively.
