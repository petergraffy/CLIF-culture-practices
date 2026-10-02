# Culture practices: analysis definitions

## Population and counting units

The practices denominator is every valid ICU stay represented in ADT, including stays without culture. Culture results and organism distributions are summarized among observed ICU cultures. No adult-only restriction is imposed by the current project.

ICU stays merge overlapping and contiguous ICU ADT rows within patient/hospitalization. Missing ICU exit time is replaced by hospital discharge; the frequency of this imputation is exported by QC. ICU intervals are half-open `[entry, exit)`: a specimen at exit does not belong to the departing stay. Study dates are inclusive calendar dates, implemented using an exclusive midnight after the end date.

One culture event is patient/hospitalization/merged ICU stay/order time/collection time/source fluid name/source method name. Standardized categories are descriptive attributes, not event identifiers. Latest available `result_dttm` is retained within event/isolate; when `organism_id` is absent, organism labels provide the fallback isolate key. Named organisms remain separate within polymicrobial events. Conflicting AST results cannot be resolved by chronology because the CLIF susceptibility schema has no result timestamp.

The event definition cannot distinguish simultaneous separate specimens with identical fields when the source has no specimen identifier. Sites should validate this counting rule against their source laboratory system.

## Results and yield

Rows are classified into positive, negative/no growth, mixed/contaminated, or indeterminate. Explicit negative text and missing/pending/cancelled results are not organisms. Mixed or normal flora and explicit contamination are separated from named isolates. Event precedence is positive, mixed/contaminated, indeterminate, negative/no growth. A named isolate plus mixed flora is positive; no growth plus an unresolved report is indeterminate.

Positivity uses `positive / (positive + negative/no growth)`. Exports retain counts of all four categories and interpretable events. Mixed/contaminated and indeterminate events are excluded from the yield denominator, and should be reported alongside yield.

## Denominators

Calendar-month ICU days are exact duration overlaps with that month, including stays beginning before the study window. Admission counts use actual ICU entry dates; carry-in stays contribute days without becoming new admissions. Calendar-month events per 100 admissions describe utilization relative to incoming volume, not the fraction of admissions cultured.

`monthly_admission_cohort_culture_proportions` assigns culture availability to the ICU admission month and includes uncultured stays. It reports stays with follow-up truncated at study end. Timing analyses include stays entering within the study window, with observation censored at study end. Counts of active stays cultured in a collection month are labeled as active-stay counts, without interpreting them as admission-cohort proportions.

## Specimen continuity

`10_quality_checks.R` reports category observation ranges and raw fluid-name/category changes. Categories observed for at least 12 months but absent for the final three or more months are flagged for review. This is a screening flag, not proof of an ETL error. Rare or newly introduced categories require manual review too.

At UCMC, the existing exploratory outputs show lower respiratory and pleural categories disappearing after April 2023. No replacement mapping has been inferred: source fluid names can be test order names and need not identify the actual specimen. Confirm source ETL before interpreting that period.

Verified corrections can be entered in `config/specimen_category_overrides.csv`: exact fluid name, original category, inclusive dates, replacement category, documented reason. Original categories are preserved privately; overlapping rules fail. If the mapping cannot be established, restrict `study_start_date`/`study_end_date` to a validated period and rerun. Do not treat a missing category as evidence of zero clinical sampling.

## Temporal analysis

Detection counts use negative-binomial generalized additive models (`mgcv::gam`, REML), a penalized cubic smooth of elapsed calendar time (4–8 basis dimensions, based on available months), annual sine/cosine seasonal terms, and log ICU-day or admission offsets. ICU-day models are preferred for calendar-time detection burden; admission models are companion utilization screens.

This handles extra-Poisson variance and permits reversals or plateaus. The main descriptive contrast compares season-standardized fitted endpoints and annualizes their ratio over the observed period. It is a net change, not a constant annual slope or a claim that every intermediate month increased. Smooth-term temporal p-values are exported separately from endpoint-contrast p-values. Both endpoint intervals and monthly fitted curves are exported. Fitted curves and their pointwise 95% confidence intervals hold the annual sine/cosine terms at zero (their cycle means), displaying the season-adjusted long-term component. Observed monthly bars or points retain seasonal variation. This is standardization on the model link scale, not an arithmetic average of seasonal response-scale rates. Seasonal terms remain in model fitting and residual diagnostics. Months with no observed culture events are excluded from organism temporal inference, since source coverage cannot be distinguished from absence of clinical sampling. At least 24 usable months and 25 detections are required; otherwise models are explicitly not estimated.

BH FDR applies to endpoint contrasts separately within the organism or targeted-text screen and denominator. Direction labels require FDR < 0.05 and no residual-dependence flag. Pearson residual correlation between adjacent calendar months above 0.3 in magnitude is flagged (at least 12 adjacent pairs are required), and those results remain exploratory. The GAM does not itself model serial dependence; persistent correlation requires a time-series extension or block-bootstrap inference before confirmatory claims. Endpoint confidence intervals do not correct for serial dependence. Rank plots select organisms by ICU-day endpoint contrasts, show GAM curves with intervals, and display estimated net change even when not significant; consult model status, FDR, and residual diagnostics.

Reference: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/negbin.html

## Susceptibility analysis

Sites with both microbiology tables run `09_susceptibility_trends.R`. Missing susceptibility tables produce an availability record and a clean skip. Raw organism resistance names remain a separate text-only screen in script 08; they do not substitute for susceptibility tests.

Pinned mCIDE files in `config/mcide` were obtained from the CLIF main branch on October 2, 2026:
https://github.com/Common-Longitudinal-ICU-data-Format/CLIF/tree/main/mCIDE/microbiology_susceptibility

The authoritative fields are `organism_id`, `antimicrobial_category`, and `susceptibility_category`. `antimicrobial_name`, `sensitivity_name`, and `susceptibility_name` may be retained in the site's source table but are not used to reinterpret a canonical category. No MIC breakpoints are inferred. Missing standardized interpretation is unavailable even when raw text says resistant. Categories are susceptible, non_susceptible, indeterminate, or unavailable; unexpected categories and drugs are counted in QC.

Only positive, identifiable ICU isolates are linked by `organism_id`. An ID mapping to multiple event/category combinations fails rather than allowing cross-linkage. Repeated records for one isolate/drug collapse to one observation. Conflicting susceptible and non-susceptible reports become indeterminate and are audited. A known S or NS report with duplicate unavailable/indeterminate entries remains interpretable unless S and NS conflict.

Analyses are organism–antimicrobial specific, overall and stratified by specimen. An organism has no universal susceptible status: it may be susceptible to one drug and non-susceptible to another. The analysis exports:

- Susceptible and non-susceptible isolate detections per 100 ICU days, with NB GAM net-change contrasts.
- Non-susceptible fraction among interpretable tests, with a quasi-binomial GAM using the same flexible time and seasonal terms; its endpoint effect is an annualized odds ratio, not a rate ratio.
- Testing coverage: interpretable tests / all ID-linked positive isolates of that organism, plus unavailable, indeterminate, conflicting, and absent tests.
- Fitted curves, confidence intervals, model status, endpoint FDR, and residual flags.

Rate models include only months with at least `AST_MIN_TESTING_FRACTION` interpretable coverage (default 0.5). Months with zero interpretable tests have missing displayed detection rates and missing tested fractions, not apparent absence of non-susceptibility. Fraction models use months with interpretable tests; at least 24 usable months, 30 interpretable tests, and 10 observations in each outcome category are required. FDR families are outcome × specimen stratum within site.

These describe clinically tested isolates, not population resistance prevalence. Testing selection, repeat cultures, changing panels, breakpoint revisions, and specimen/case mix can affect trends. Pooling requires a common observation window, stable mappings, comparable panels, and site-level QC. Sites should consider first-isolate and stable-panel sensitivities before inferential pooling.

## Reproducibility and exports

Run `Rscript code/00_run_pipeline.R` from the repository root. Every run has its own `output/runs/<run_id>/` aggregate tree and ignored `data/intermediate/runs/<run_id>/` private tree. All downstream scripts use that run's intermediates. Standalone scripts remain available, but reject legacy intermediates without the shared event identifier.

The run manifest records config/code/mapping hashes, R/package versions, date window, and completion status. Exact source, mapping, lockfile, and config snapshots are retained with the private run provenance, including uncommitted code. Input metadata records filenames, sizes, and modification times; source table content is not hashed automatically. Exact local paths/environment overrides are private. Export audit rejects event/organism/patient/stay identifiers and clinical timestamps, including the new event IDs. It is a column-name audit, not small-cell suppression or a replacement for site release review.

Use `renv::restore()` to install pinned dependencies. Tests: `Rscript tests/test_core.R` and `Rscript tests/run_integration.R`. UCMC cannot validate real susceptibility ETL; synthetic tests exercise the standardized schema and optional-table behavior. A site with real AST is still needed for clinical validation.

## Cross-site pooling

The central workflow in `code/11_pool_site_trends.R` refits site models over common usable months, pools log endpoint contrasts, and fits secondary joint season-adjusted trajectories. See [pooling definitions](pooling.md) for eligibility, equal-site weighting, heterogeneity, inference limits and registry configuration. Site model exports now include the log endpoint effect/SE, NB theta, long-term EDF and smoothing parameter.
