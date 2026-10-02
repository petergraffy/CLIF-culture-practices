# Cross-site organism and susceptibility trends

This is a central, aggregate-only workflow, run separately from each site's pipeline. It produces a primary random-effects pooled net change and secondary joint season-adjusted trajectories. Current UCMC outputs alone cannot produce a cross-site estimate. Synthetic examples validate execution, not clinical generalizability.

## Prepare site exports

Run `Rscript code/00_run_pipeline.R` at each site using matching analysis code and pinned mCIDE files. Use a completed, clean-audit `output/runs/<run_id>` directory. Do not mix files from different runs.

Script 08 now exports `monthly_all_organism_detection_counts_*.csv` with standardized organism-category keys, site, month and counts. This is independent of local top-N ranks. Pooling rejects older runs missing this export: a category absent from a top-N list cannot be interpreted as zero.

Required aggregate inputs are the all-organism counts, monthly ICU denominators/source coverage, public run manifest and clean export audit. Sites with AST also supply script 09's monthly organism–antimicrobial–specimen aggregates. Private culture/AST rows are never read by the central script.

Copy `config/pooling_sites_template.csv` to the ignored local file `config/pooling_sites.csv`. Enter exactly one completed run per site:

| Column | Meaning |
| --- | --- |
| `site_name` | Must match the run manifest; unique across rows |
| `run_dir` | Absolute path to the received aggregate run directory |
| `validated_start_date`, `validated_end_date` | Inclusive dates within that run's study window with validated source coverage and mappings |
| `culture_qc_pass` | TRUE only after site ETL, organism definitions and relevant specimen mappings are reviewed |
| `ast_qc_pass` | TRUE only for validated AST linkage, categories, testing panels and interpretation comparability; FALSE for sites without usable AST |

These flags are scientific eligibility metadata, not an automatic validation result. UCMC's unresolved specimen discontinuity still requires a validated window or mapping correction. A clean privacy audit does not establish clinical QC or authorize data release.

Run from the repository root:

```sh
Rscript code/11_pool_site_trends.R config/pooling_sites.csv
```

An optional second argument chooses a new/empty output directory. Default output is an isolated `output/pooling/<timestamp>_<pid>` directory. `renv::restore()` installs the pinned dependencies, including `metafor`. `AST_MIN_TESTING_FRACTION` defaults to 0.5 and `AST_MIN_LINKAGE_FRACTION` to 0.9, as in the site analysis. Central screening recomputes eligibility from exported counts; the manifest records both central thresholds.

## Shared observation window and eligibility

Only complete calendar months within each site's validated interval are retained. For each organism/outcome/denominator/specimen screen, models are refit over the intersection of usable months across participating sites. The endpoints and month set are identical across sites; gaps are allowed and reported by the common-month count. Original full-window effect estimates are not combined.

Months need positive denominators and observed culture-source activity. AST months containing organisms need adequate linkage completeness; detection months additionally need interpretable testing coverage using all observed isolates, while tested-fraction months need interpretable tests. Monthly specimen-specific source activity is preserved rather than replaced by overall culture counts. Months without source activity are unavailable, not presumed zero. An organism category absent from a complete all-organism export can be filled with zero using validated denominators. Absent AST pairs are never filled as susceptible or as zero resistance. For an observed pair, genuine zero-organism months remain eligible for rate models only with the site's explicit `culture_coverage_validated` declaration, specimen-source activity, and complete organism identification; they never contribute to tested fractions. Older AST exports lacking observation/linkage fields are rejected and require a fresh site run.

All sites with data for an AST pair enter its common-month selection, but pairs missing at a site omit that site entirely. The registry audit records AST availability. Per-screen site effects expose sparse-data exclusions and the actual number of common months. After common-month alignment, each site must meet the existing minimum counts/months/outcome-variation requirements. These exclusions do not depend on p-values or trend direction. The estimand applies to eligible sites, not necessarily every CLIF site.

## Primary analysis: random-effects meta-analysis

Each site's GAM exports the annualized endpoint log ratio and its standard error directly, plus theta, long-term effective degrees of freedom and smoothing parameter. Ratios are season-adjusted endpoint contrasts, not constant slopes.

For each screen, `metafor::rma.uni` pools the log ratios using REML between-site variance estimation and modified Knapp–Hartung inference (the standard-error adjustment is floored at 1 to prevent narrower-than-unadjusted intervals). Outputs include the pooled ratio and 95% CI, p-value, tau-squared, I-squared, heterogeneity-test p-value, and a prediction interval when at least three sites are eligible. Forest plots show every estimable site, including nonsignificant results.

At least two sites are needed for a pooled estimate; two-site results are explicitly exploratory. Prediction intervals and heterogeneity estimates are uncertain with few sites. Susceptibility detection-rate ratios and tested-fraction odds ratios remain separate. BH FDR is applied to pooled endpoint p-values within screen × outcome × denominator × specimen stratum.

Site serial-dependence flags propagate to pooled inference. If any contributing site is flagged or lacks an adequate diagnostic, the pooled result is labeled exploratory and is not assigned a significant increasing/decreasing label. The pipeline does not repair serial dependence; block-bootstrap or a validated correlation model remains necessary before confirmatory inference for those screens.

## Secondary analysis: joint trajectories

At least three eligible sites are required. The joint NB GAM uses:

```r
y ~ site + s(time_years, bs = "cr", k = k) +
  s(time_years, site, bs = "sz", k = k) +
  site:season_sin + site:season_cos + offset(log_denominator)
```

Site baseline rates are fixed effects. Penalized sum-to-zero site smooth deviations allow different trajectories around the shared long-term smooth. Seasonal coefficients are site-specific. Basis size is based on common months (4–8), rather than total site-month rows; REML estimates smoothing. Count models have one joint NB dispersion parameter, while the primary site models estimate dispersion separately. Fraction models use a quasi-binomial response `cbind(non_susceptible, tested - non_susceptible)` and omit the offset.

Predictions set seasonal terms and count offsets to zero, then transform back to per-100 rates or fractions. The pooled curve is the arithmetic mean of site predictions with fixed **equal site weights**. It describes an average of participating sites; it is not an ICU-day-weighted network burden curve or the curve for a new site. A large site's higher precision still influences joint coefficient estimation.

Pointwise 95% CIs use the full coefficient covariance and a delta-method gradient for the response-scale mean, then a log/logit transformation to maintain valid bounds. They include smoothing-parameter uncertainty but are not simultaneous bands or between-site prediction intervals. Site residual diagnostics are exported. Curves remain secondary/descriptive, particularly with residual dependence, sparse sites or discordant trajectories.

## Outputs and provenance

- `pooled_effects.csv`: refitted site effects, counts/status, common endpoints/month count and model diagnostics.
- `pooled_meta.csv`: random-effects results, heterogeneity, prediction intervals, pooled FDR and inference status.
- `pooled_joint.csv`: joint fit status, dispersion, long-term EDF and residual flags.
- `pooled_curves.csv`: site and equal-site joint predictions with pointwise CIs.
- `pooled_diagnostics.csv`: site residual diagnostics within each joint screen.
- `figure_index.csv`: screen keys mapped to forest and joint figure filenames.
- `site_registry_audit.csv`, `input_file_hashes.csv`, `pooling_manifest.json` and source snapshots: run identifiers, eligibility, aggregate input hashes and exact central code/dependency provenance.

Empty curve/diagnostic outputs are omitted when no joint model is estimable. The central workflow rejects identifier/timestamp columns, duplicate site/month rows, mismatched analysis/mCIDE versions, incomplete runs, invalid intervals and ambiguous exports. Outputs remain aggregate; small-cell release review still belongs to each institution.

Validation: `Rscript tests/test_pooling.R` exercises four synthetic sites, rate and susceptibility outcomes, calendar alignment, equal-site averaging, eligibility, input checks, optional AST, figures and the central CLI. `POOLING_TEST_OUTPUT` may specify a fresh directory to retain synthetic fixtures and figures after the R process exits.

References: [metafor random-effects models](https://wviechtb.github.io/metafor/reference/rma.uni.html) and [mgcv factor smooths](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/factor.smooth.html).
