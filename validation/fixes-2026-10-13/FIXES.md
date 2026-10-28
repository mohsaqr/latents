# Audit repairs — 2026-10-13

All 27 confirmed findings in the [original audit](../audit-2026-10-13/REVIEW.md) have repairs and regression checks. These changes preserve the user's existing uncommitted work. Version 0.9.23 is unchanged; nothing was committed, pushed, released or submitted.

The original audit records the pre-repair state. Its characterization tests intentionally reproduce defects and must not be treated as correctness specifications for repaired code. New desired-behavior regressions live in the package test directories; final verification evidence lives here.

## Finding closure

| Finding | Resulting behavior | Source and regression evidence |
|---|---|---|
| SV01 | Bootstrap validates original indicators, row/group/time alignment and normalized sampling weights before resampling; bad inputs preserve RNG state. A likelihood check also protects older fits. | `R/inference.R`; service unit suite |
| SV02 | Gaussian residual moments and categorical observed/expected tables use fitted sampling weights. Unsupported weighted test statistics, degrees of freedom and p-values are `NA`. | `R/bivariate-residuals.R`, `R/get-data.R`; service unit/equivalence suites |
| SV03 | Rubin pooling handles zero within/between variance: exact constants have zero relative variance increase/FMI and infinite degrees of freedom; zero within variance with positive between variance has the proper infinite/one limits. | `R/pool-imputations.R`; service unit/equivalence suites |
| SV04 | Descriptives reject matrix/array-valued indicator, ID and grouping columns with a classed data error. | `R/descriptives.R`; service unit suite |
| SV05 | R3STEP scales predictor columns for fitting/rank checks and transforms coefficients/covariance back to original units; ordinary and robust results agree across unit changes. | `R/three-step.R`; service equivalence suite |
| RA01 | Classification follows comparison filtering/order; empty comparisons give a typed empty table. | `R/compare-models.R`; recent API suite |
| RA02 | General LTA grid comparisons use appropriate criteria and posterior summaries, including occasion-varying transitions. | `R/compare-models.R`, `R/enumeration.R`; recent API suite and weighted/categorical/count probes |
| RA03 | Homogeneous/general transition and membership-covariate fits compare directly or through their grids. | `R/compare-models.R`; recent API suite |
| RA04 | Smallest-profile curves distinguish covariance structures. Noise mass stays in Gaussian-share denominators; older membership fits reconstruct missing effective counts. | `R/selection.R`; recent API suite and structure HTML |
| RA05 | Regression-comparison `"all"` requests only supported tables. | `R/get-data.R`, `R/compare-models.R`; recent API suite |
| RA06 | Family/transition enumeration `"all"` returns candidates and best; all-failed grids return available candidates. | `R/additive-compare.R`, `R/lta-compare.R`; recent API suite |
| RA07 | Same-data checks include raw continuous values, centering offsets/scheme, categorical labels/codes, extra indicators, groups/order, occasions and weights. Raw shifts/category relabeling are refused. | `R/compare-models.R`; recent API suite |
| RA08 | The calibrated statistic is accurately described as a maximum perpendicular chord-distance knee, a simplified Kneedle-inspired rule without smoothing/sensitivity thresholds. | `R/selection.R`, `R/get-data.R`, enumeration documents; original calibration test and renders |
| RA09 | Noise entropy uses the complete posterior matrix, even with one Gaussian profile plus noise. | `R/compare-models.R`, `R/enumeration.R`; recent API suite |
| RA10 | Covariate enumeration tolerates absent replication counts, reporting `NA`, and resolves available covariance metadata. | `R/enumeration.R`; recent API suite |
| MC01 | Covariate centering propagates through fitting, stored indicators/offsets, reconstruction, printing and Wald inference. Fit-table centering is appended without displacing existing columns. | `R/multilpa.R`, `R/covariates.R`, `R/covariate-inference.R`, `R/methods.R`; core suite |
| MC02 | Active covariance floors preserve volume, shape and orientation constraints; previous-covariance fallbacks must satisfy the floor and family. Active prior floors are repaired too. | `R/covariance-structure.R`; core suite checks all 14 structures, objective ascent, malformed starts and seven public fits |
| MC03 | Growth likelihood, variance updates and analytic scores use stable conditional residual/random-effect calculations when Woodbury subtraction loses precision; ordinary paths retain established arithmetic. | Growth engine/inference; core suite, independent scalar Gaussian/score references and existing growth equivalence suites |
| MC04 | Prior help states intentional EII scale, EEI mean/shrinkage, VEV degrees-of-freedom and mean-prior denominator conventions; established mclust behavior is retained. | `R/prior.R`, generated help; prior suites and active-floor tests |
| PD01 | Named noise diagnostics enforce generic restrictions: default `"all"` gives supported sizes; explicit entropy/posterior/AvePP requests are refused. | `R/plot-named.R`; plot/prediction suite |
| PD02 | Weighted tree shares and flows use weights at their fitted independent-unit level. | `R/plot-gg.R`; plot suite and weighted tree/table HTML |
| PD03 | Released enumeration/tree abbreviations and explicit `NULL` are normalized before criterion dispatch. | `R/plot.R`; plot suite covers `e`, `enum`, `t`, `tr`, `NULL` |
| PD04 | Shared help restores every released alias and all 14 plot-method usages/defaults; prediction/simulation/result-method aliases are preserved too. | Roxygen method annotations; documentation gate checks all 199 released aliases and pkgdown references |
| PD05 | Aggregated-binomial prediction consistently returns expected successes. New trials come from response columns or named-only `trials`; `trials = 1` gives probabilities. Missing/inconsistent exposures raise classed errors. | `R/mixture-regression-predict.R`, prediction help; plot/prediction suite |
| PD06 | Article states bootstrap p-value `(b + 1) / (B + 1)` and explains that invalid replications withhold the test. | Enumeration article render and bootstrap suites |
| PD07 | Help distinguishes explicit-seed restoration from `seed = NULL` advancing RNG. Regression new-data evaluation failures become `latents_bad_data`. | Prediction docs/implementation; plot/prediction suite |
| PD08 | Formerly base regression, pooled and bootstrap plots draw once when called, including inside loops, returning reusable ggplots invisibly. The return-type migration is explicit. | Plot helper and three methods; plot suite, drawing-integration suite, HTML previews |

Additional audited edges: named plots/report views walk the complete S3 inheritance chain; `intervals = FALSE` avoids inference/error bars; categorical-label coercion has a local warning handler; course-engagement prose describes continuous log-count-scale simulation followed by standardization/rounding. README and generated help were rebuilt from their sources.

## Scientific and compatibility decisions

Weighted residual associations are descriptive: a general sampling-weight design does not justify the previous unweighted reference distribution, so inference columns are withheld. Unweighted inference retains its established behavior.

Aggregated-binomial response/class-response predictions use expected-success units. For new data without response counts, `trials` accepts a scalar, row vector or column name. This argument follows `...`, preserving released positional arguments. Posterior exposures must agree with the observed outcome.

The knee retains the owner's calibrated distance. Tiny positive bends may yield a knee; corrected wording does not claim the full Kneedle algorithm or certainty about model selection. Active covariance floors can take longer than clipping; bounded common-orientation sweeps remain generalized EM updates. Growth reference tolerances accommodate established relative Cholesky jitter. Tests verify the reproduced cases and independent references, rather than proving correctness for every possible model/input/platform.

## Verification

- Source unit suite: **6,230 passed expectations / 931 cases**.
- Independent equivalence suite: **402 passed expectations / 62 cases**.
- Both suites: zero failures, errors, warnings or skips; **6,632 expectations / 993 cases total**, including 283 additional expectations since the audit.
- All 199 released help aliases and 14 actual plot-method signatures preserved; pkgdown reference checks pass.
- Configured lint: zero findings. Whitespace/diff check passes.
- README regenerated and updated enumeration article rendered; build/check logs record installed vignette verification.
- Native `R CMD build` succeeds. `R CMD check --as-cran --no-manual` completes with **zero errors, zero warnings, two environment notes**: DNS-restricted URL/repository checks and inability to verify current time. Package code/help checks, ordinary and `--run-donttest` examples, installed tests and vignette rebuilding pass.
- Installed CRAN-mode tests: **3,601 passed expectations, zero failures/warnings, 182 intentionally skipped cases**. The complete source suites above run with `NOT_CRAN=true` and exercise those local-only tests without skips.

`test-run.log`, suite CSV/RDS files and `SESSION.txt` record reproducible test/environment evidence. `help-check-final.log`, `help-alias-check.csv`, `lint.log`, `readme.log`, `article-render.log`, `build.log` and `check.log` record their respective gates. Native archive and installed-check output: `tmp/fixes-2026-10-13-build/`.

From the repository root, run `Rscript validation/fixes-2026-10-13/run-tests.R` and `Rscript validation/fixes-2026-10-13/check-documentation.R` to repeat the source and documentation gates.

New regressions: `test-audit-services-fixes.R` in both test directories; `test-audit-recent-api.R`, `test-audit-plots-predictions.R`, `test-review-core-audit.R`, `test-audit-drawing-integration.R` under `tests/testthat`; independent `test-review-growth-quadratic.R` under `tests/equivalence`.

Visual previews are retained: [plots/predictions](../../tmp/audit-plot-prediction-fixes.html), [structure curves](../../tmp/audit-fix-previews/smallest-structures.html), [pooled/bootstrap](../../tmp/audit-fix-previews/services.html), [enumeration article](../../tmp/audit-fix-previews/class-enumeration.html). The user subsequently requested “knit and show me”: both plot Rmd previews re-rendered successfully and all three plot previews were opened in the browser. Numerical/contract checks and presenting previews do not imply user visual acceptance. The user then requested the final repository handoff, now recorded in `HANDOFF.md`.

Specialist records: [core](../audit-2026-10-13/model-core/FIXES.md), [recent API](../audit-2026-10-13/recent-api/fix-closure.txt), [plots/predictions](../audit-2026-10-13/plots-docs/fixes.txt). References there to pending root integration are superseded by final verification here.
