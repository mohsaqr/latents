# Latents November review fixes

The repairs are implemented locally as **latents 1.1.17**. All fifteen major findings left open in the [5 November review](../review-2026-11-05/REPORT.md) are addressed by numerical corrections, truthful diagnostics or explicit refusals of unsupported inference. The documentation registry and test portability failures are repaired. Existing uncommitted work is preserved; nothing has been committed or published.

The existing 1.1.16 BH default for residual tables is retained. This session does not introduce that default. `adjust = "none"` continues to return raw p-values.

## Changes that affect results

| Area | Repaired behavior | Effect on results |
|---|---|---|
| Grand centering and scaling | Moments use positive sampling weights. Grand SD uses the population denominator for weighted and unweighted data. | Weighted constrained fits now agree with integer-expanded data, including missing cells. Unweighted grand-scaled likelihoods and scales also change slightly from the former sample-SD convention; refit when exact values matter. Person scaling retains sample SD. |
| Person centering | Any positive-weight group with exactly one observed value of an indicator is refused. | Artificial zero deviations can no longer support a spurious profile. The error identifies the affected indicators and groups; no observations are silently altered. |
| OPG inference | Transform sensitivity uses inverse observed information; grand nuisance equations use weighted population moments. | Probability SEs respect indicator translation. The original counterexample's roughly 8–12% OPG changes disappear. |
| Staged three-step correction | `step_one = "estimated"` refuses staged/measurement-held fits whose first-stage uncertainty is unavailable. | Incomplete corrections are no longer labelled as full estimated-measurement inference. Supported joint fits and conditional `step_one = "fixed"` remain available. |
| External staged bootstrap | Measurement provenance must match the multiset of fitted indicator observations, including missingness and duplicate multiplicity. | Previously misleading intervals are refused. Equivalent row orders and metadata names do not make same-sample measurement external. |
| Three-step membership regression | Separation/boundary checks precede pooled and transition-chain Wald inference. Private QR charts center and scale covariates. | Separated logits no longer produce extreme apparent significance. Affine predictor changes, including shifts by 1e8, preserve supported probabilities and uncertainty. |
| Regression factor levels | Unused levels are removed from membership, group-membership and random designs, retaining coding for prediction. | Unused levels no longer inflate parameter counts/BIC or cause spurious singular inference. |
| Identifiers | Output collisions use unique `id_` or `cluster_` prefixes across tables and their readers. | Original IDs survive when an identifier is called `class`; assignments can be joined back to data. Estimates are unchanged. |
| Growth identification | Covariance-design rank and repeated categorical capacity are checked before fitting. | Unsupported two-wave covariance decompositions are refused. Identifiable diagonal, staggered-occasion and proportional-covariance configurations have positive controls. |
| Transition identification | Known short-panel ridges and a necessary latent-path dimension bound are checked without covariates. | Arbitrary two-wave stayer/group-chain decompositions are refused. This is a scoped guard, not a general identification proof for covariate models. |
| Weak transition origins | Tables report `origin_expected_count` and `supported`; fewer than one expected outgoing move warns and withholds Wald inference. | Nonzero rows retain their formally free parameters and IC penalty. Exactly empty rows retain fixed/unestimated treatment. The diagnostic threshold signals weak information rather than deleting parameters at an arbitrary numerical tolerance. |
| Categorical bounds | A separate reporting predicate covers response and count bounds in checks, analytic-LRT warnings and candidate trust summaries. | Bounds no longer receive a clean diagnostic. Ordinary LCA's stored bootstrap flag is preserved; homogeneous transition probability flags are corrected to match the general engine. |
| Categorical bootstrap LRT | Simulation preserves exact joint-item and non-syntactic names. | Tests that previously failed every replicate can return calibrated bootstrap outputs; the new fixture has three usable draws out of three. Three draws verify execution, not p-value accuracy. |
| Ragged transition bootstrap | Schema checks precede alignment; extraction failures are counted and their reasons retained. | The reproduced crash becomes eight usable and seven explicitly rejected draws out of fifteen, with a warning that intervals condition on retained resamples. No coverage claim is made. |
| Cross-level prevalence | Both forms report the weighted posterior share of individuals. Fit tables and printing use positive-weight sample counts. | Full-model prevalences change on unequal groups; zero-weight groups/rows disappear from displayed sample counts. The warning against choosing form/class count from cross-level criteria remains. |

Under explicit `drop_empty = FALSE`, homogeneous and general transition engines now both use all positive-weight rows for individual BIC. `n_informative` retains its distinct meaning: positive-weight rows containing observed indicators. The corrected homogeneous BIC changes by `k * log(n_all / n_observed)`.

Standardized profile references use population moments for weighted and unweighted fits, so all-one weights and integer expansion agree. This also updates standardized plots. `descriptives()` continues to report unweighted recorded-row summaries and ICC, but now warns explicitly for a weighted fit; a new weighted ICC estimator has not been invented.

## Other verified repairs

- Categorical three-step adjustments count complementary binary hypotheses once under every supported multiplicity method. Unavailable modal overall tests warn; zero-SE profile/covariate Wald tests are unavailable.
- Second-order transition tables identify the first move at occasion 2 and document averaging over applicable occasions. Empty-origin counts use the correct moves. Trailing-dot terms retain names and SEs; isolated extreme predictor rows do not falsely mark finite interior fits as bound-active.
- Valid EII-to-EEI/EEE transition nesting is accepted; general `max_iter = 0` fits do not emit a false convergence warning.
- Correlated growth replication describes the fitted model rather than its initializer. One-class growth share intervals are `[1,1]`. Ordinal fitted/observed values, intercept ratio labels, row/group labels and recovery row-order checks are corrected.
- Classification gaps compare weighted modal and estimated shares. Residual details identify the smallest adjusted p-value accurately. Unconverged sensitivity refits are not ranked as a best optimum.
- Small coefficients retain digits beside large intercepts; predictors named `all` or `common` keep their names. Zero-weight missing Gaussian values are handled consistently with excluded categorical values.
- Articles correct swapped ESM class labels, bootstrap valid-draw denominators and withholding rules, AIC/entropy interpretations, LCGA bootstrap advice, a categorical plot example and residual-z wording. Numerical class descriptions have guards.

## Documentation and test repairs

All **44** documentation-registry rows remain. Twenty-three stale phrases now anchor surviving statements of the same behavior; reference metadata, backing tests, test names and notes are retained. Removed reference lists have not been restored merely to satisfy literal matching.

The retained Mplus STANDARDIZE comparison now tests population scaling directly, without the old n-versus-n-minus-one likelihood offset. Its likelihood and posterior checks pass.

Dispatch unit tests mock processor availability locally when testing a specified worker budget. Live multicore tests skip when the platform cannot report available processors. Unknown-core serial fallback remains independently tested. Stale BH, summary, schema and ratio-display expectations were updated to assert the corrected behavior.

Public help was regenerated from roxygen; all 77 Rd pages pass parsing and help checks. NEWS documents numerical changes and new refusals; README was knitted from its source. Detailed area reports describe implementation decisions and source paths:

- [Inference and registry](inference/REPORT.md)
- [Regression, growth and article corrections](regression/REPORT.md)
- [Transition repairs](transitions/REPORT.md)
- [Independent inference checks](independent/inference-review.md)

## Verification

| Validation | Passed | Failed / errors / warnings | Skipped |
|---|---:|---:|---:|
| Complete development sweep plus latest whole-file reruns: 258 files, 2,117 cases | 15,105 | 0 / 0 / 0 | 28 |
| Four new regression files, fresh final installed package | 249 | 0 / 0 / 0 | 0 |
| Documentation claims and retained Mplus STANDARDIZE comparison | 186 | 0 / 0 / 0 | 0 |
| `R CMD check --no-manual`, including examples, tests and vignette rebuild | **Status: OK** | 0 errors / 0 warnings / 0 notes | — |

The first full sweep completed with 15,064 passes, 18 failed expectations, one error and 16 warnings. These were stale sample-SD/diagnostic-text expectations, the newly refused unidentified setup, and expected boundary/cross-level warnings escaping older fixtures. Complete affected files were corrected and rerun; the table uses their latest outcomes. The original sweep is retained, and [per-file provenance](validated-by-file.csv) identifies each replacement. This is a complete sweep followed by whole-file rechecks, rather than a claim that the original snapshot passed unchanged. [Raw sweep](full-suite.log), [original conditions](full-suite-conditions.log), [reconciliation](reconciliation.log), [validated results](validated-results.csv).

Two production edits followed the full-sweep snapshot: inspection of nested general-LTA measurement bounds and qualification of `utils::tail`. Their affected diagnostics, enumeration, transition, identifier and regression files passed fresh installed-package rechecks. The checked source archive matches all 191 current R/help/namespace files byte for byte; DESCRIPTION source fields agree after normalizing build formatting. [Package check](root/00check-final.log), [source comparison](root/check-source-verification.log).

The 28 skipped cases are platform exclusions: 15 unsafe Accelerate fork tests, seven unknown-processor checks, four unavailable socket-port checks and two worker-limit checks. A separately rerun dispatch file exercises its mocked worker-budget tests without the CRAN core limit. Live multicore equivalence is not certified in this environment.

The new regression files include independent weighted-moment/replication checks, positive controls for supported configurations, invalid-input refusals and exact schema/identifier checks. Counts from overlapping targeted runs are not additive.

Independent checks found affine probability errors below 1e-9, a multi-class likelihood mapping error about 2e-9, weighted FIML nuisance equations below 1.5e-13 and nuisance replay error below 4.5e-16. These are numerical checks of the repaired equations, not new coverage simulations.

The previous review's counterexamples remain as dated evidence. New regressions assert the corrected behavior, and repair logs retain before/after probes. There was no attempt to make invalid Wald results appear valid by fabricating uncertainty or suppressing errors.

## Rendered previews

The [standardized-profile and categorical diagnostics preview](../../tmp/fixes-1105-preview.html) and [identifier/trajectory preview](../../tmp/regression-identifiers-1105.html) each contain two rendered images. They are retained for visual review. Article render statuses and links are listed in the regression report.

## Scope and remaining limits

The changes resolve the findings enumerated in the consolidated November review and several nearby verified minor defects. The older reviews' unconfirmed edge cases are not a claim of package-wide completeness. No new full Monte Carlo coverage study, entire external-software matrix or live worker stress test was performed. Expensive article bootstraps were not rerun merely to verify prose corrections; completed renders and lightweight previews are identified in the regression report.

