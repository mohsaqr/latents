# Review 2026-11-04: inference-covariates

Probes: `validation/review-2026-11-04/probes/inference-covariates-p*.R`, all run through the two-slot wrapper on the working tree (1.1.15). Nothing in R/, man/ or tests/ was changed; fixes marked "verified" were patched in memory inside a probe only.

## 1. What was reviewed

| file | depth |
|---|---|
| R/covariate-inference.R | every line read; scores, Jacobian and SEs probed |
| R/covariates.R | likelihood, designs, M-step, fit flow and assembly read; likelihood rebuilt from scratch in base R. Print methods skimmed |
| R/inference.R | `parameter_inference.multilpa`, `.multilpa_response_bound`, `vcov`, `confint`, regularity, documentation read. Label helpers and `.multilpa_inference_jacobian` not re-derived: the 2026-11-01 audit checked them numerically and I relied on that |
| R/inference-service.R | every line read; sandwich, OPG, grand transform and p adjustment probed |
| R/robust.R | read; per-group scores compared row by row with numerical derivatives |
| R/bootstrap-inference.R | refit, replicate and table logic read and probed; plot layout read only |
| R/pool-imputations.R | read; Rubin's rules recomputed by hand; alignment probed |
| R/structure-inference.R | not re-derived; exercised through VEI, EVE, EEV and VVE fits |
| R/gaussian-diagnostics.R | read only (hand-checked by the 2026-11-01 audit); caller weighting checked by reading |
| R/simulation-service.R | read; cluster resampling probed. BLRT p-value only read |

Not reviewed: the `parameter_inference()` methods of the transition, growth and additive families, mixture-regression inference, `three_step()`, `bootstrap_lrt()` beyond the shared helpers, plots.

## 2. Findings

### F1. MAJOR, CONFIRMED. A held categorical block is treated as free in a covariate fit's information
`R/covariate-inference.R:790`, used at `:560`.

`.multilpa_cov_held_parameters()` returns `"response"` for `fixed = "response_probabilities"`. The test at `:560` compares it with the estimation-scale labels, where that block is called `"response_logit"`, so it never matches. The response logits therefore stay free in the Hessian, although the fit holds them. This path has no `free == n_parameters` check, which `parameter_inference.multilpa()` has (`R/inference.R`).

The documentation says the opposite: `?multilpa` (*Covariates*) calls this "the two-step estimator, whose criteria and standard errors condition on the held measurement", and `?parameter_inference` says the same.

Evidence, two-step LCA, n = 500, five binary items, weak separation (`p12-heldfix.R`, `p6-heldlca.R`). `n_parameters` is 2. The reference is a base-R likelihood of the two logits, with `numDeriv::hessian`:

| term | package SE | conditional SE by hand | ratio | one-step SE |
|---|---|---|---|---|
| (Intercept) | 0.2706 | 0.1383 | 1.96 | 0.2581 |
| z | 0.1815 | 0.1727 | 1.05 | 0.1803 |

- A second data set with Gaussian and categorical indicators (`p2-fixedcov.R`) gives ratios of 1.029 and 1.004. There the package SEs equal a hand Hessian with the 10 response logits free (0.1270933, 0.1451256), and the Gaussian block correctly held. That identifies the cause.
- `vcov(fit, scale = "unconstrained")` is 12 x 12 with the held logits in it, where it should be 2 x 2.
- The Hessian is also taken where the score of those logits is not zero, so a refusal for singular information is possible. I did not produce one.

Affected: `fixed = "measurement"` or `"response_probabilities"` with `profile_covariates` or `group_covariates` and any categorical indicator, under every `vcov_type`. Gaussian-only two-step fits are correct.

Smallest fix: `:790` becomes `c("response", "response_logit")`. Verified in memory: the SEs become 0.1382866 and 0.1727204, equal to the hand values, and the unconstrained vcov is 2 x 2.

**Changes reported numbers** for those fits; the owner decides.

### F2. MAJOR, CONFIRMED. `vcov_type = "opg"` under a grand transform breaks invariance
`R/inference.R:1165`, `R/covariate-inference.R:595`.

The added term for an estimated grand mean or SD uses `L = bread C diag(1/slope)`. That influence is the first-order expansion only when `bread` is the inverse observed information. With the OPG inverse in its place, `L` is non-zero for parameters the transform cannot affect.

Mixing probabilities and membership logits are exactly invariant to a grand location or scale change of the indicators, so their SEs must be too. Observed and robust keep this exactly (ratio 1 to printed digits); OPG does not.

Evidence (`p4-opg-rows.R`, `p3-matrix.R`):

| quantity (OPG SE) | no transform | grand centring | centring + scaling |
|---|---|---|---|
| group-class probability | 1.1752 | 1.5222 | 1.639 |
| profile probability, class 2 | 1.0169 | 1.3056 | 1.426 |
| profile probability, class 1 | 0.1544 | 0.1775 | 0.190 |

In a covariate fit the coefficient SE ratios against the untransformed fit are 0.986 to 1.331 (centred), 1.008 to 1.311 (scaled) and 1.002 to 1.475 (both).

Affected: OPG with `centering = "grand"` or `scaling = "grand"`, a narrow combination. The 2026-11-01 stacked check covered robust and observed only.

Smallest fix: pass `information$inverse` as `bread` for every type at both lines. Verified by hand on the covariate-free path (`p7-bound.R`): the probability SE goes from 0.04934 back to the untransformed 0.04199, and the variance rows equal the untransformed OPG values. The covariate path was not patched and rerun; the line there has the same form.

**Changes reported numbers** for that combination only. The stale comments at both lines ("the inverse the model's own covariance was formed with") and in the `.multilpa_transform_covariance` roxygen would need rewording.

### F3. MAJOR, CONFIRMED. A staged bootstrap with a measurement from other data gives intervals that do not describe the estimate
`R/bootstrap-inference.R:137` and `:196` onward; documented in `?fit_staged`.

`estimate` stays the fit's value, held at the supplied measurement, while every replicate re-estimates the measurement on the resampled current data. Nothing checks that the measurement came from these data, and nothing warns. The non-staged branch of the same function (`:199`) refuses a held fit for exactly this reason ("the held values came from another fit, and resampling these data does not resample them"). `?fit_staged` describes supplying a measurement "to reuse established profiles", which is the case where it comes from elsewhere.

Evidence (`p5-boot-pool.R`, `p10-units.R`), count of rows whose estimate lies outside its own percentile interval:

| measurement fitted on | rows outside | example |
|---|---|---|
| another sample, `y1` shifted by 3 | 11 of 18 | mean 1.450, interval (-1.670, -1.358) |
| another sample, same population | 3 of 14 | mean -1.713, interval (-2.107, -1.809) |
| same data (control) | 0 of 14 | |

Monte Carlo note: 20 and 30 replicates, so the same-population count is rough; the shifted case is unambiguous.

Smallest fix: refuse with `latents_unsupported_inference` when the stored first-stage data differ from the second stage's. `identical(stage_one$indicator_data, indicator_data)` printed FALSE here; I did not test that comparison on all same-data cases, for example a supplied same-data fit with rows in another order or an all-categorical fit. No estimate changes; the affected table becomes a refusal.

### F4. MINOR, CONFIRMED. `p_adjusted` depends on the origin of the indicators
`R/inference-service.R:624`; `R/inference.R` "bounded" rows; `R/covariate-inference.R` frame.

Means are tested against zero and enter the correction family with the slopes. Adding a constant to the indicators changes those tests, and so the adjusted p-value of a covariate slope, whose own p-value does not move. This is as documented ("taken over the tests the table reports") and `adjust` defaults to `"none"`, so it is a matter of what the number means rather than a broken formula.

Evidence (`p3-matrix.R`), indicators shifted by +1000, family of 11 tests of which 6 are means, slope `w` raw p = 0.1055 in both:

| method | original | shifted |
|---|---|---|
| Holm | 0.4221 | 0.3166 |
| BH | 0.1321 | 0.1189 |
| Bonferroni | 1 | 1 |

Smallest fix: leave `mean` rows out of the family, or give them no test, by the reasoning already applied to response probabilities. **Changes `p_adjusted`**, or `p_value` for means under the second option; the owner decides.

### F5. MINOR, CONFIRMED. The covariate path has no check that the score is small
`R/covariate-inference.R:383` onward, against `R/inference.R:1198`.

The covariate-free path warns `latents_unconverged` when the score is worth more than 1% of a standard error. The covariate path has no such check.

Evidence (`p7-bound.R`), largest |score| x SE, in % of an SE:

| tol | covariate fit | warned | plain fit | warned |
|---|---|---|---|---|
| 1e-8 (default) | 0.14 | no | 0.005 | no |
| 1e-4 | 6.2 | no | 0.12 | no |
| 1e-2 | 28.5 | no | 28.0 | yes |

Only loose `tol` is affected. Fix: the same warning. No number changes.

### EDGE
- Bootstrap with a rare ordinal level (`p9-rareord.R`): a resample that misses the level keeps the replicate with a diverging intercept, where a categorical item drops it. SE 11.07 against Wald 0.76, interval (-26.3, -0.63), 40 of 40 valid. The replicate fits do warn that a declared level is unobserved, and the caller sees that warning.
- Person centring with two group classes returned no table for any `vcov_type` in `p8-roworder.R`. The harness swallowed the error, so I did not read it. Between-group differences are removed by the centring, so a singular information refusal is the likely reason, but that is unconfirmed.
- Under `weights`, the grand offset is the unweighted mean of positive-weight rows (0.2074, weighted 0.2936). The stacked equations use the same definition, so the SEs are consistent (coefficient ratio 1). Whether that is the wanted offset is an interpretation question.

## 3. Verified sound

| item | check | result |
|---|---|---|
| One-step likelihood, reference class, logit coding | base-R likelihood from scratch, 3 profiles x 2 group classes, group covariate (`p1`) | -1416.24529003 both |
| Covariate observed, robust, OPG SEs | independent `numDeriv` Hessian and per-group Jacobian, with G/(G-1) | max rel. diff 2.1e-8, 4.4e-8, 8.6e-10 |
| Weighted covariate sandwich | same, by hand (`p10`) | 4.6e-9; default type robust |
| Covariate analytic scores | against `numDeriv` at a displaced point: full, FIML full, FIML + categorical, weighted slopes by class, ordinal + count, VEI, EEV, VVE, equal full, categorical pair (`p3`) | 4.5e-8 to 4.3e-6 |
| Delta Jacobian, covariate path | against numerical Jacobian, same configurations | at most 7.3e-9 |
| Group scores row by row, covariate-free | Jacobian of group log likelihoods: diag, full, FIML, weighted, mixed, VEI, VVE, local dependence, direct (`p4`) | at most 1.3e-6; previously only column sums were checked |
| Reference-class rebasing | permuted fit against hand contrasts of V (`p1`) | slope -0.6287 SE 0.16209; 0.7742 SE 0.16473, equal by hand; group slope sign flips |
| `boundary = "fix"` | base-R conditional Hessian, floor on last and on interior category, both paths (`p7`) | at most 5.3e-7 |
| Row-order invariance | shuffled, interleaved groups; 9 configurations under all types, plus weighted robust (`p8`) | SE ratio 1 to printed digits |
| Units and origin, covariate path | y x 1e4, y / 1e3, z x 1000 + 5e4, w / 100 - 7 (`p10`) | z statistics agree to 3e-9, mean SE scales exactly |
| Grand transforms, observed and robust | membership and mixing SEs invariant, weighted too (`p3`, `p4`, `p5`) | ratio 1 |
| Rubin's rules | by hand, m = 3, with a relabelled imputation in the run (`p5`) | estimate, SE, df, interval, p all differ by 0 |
| Structural constant in pooling | one group class (`p10`) | within = between = riv = fmi = 0, df Inf |
| Label alignment | `.multilpa_pool_align` undoes a permutation of profiles and group classes (`p5`) | coefficient difference 0 |
| Refit arguments | same data and arguments return the same fit: VVI, VEI, EVE, FIML, weighted, mixed, NB, local dependence, direct, grand, person (`p5`) | LL diff 0; counts and structures equal |
| Categorical pair bootstrap | runs through stored data (`p6`) | 8 of 8 valid |
| Resampling unit, failure accounting | read | groups redrawn as new groups, draws fixed before any refit, failures counted and warned |
| Refusals | prior and direct effects with covariates refused (`p3`, `p9`) | |

I found no defect in Rubin's formulas, the CR1 factor, the percentile intervals, or the logit and log interval scales.

## 4. Test gaps

Sabotage was applied in memory (`p11`, `p13`, `p14`), each against the files named.

| break | files run | result |
|---|---|---|
| Held response block (F1) | `test-covariate-fixed.R` by reading | Gaussian-only fixture; would pass |
| OPG under grand transform (F2) | grep | no test combines them |
| External measurement bootstrap (F3) | `test-staged-bootstrap.R` by reading | supplied-measurement tests check controls, not the interval |
| Response-simplex Jacobian uses wrong reference | `test-covariate-measurement`, `-inference`, `-structure-inference`, `test-categorical-inference`, `test-audit-1101-covariates-diagnostics` | **206 tests, 0 failed**. Binary items cannot show it (the sign flips, the variance does not). Other test files not run |
| Rubin `(1 + 1/m)` removed | `test-pool-imputations.R` | 0 of 39 failed, although its header announces Rubin's rules |
| same | `test-audit-1101-accessors.R`; `tests/equivalence/test-pool-imputations.R` | 4 of 49 failed (baseline 0); 10 of 40 failed. Protected, but not by the file named for it |
| Rubin df not squared | `test-pool-imputations.R`; accessors | 0 of 39; 4 of 49 |
| Pooled interval uses `qnorm`; fmi formula wrong | `test-pool-imputations.R` only | 0 of 39 each; other files not run with these two |
| CR1 factor removed | `test-cluster-robust-cr1.R` | 3 of 6 failed: caught |
| Covariate variance delta removed | `test-covariate-inference.R` | 1 of 40 failed: caught |
