# Round-2 mathematical review: inference (2026-10-26)

## Scope

Files: `R/inference.R`, `R/inference-service.R`, `R/bootstrap-inference.R`,
`R/pool-imputations.R`, `R/sensitivity.R`, `R/alignment-service.R`,
`R/simulation-service.R`, `R/accessors.R`. Where a check had to follow a
primitive into its caller, I also read `R/robust.R` (`.multilpa_group_scores`),
`R/covariate-inference.R` (it calls `.multilpa_response_bound()`) and
`R/mixture-regression-inference.R` (`.mixture_invert`, the correlation-form
inversion).

Method: I loaded the package with `pkgload::load_all()` (NOT_CRAN=true, OMP and
VECLIB at one thread) and checked it against independent oracles:

- `numDeriv::grad` and `jacobian` of the package's own log likelihood and of
  per-group log likelihoods;
- `sandwich::vcovCL`;
- `mice::pool.scalar`;
- `lm`;
- exact 200-bit inverses (`Rmpfr`);
- Monte Carlo studies.

Scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-inference/`.
In the tables, `s*.R` names a script in that folder. No file under `R/`,
`tests/`, `equivalence/` or `man/` was edited.

A note on the oracle. `numDeriv::hessian` of a log likelihood of size about
2e3 cannot serve as the reference Hessian: its Richardson steps shrink to about
3e-5, where rounding dominates. It was 0.3% to 96% off the package (`s2b`).
Instead, the analytic score was first verified against `numDeriv::grad` and
central differences of the log likelihood. The reference Hessian is then
`numDeriv::jacobian` of that verified score. Against that reference, the
package's Hessian agrees to 5e-10.

## Findings

| id | sev. | file:line | what the code does | what it should do | reproduction | measured effect | proposed fix |
|---|---|---|---|---|---|---|---|
| INF2-1 | **high** | inference.R:930-933 and 1048-1050; the `centering` and `scaling` docs at multilpa.R:1009-1035 | Under `centering = "grand"` (and `scaling = "grand"`), the Wald path differentiates the likelihood of the transformed indicators. It treats the sample grand mean (and SD) as known constants, so the reported centred mean `mu_k - xbar` gets the standard error of `mu_k`. The estimates satisfy `sum_k pi_k mu_kj = 0` exactly (max 8.3e-15 over 300 fits), yet the reported covariance gives that identically-zero contrast a positive variance. | Report the SE of the reported estimate. With one group class and complete data, `xbar = sum_j pi_j mu_j` exactly at the MLE, so the delta method through `g(theta) = mu_k - sum_j pi_j mu_j` is exact to first order. In general, stack `xbar` (and the SD) as extra estimating equations and use the sandwich, or point to the bootstrap, which re-centres each resample. | `s12-centring-mc.R`, `s13-grand-fix.R` | Monte Carlo, 300 replications, 40 groups of 6, 2 profiles. Reported SE divided by the empirical SD of the reported mean: **0.779, 0.690, 0.789, 0.945** for grand centring; 0.910, 0.722, 1.004, 0.994 for grand centring plus scaling. With `centering = "none"` on the same design: 1.014, 0.983, 0.986, 0.974. The Monte Carlo SD of each ratio is 0.041. So SEs are up to 31% too small, and z tests of "profile mean = grand mean" are anti-conservative. | Apply the delta method on the centred-mean map. It restores the ratios to 0.961, 0.959, 0.993, 0.961 (`s13`); a 200-resample cluster bootstrap on one data set gave 0.94, 0.86, 0.90, 0.97. If that is not done, document that the SEs are conditional on the centring constants and steer users to `method = "bootstrap"`. |
| INF2-2 | **high** (default SEs); medium (refusal) | inference.R:1029-1031 (`.inference_require_rank` on the group scores), together with the default `vcov_type = "observed"` | **Observed SEs too small.** Under `centering = "person"`, the observed-information SEs treat the centred rows as independent draws. They are not: within a group they sum to zero, and each row's offset depends on the group's latent mix. **Robust path refused.** The robust path would correct this, but it is refused as "rank deficient". Every group's scores satisfy one exact linear relation per centred indicator `j`: `sum_k sigma2_kj s_g(mu_kj) + sum_{k<K} (mu_kj - mu_Kj) s_g(logit_k) = 0`. That relation is the per-group form of `sum_k pi_k mu_kj = 0`, so the meat is singular in exactly those directions, and correctly so. | Accept the sandwich, which never inverts `B`: keep the `G > p` cluster rule for `"robust"` and apply the rank refusal to `"opg"` only. Alternatively, exempt the known constraint directions. Then recommend `vcov_type = "robust"` (or make it the default) for person-centred fits. | `s10b-diag.R`, `s11-centred-mc.R`, `s14-person-sandwich.R` | The rank deficiency equals the number of centred indicators (2) in 300/300 fits, and `vcov_type = "robust"` was refused in 300/300. In the Monte Carlo (300 replications), observed SE divided by the empirical SD was **0.838, 0.777, 0.780, 0.929**. The refused CR1 sandwich, computed by hand, gave 0.967, 0.961, 0.909, 1.012 (Monte Carlo SD 0.041). On one fit, the identically-zero contrast was reported with SE 0.0858, against mean SEs of 0.096 and 0.200. | As in the "should do" column. |
| INF2-3 | medium | accessors.R:530 (profile table, `boundary = "error"`) vs accessors.R:706 (response table, `boundary = "fix"`). The same default holds at extra-indicators.R:730, 764 and direct-effects.R:385; those I read but did not run. | A fit where some category is unused by one class (round 1 called this the most common real-data case, INF-6) gets `NA` mean and variance SEs in `get_results(fit, "profiles")`, with a notice. `get_results(fit, "responses")` on the same fit reports conditional SEs. `get_results()` cannot pass `boundary`. | Use the same policy in every measurement table: `boundary = "fix"`, whose conditional covariance is verified below to 1e-9. | `s6-accessor-se.R` | Profile table: 4 of 4 rows `NA`. Response table: 5 of 6 rows have SEs; the held probability is `NA`. | Pass `boundary = "fix"` at accessors.R:530 (and at the three sibling calls). |
| INF2-4 | low | inference-service.R:53-54 (doc) vs 70-72 (code) | For `h_rule = "scaled"`, an unusable scale falls back to the width `1e-5 * (1 + abs(theta))`. | The doc says the fallback is the `"plus1"` width, `step * (1 + abs(theta))`. | `s7-service.R` | At `step = 1e-3`, `theta = 10`: the step used was 1.1e-4; the documented width is 0.011. Harmless with an analytic score. | Use `step * (1 + abs(value))`, or correct the doc. |
| INF2-5 | low | sensitivity.R:41-42 (`@param tolerance`: "Two log likelihoods within this of each other count as the same optimum") vs sensitivity.R:273-283 (single-linkage chaining) | A chain of maxima, each within `tolerance` of the next, is one optimum however far it spans. | Measure from the best value of the basin (leader clustering), or document the chaining in `@param tolerance`. | `s5-bootstrap-pool-sensitivity.R` | 30 log likelihoods 0.9e-4 apart span 2.6e-3, that is 26 times the tolerance of 1e-4. All 30 get `optimum = 1` and `best = TRUE`. | Start a new optimum when `value < leader - tolerance`. |
| INF2-6 | low | simulation-service.R:115-116 | The doc says "replicate `i` sees the same stream whatever the refit consumes". | In fact a refit drawing from the global stream (a refit gets no seed) shifts every later resample. Reproducibility under one `seed` is unaffected. | `s15-stream.R` | Replicate 1 is identical when the refit consumes 5 extra uniforms; replicate 2 is not. | Reword the doc, or give each replicate its own derived seed. |

Observations (not defects):

- **O1. The `"scaled"` step uses the conditional SE.** Its scale is
  `1 / sqrt(colSums(scores^2))`, the conditional (OPG-diagonal) standard error.
  For collinear coordinates this is far below the marginal SE: 1e-8 for a slope
  at origin 1e4. In a linear regression at origin 1e4, the SE error against an
  exact 200-bit inverse was 2.5e-5 for `"scaled"` and 5.9e-9 for `"plus1"`
  (`s7b`). Both are negligible. The scaled rule's benefit lies in nonlinear
  scores (MR-F6).
- **O2. INF-3 is still present.** Round 1 reported it as a by-design behaviour
  change, and REPORT.md does not list it. `parameter_inference(fit, data =
  read.csv(...))` and the bootstrap with CSV data are refused for every
  centring, while `sensitivity()` accepts the same CSV (`s7c-csv.R`). The two
  verbs disagree on what counts as the same data.
- **O3. A refusal that is correct.** On my first ordinal fit the ordinal item
  was a rounding of indicator `a`, so it was quasi-separated (intercepts
  45-69). The information matrix had two eigenvalues near 1e-8, and the
  refusal was correct (`s3b`).

## Verified correct

| Algorithm | Check | Max discrepancy |
|---|---|---|
| Analytic two-level score (`.multilpa_score`) for 22 fits: VVI, EEI, VVV, EEE, VEI, EVE and VVE charts; mixed and all-categorical (3- and 4-category); sampling weights (diagonal and VVV); ordinal plus Poisson count; negative binomial; direct effects; local dependence; FIML (continuous; continuous with categorical codes missing); person centring; grand scaling; single level | `numDeriv::grad` of the package log likelihood, and central differences (h = 1e-5), both at the estimate and 0.05 away (`s1`, `s1b`, `s3`, `s3c`, `s10`, `s10b`) | at most 9.3e-7 relative with numDeriv at d = 1e-4 (noise-limited), except FIML-categorical at 3.6e-5, which is 7.0e-8 absolute by central differences; at the best step, 6.9e-9 relative (VVI) |
| Observed information by `optimHess` of the scaled score (`.multilpa_observed_hessian`), with back/forward transform for direct effects | `numDeriv::jacobian` of the verified score (`s2`, `s3`, `s10`) | relative 4.2e-8 (EVE); 5.1e-10 to 1.5e-8 otherwise |
| Unconstrained covariance (inverse information) | same | 3.0e-7 (EVE), otherwise at most 1.4e-7 |
| Delta method to natural parameters (softmax rows, log variance, log-Cholesky, chart Jacobians, ordinal, count, NB, direct-effect shear) | `numDeriv::jacobian` of an independent decode-to-natural map, `J V J'` | SE relative 2.1e-7 |
| Robust CR1 `G/(G-1) A^-1 B A^-1` (weighted pseudo-likelihood scores `w_g s_g`; G = positive-weight, informative groups) | per-group `numDeriv` scores; `sandwich::vcovCL(type = "HC0", cadjust = TRUE)` (`s2`, `s8`) | vs numDeriv 1.4e-6 (EVE), at most 2.8e-7 otherwise; vs vcovCL at most 3.4e-8 (5 fits, including weighted and single level with G = 320) |
| OPG (BHHH) | per-group numDeriv scores | 2.2e-6 (EVE), at most 3e-7 otherwise |
| MLR scaling correction `tr(A^-1 B)/q` | hand computation | 2.7e-8 absolute |
| `.multilpa_response_bound()`, `boundary = "fix"`. Designs: one floored category in a 3- and a 4-category item; two floored categories including the last; a row with a single live category; a binary item. Each design was run in 3 category orders. | (A) numDeriv in the package's own free coordinates; (B) a model with the floored probabilities fixed at their estimates and the rest `(1 - fixed) * softmax`; (C) SEs compared across category orders (`s4`, `s4c`) | (A) 1.7e-9; (B) 5.6e-9; (C) 3.7e-10. The held-largest-logit rule is equivalent to holding `p_last` at the floor, to within O(floor). |
| `boundary = "fix"` with `vcov_type = "robust"` | conditional numDeriv sandwich (`s9`) | 5.1e-8 |
| Covariate family consuming `.multilpa_response_bound()` with labels `x:1`, `y, 2`, `z[3]` | held set and SE table (`s4b`) | names match; the held probability has SE `NA`; the others are finite |
| `fit_staged()` conditional covariance; staged bootstrap resamples both stages | numDeriv over free coordinates (`s9`) | 2.4e-9; `stages_resampled = 2` |
| Bootstrap percentile interval and SD | `quantile(type = 7)` and `sd` of the returned replicates (`s5`) | 2.2e-16 and 0 |
| Label alignment and relabelling | a known permutation (profiles `c(3,1,2)`, group classes `c(2,1)`) recovered (`s5`) | coefficients identical (0) |
| BLRT p-value `(1 + #{T* >= T})/(B + 1)` | hand computation | exact (0.5) |
| Rubin's rules, end to end (estimate, `T`, df `(m-1)(1+1/r)^2`, riv, fmi, t interval, pooled vcov) | `mice::pool.scalar(n = Inf)` on 26 rows with `B > 0` (`s5`) | qbar 4.4e-16, T 2.8e-17, df relative 8.4e-16, riv 4.3e-19, fmi 4.3e-19, interval 0; vcov diagonal equals SE squared |
| `get_results(, "profiles")` SE plumbing, raw and standardized (variance, covariance diagonal, `shared` outcome) | the `parameter_inference()` rows (13 fits, `s6`) | 0 |
| `p.adjust` family excludes untestable (`NA`) rows | hand computation (`s7`) | exact |
| `.inference_numeric_jacobian` with a scale | `exp` at scales 1e-4, 1, 1e3 | relative 2.0e-8 |
| `.inference_score_information` `"scaled"` and `"plus1"` | analytic normal-regression Hessian, exact 200-bit inverse, origins 0, 2000, 1e4 (`s7b`) | SE relative 2.5e-5 (scaled, origin 1e4), otherwise at most 9e-8 |
| Correlation-form inversion `.mixture_invert` | one-class Gaussian `mixture_regression` vs ML `lm` SEs, predictor origins 0, 2000, 1e5 (`s16`) | 1.4e-12, 5.3e-8, 1.6e-5 (condition-limited; the raw-matrix eigen ratio 3.3e-19 would have refused) |
| `.multilpa_optimum_index` on separated values | hand computation | exact |
| Related test files on the current tree | test-review-1026-inference, test-inference, test-pool-imputations, test-bootstrap-inference, test-sensitivity, test-review-1026-categorical | 72 tests, 391 passes, 0 failures or errors, 1 skip |

## Not checked, and why

- Bootstrap percentile coverage, and the bias from dropping failed resamples:
  this needs a Monte Carlo study of the bootstrap itself, which is beyond this
  round's budget. Only one-data-set ratios are reported, in INF2-1.
- The general-case remedy for INF2-1: with several group classes, unequal
  group sizes or FIML, `xbar` is not a function of the parameters alone. I
  verified the delta fix only for one group class with complete data.
- Whether person centring with 2 or more group classes also makes the group
  scores rank deficient: derived only for one group class. The INF2-2 Monte
  Carlo used one group class.
- `extra-indicators.R` and `direct-effects.R` measurement tables (INF2-3
  siblings): read, not run.
- Greedy alignment above 7 classes, and alignment under near-duplicate
  profiles: not stress-tested.
- The Mplus MLR claim (no `G/(G-1)` factor): Mplus is not installed.
- `pool_imputations()` for covariate fits and with a floored response in some
  imputations: not run. By reading, the pooling refuses such a fit, because
  the inference tables would differ.

## Fixes applied (same day, on the coordinator's decisions)

| id | change | files |
|---|---|---|
| INF2-1 | **Joint estimation.** The grand mean (and the grand standard deviation) are now stacked with the model's score as estimating equations. Each unit contributes `sum_i (x_ij - m_j)` and `sum_i ((x_ij - m_j)^2 - s_j^2 (n_j - 1)/n_j)`, over observed cells of carrying units. The added covariance is `c (bread B_SE L' + L B_ES bread) + c L B_EE L'`, with `L = bread C diag(1/slope)`. `C` is the derivative of the summed analytic score with respect to the transform, by central differences on re-transformed data. `c` is the CR1 factor for robust and 1 otherwise. This works for any number of group classes, under FIML, and with weights. **Negative variances.** A variance negative beyond 1e-8 of the model's own raises `latents_singular_information`; within rounding it is zero. | inference-service.R (`.multilpa_grand_transform`, `.multilpa_transform_covariance`); inference.R (called after the unconstrained covariance is formed) |
| INF2-2 | **Robust allowed.** The rank refusal now applies to OPG only; robust keeps the rule of more groups than parameters. **Person transforms.** Person-centred or person-scaled fits default to `"robust"`. `"observed"` and `"opg"` raise `latents_unsupported_inference`. `pool_imputations()` resolves the type the same way. | inference.R; inference-service.R (`.multilpa_transformed_vcov`); pool-imputations.R |
| INF2-3 | Every measurement table that reports SEs uses `boundary = "fix"`: profiles, ordinal, count and direct effects. | accessors.R; extra-indicators.R (2 calls); direct-effects.R (1 call) |
| INF2-4 | The documentation of the fallback width now matches the code. | inference-service.R |
| INF2-5 | Optimum counting uses the anchor rule: each optimum spans at most `tolerance`. Documented in `@param tolerance`. | sensitivity.R |
| INF2-6 | Every replicate's resample is drawn before the first refit. | simulation-service.R |
| INF-3 | **One identity rule.** A single rule (`.multilpa_indicator_reproduced`) is shared by the row-order check and by `.multilpa_inference_matrix`. A supplied column is accepted when it is exact, when the fit's own transform maps it onto the stored values, or when it is a CSV round trip of the original or of its reconstruction. The stored values are then used. **Direct effects.** Direct-effect covariates are also accepted as a CSV serialization. A changed value is still refused. | accessors.R; inference.R |

**Monte Carlo after the fix.** Values are the reported SE over the empirical SD, 200 replications, Monte Carlo SD 0.05 (`f1-transform-mc.R`). Each cell lists the four profile means.

| Design | Before | After |
|---|---|---|
| Grand centring | 0.800, 0.701, 0.733, 1.018 | 0.993, 0.975, 0.908, 1.048 |
| Grand centring and scaling | 0.917, 0.725, 0.925, 1.085 | 0.974, 0.974, 0.897, 1.042 |
| Grand scaling only | 0.521, 0.524, 0.584, 0.924 | 1.073, 1.085, 0.995, 1.045 |
| Grand centring, robust | not reported | 1.001, 0.974, 0.918, 1.051 |
| Grand centring, two group classes | 0.623, 0.509, 0.552, 0.810 | 1.050, 0.988, 0.918, 1.058 |
| Grand centring, FIML | 0.815, 0.683, 0.750, 0.977 | 0.997, 0.946, 0.921, 0.993 |

**Tests.** `tests/testthat/test-review-1026-inference-r2.R` has 8 tests with 33 expectations, all passing. All 8 fail on the pre-fix code. I checked that on a scratch copy of the tree with each of my edits reversed and `R/inference.R` as it stood before the patch (`old/latents`): 2, 1, 2, 0, 3, 2, 1 and 0 failures, and 0, 0, 1, 1, 0, 0, 0 and 1 errors, per test.

Three existing tests asserted the old contracts and were updated:
- `test-sensitivity.R`, the chained optimum;
- `test-audit-inference-repairs.R`, CSV refused by `parameter_inference()`;
- `test-centering.R`, observed SEs on a person-centred fit with 20 groups and 21 parameters, which now expects `latents_too_few_groups`; inference is checked on 40 students.

### Covariate fits (follow-up)

INF2-1 and INF2-2 are now applied in `R/covariate-inference.R`, reusing the
helpers.

- **INF2-1 (grand transforms).** `.multilpa_transform_covariance()` takes
  `score_of` and `group_scores_of`, so the covariate parameter vector (with its
  `moves` coordinates) gets the same stacked-equation term. The bread is the
  inverse information, or the OPG inverse.
- **INF2-2 (person transforms).** Person-transformed fits default to
  `"robust"`, and `"observed"` and `"opg"` are refused. This applies in both
  `parameter_inference()` and `vcov()`. The rank refusal now applies to OPG
  only.

Monte Carlo, 200 replications, mean over the four profile means of SE divided
by the empirical SD (Monte Carlo SD 0.05; `f3-covariate-mc.R`). The "before"
column was run on the old-code copy.

| Design | Before | After |
|---|---|---|
| No transform | 0.979 | 0.979 (unchanged) |
| Grand centring | 0.855 | 1.018 |
| Grand scaling | 0.609 | 0.983 |
| Person centring | 0.938 (observed); the sandwich was refused (`latents_too_few_groups`) | 1.086 (robust, now the default) |

Two calibration tests were added: 60 replications each, `skip_on_cran`. Both
fail on the old code.
