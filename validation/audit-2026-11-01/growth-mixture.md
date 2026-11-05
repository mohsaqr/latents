# Audit 2026-11-01 — slice `growth-mixture`

Files: R/growth-mixture.R, R/growth-mixture-engine.R, R/growth-mixture-correlated.R,
R/growth-mixture-inference.R, R/growth-mixture-methods.R,
R/growth-mixture-parameter-inference.R, R/growth-mixture-plot.R,
R/growth-mixture-predict.R, R/growth-mixture-simulate.R, R/trajectory-views.R.
133 functions. Probe scripts (synthetic data, one R process, 1 thread) are in the
session scratchpad (`gm/p1.R` … `gm/p11.R`); their key outputs are quoted below.

Reference checks run (all passed, no change needed):

- One-class marginal log likelihood (random intercept + slope, 10% rows dropped):
  latents -388.620609285, `nlme::lme(method = "ML")` -388.620609291, hand-integrated
  normal at latents' estimates -388.620609285. Fixed effects and VarCorr agree to 4+ digits.
- Analytic scores vs `numDeriv::grad`, and stored vcov vs `solve(-numDeriv::hessian)` in the
  data's coordinates (so the standard-coordinate Jacobian is checked too), for
  varying / equal+equal-variance / proportional / diagonal covariances with a common
  covariate and a membership covariate: max relative gradient error 8.6e-7, max
  relative SE error 1.7e-7.
- Correlated residuals (AR and Brownian): fitted LL equals a hand dense likelihood
  (mixture over classes, `Z G Z' + W + sigma^2 I`) on scrambled rows.
- `predict(..., type = "class_response")` equals `X beta` (+ offset) exactly;
  `predict(newdata = fitted data)` equals the fitted-data prediction for
  conditional, posterior and random_effects types (uncorrelated, AR, Brownian,
  multilevel, offset, factor, no-intercept random design): differences 0.
- Predicted random effects equal `G Z' V^-1 (y - X beta)` by hand (4e-16); trajectory
  table equals `X beta` with common covariates at their mean (9e-16); individual
  table rows align with scrambled character ids.
- `simulate()` then refit recovers coefficients, residual and random variances.
- Growth bootstrap tables carry `replicate_identity` (block, class, term, estimate) and
  `as.data.frame(what = "replicates")` works; a reordered copy is refused. A Wald/bootstrap
  share-SE ratio of 0.62 (30 replicates) was traced to the drawn resamples: the SD of
  the true class composition over those exact draws is 0.025112, equal to the bootstrap SE.

## Functions

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.growth_correlation_names` | R/growth-mixture-correlated.R:15 | ran (ar, brownian) | ok |
| `.growth_is_correlated` | R/growth-mixture-correlated.R:24 | ran | ok |
| `.growth_process_covariance` | R/growth-mixture-correlated.R:34 | ran; dense LL equals hand AR / Brownian likelihood (p7) | ok |
| `.growth_dense_patterns` | R/growth-mixture-correlated.R:50 | ran (scrambled row order, p7) | ok |
| `.growth_dense_covariances` | R/growth-mixture-correlated.R:67 | ran (p7) | ok |
| `.growth_dense_terms` | R/growth-mixture-correlated.R:86 | ran; LL vs hand dense likelihood, ar and brownian (p7) | ok |
| `.growth_dense_scores` | R/growth-mixture-correlated.R:138 | covered by test-growth-correlated (numerical gradient) | ok |
| `.growth_correlated_bounds` | R/growth-mixture-correlated.R:214 | read; ran (p7) | ok |
| `.growth_correlated_refit` | R/growth-mixture-correlated.R:251 | ran ar + brownian (p7) | logged |
| `.growth_newton_finish` | R/growth-mixture-correlated.R:348 | read; ran in every fit | ok |
| `.growth_residual_time` | R/growth-mixture-correlated.R:431 | read; covered by test-growth-correlated refusals | logged |
| `.growth_statistics` | R/growth-mixture-engine.R:23 | ran via LL checks (p1,p7) | ok |
| `.growth_patterns` | R/growth-mixture-engine.R:64 | ran; unbalanced (10-15% rows dropped) data in every probe | ok |
| `.growth_coefficients` | R/growth-mixture-engine.R:80 | ran | ok |
| `.growth_residuals` | R/growth-mixture-engine.R:89 | ran via LL checks | ok |
| `.growth_direct_residuals` | R/growth-mixture-engine.R:123 | read; covered by test-review-1026-growth (large-location cases) | ok |
| `.growth_class_terms` | R/growth-mixture-engine.R:158 | ran; one-class LL vs nlme::lme ML (diff 6e-9) and hand normal (p1) | ok |
| `.growth_covariances` | R/growth-mixture-engine.R:202 | ran (proportional scales, p2) | ok |
| `.growth_factor` | R/growth-mixture-engine.R:225 | ran | ok |
| `.growth_floor_covariance` | R/growth-mixture-engine.R:253 | read; ran via EM | ok |
| `.growth_structure_terms` | R/growth-mixture-engine.R:270 | ran | ok |
| `.growth_expectation` | R/growth-mixture-engine.R:294 | ran; LL vs hand (p1,p7) | ok |
| `.growth_maximization` | R/growth-mixture-engine.R:318 | read (closed forms checked against ECM derivation incl. proportional); EM monotone in test-growth-mixture | ok |
| `.growth_structure` | R/growth-mixture-engine.R:446 | ran (observation and two-level) | ok |
| `.growth_structure_names` | R/growth-mixture-engine.R:461 | ran | ok |
| `.growth_observed` | R/growth-mixture-engine.R:475 | read | ok |
| `.growth_em` | R/growth-mixture-engine.R:483 | ran | ok |
| `.growth_names` | R/growth-mixture-inference.R:11 | ran; names == pack names, length == n_parameters (p2) | ok |
| `.growth_pack` | R/growth-mixture-inference.R:37 | ran; round-trip in test-growth-mixture | ok |
| `.growth_unpack` | R/growth-mixture-inference.R:59 | ran | ok |
| `.growth_scores` | R/growth-mixture-inference.R:106 | ran; analytic vs numDeriv gradient, 4 configs, max rel err 8.6e-7 (p2) | ok |
| `.growth_quasi_newton` | R/growth-mixture-inference.R:278 | read; ran in every fit | ok |
| `.growth_inference` | R/growth-mixture-methods.R:23 | ran; vcov vs solve(-numDeriv Hessian), 4 configs (p2) | ok |
| `.growth_unreported_held` | R/growth-mixture-methods.R:57 | ran (AR boundary fit, random table NA SE, p7) | ok |
| `.growth_inference_at` | R/growth-mixture-methods.R:71 | ran (observed); robust/opg covered by test-growth-mixture | ok |
| `.growth_standard_jacobian` | R/growth-mixture-methods.R:137 | ran; via p2 data-coordinate vcov check | ok |
| `.growth_resolve_inference` | R/growth-mixture-methods.R:202 | ran | ok |
| `.growth_class_names` | R/growth-mixture-methods.R:216 | ran | ok |
| `.growth_delta_table` | R/growth-mixture-methods.R:224 | ran | ok |
| `.growth_coefficient_table` | R/growth-mixture-methods.R:240 | ran; vs nlme fixef (p1) | ok |
| `.growth_random_table` | R/growth-mixture-methods.R:254 | ran; vs nlme VarCorr (p1), process rows (p7) | ok |
| `.growth_class_is_degenerate` | R/growth-mixture-methods.R:318 | ran (p7) | ok |
| `.growth_class_table` | R/growth-mixture-methods.R:325 | ran; shares sum 1, expected persons sum n, weighted avg posterior (p3,p8) | ok |
| `.growth_class_prior` | R/growth-mixture-methods.R:368 | ran | ok |
| `.growth_membership_table` | R/growth-mixture-methods.R:383 | ran (p3, p6) | ok |
| `.growth_entropy` | R/growth-mixture-methods.R:405 | ran | ok |
| `.growth_random_description` | R/growth-mixture-methods.R:411 | ran | ok |
| `.growth_entropy_sum` | R/growth-mixture-methods.R:421 | ran | ok |
| `.growth_fit_table` | R/growth-mixture-methods.R:427 | ran; bic == -2LL + log(n) k (p3); zero weights (p8) | ok |
| `.growth_joint_entropy` | R/growth-mixture-methods.R:467 | read; covered by test-multilevel-growth | ok |
| `.growth_cluster_fit_table` | R/growth-mixture-methods.R:491 | ran (p6) | ok |
| `.growth_group_class_table` | R/growth-mixture-methods.R:521 | ran (p6) | ok |
| `.growth_cluster_table` | R/growth-mixture-methods.R:560 | ran (p6) | ok |
| `.growth_assignment_table` | R/growth-mixture-methods.R:580 | ran; ids, scrambled character ids (p3) | ok |
| `.growth_random_effects_table` | R/growth-mixture-methods.R:598 | ran; BLUP vs hand (p3) | ok |
| `.growth_time_variable` | R/growth-mixture-methods.R:616 | ran; bad time classed (p11) | ok |
| `.growth_design` | R/growth-mixture-methods.R:639 | ran (p9) | ok |
| `.growth_trajectory_table` | R/growth-mixture-methods.R:654 | ran; vs hand (p3, p9) | ok |
| `.growth_individual_table` | R/growth-mixture-methods.R:673 | ran; id/time/observed aligned with scrambled data, class mean vs hand (p3) | ok |
| `.growth_tables` | R/growth-mixture-methods.R:701 | ran | ok |
| `.growth_table` | R/growth-mixture-methods.R:708 | ran (all tables, 6 configs, p9) | ok |
| `.growth_present` | R/growth-mixture-methods.R:738 | ran (printing in p3/p6) | ok |
| `get_results.latents_growth_mixture` | R/growth-mixture-methods.R:850 | ran (what = all) | ok |
| `as.data.frame.latents_growth_mixture` | R/growth-mixture-methods.R:869 | ran | ok |
| `print.latents_growth_mixture` | R/growth-mixture-methods.R:877 | ran (p6) | ok |
| `summary.latents_growth_mixture` | R/growth-mixture-methods.R:927 | ran (p3) | ok |
| `print.summary_latents_growth_mixture` | R/growth-mixture-methods.R:940 | read; covered by test-growth-mixture | ok |
| `get_results.summary_latents_growth_mixture` | R/growth-mixture-methods.R:950 | read; covered by test-api-compat-0921 | ok |
| `as.data.frame.summary_latents_growth_mixture` | R/growth-mixture-methods.R:957 | ran (p3) | ok |
| `coef.latents_growth_mixture` | R/growth-mixture-methods.R:965 | ran | ok |
| `vcov.latents_growth_mixture` | R/growth-mixture-methods.R:970 | ran (p2) | ok |
| `confint.latents_growth_mixture` | R/growth-mixture-methods.R:981 | ran (p3) | ok |
| `logLik.latents_growth_mixture` | R/growth-mixture-methods.R:1000 | ran (p3) | ok |
| `nobs.latents_growth_mixture` | R/growth-mixture-methods.R:1007 | ran; zero-weight persons excluded (p8: 140 of 150) | ok |
| `.growth_weighted_persons` | R/growth-mixture-methods.R:1017 | ran (p8) | ok |
| `.growth_natural_parameters` | R/growth-mixture-parameter-inference.R:18 | ran (wald table, p4) | ok |
| `.growth_permute_params` | R/growth-mixture-parameter-inference.R:121 | read; ran via bootstrap | ok |
| `.growth_align_refit` | R/growth-mixture-parameter-inference.R:145 | ran (bootstrap p4,p5) | ok |
| `.growth_match_classes` | R/growth-mixture-parameter-inference.R:157 | ran | logged |
| `.growth_refit_arguments` | R/growth-mixture-parameter-inference.R:169 | read; ran via bootstrap | ok |
| `.growth_bootstrap_inference` | R/growth-mixture-parameter-inference.R:194 | ran iter 6 and 30; SE of share equals SD of true resample composition (0.02511 both, p5); replicate_identity present (p4, audit test) | logged |
| `parameter_inference.latents_growth_mixture` | R/growth-mixture-parameter-inference.R:260 | ran wald + bootstrap (p4,p5), boundary = fix (p7) | ok |
| `.growth_class_key` | R/growth-mixture-plot.R:14 | ran | ok |
| `.growth_scales` | R/growth-mixture-plot.R:35 | ran | ok |
| `plot.latents_growth_mixture` | R/growth-mixture-plot.R:48 | ran all views, 6 configs (p9); bad arguments (p11) | fixed |
| `.growth_observed_means` | R/growth-mixture-plot.R:80 | ran; normalised vs raw weights (p10) | ok |
| `.growth_spread` | R/growth-mixture-plot.R:111 | ran (p9) | ok |
| `.growth_plot_trajectories` | R/growth-mixture-plot.R:130 | ran (p9, p11, render) | fixed |
| `.growth_plot_individuals` | R/growth-mixture-plot.R:231 | ran; 5-class LCGA, persons = 12 (audit test, render) | fixed |
| `.growth_plot_coefficients` | R/growth-mixture-plot.R:295 | ran (p9) | ok |
| `.growth_plot_random` | R/growth-mixture-plot.R:343 | ran 1 and 2 random effects (p9) | ok |
| `.growth_plot_classification` | R/growth-mixture-plot.R:419 | ran weighted growth and LCGA fits (p8, render) | fixed |
| `.growth_new_spec` | R/growth-mixture-predict.R:18 | ran (p6, p7, p9) | ok |
| `.growth_random_design` | R/growth-mixture-predict.R:62 | ran | ok |
| `.growth_class_means` | R/growth-mixture-predict.R:84 | ran; == X beta (audit test, p1) | ok |
| `predict.latents_growth_mixture` | R/growth-mixture-predict.R:93 | ran all six types; newdata = fitted data identical (p7, p9) | ok |
| `.growth_draw` | R/growth-mixture-simulate.R:14 | ran; recovery (p4, audit test); correlated (p7); multilevel (p6) | ok |
| `simulate.latents_growth_mixture` | R/growth-mixture-simulate.R:67 | ran | ok |
| `.growth_check_arguments` | R/growth-mixture.R:8 | read; refusals ran via mixture_regression (weights+cluster, common=~1) and test-growth-mixture | ok |
| `.growth_cluster_spec` | R/growth-mixture.R:70 | ran multilevel fit + predict(newdata) (probe p6); test-multilevel-growth | ok |
| `.growth_start` | R/growth-mixture.R:114 | read; ran via every fit (p2,p3,p9 configs) | ok |
| `.growth_two_level_start` | R/growth-mixture.R:172 | read; ran multilevel fit (p6) | ok |
| `.growth_trajectory_seed` | R/growth-mixture.R:208 | read; ran (starts table kind 'trajectories', p3) | ok |
| `.growth_carrying_spec` | R/growth-mixture.R:251 | ran weighted fit with zero-weight persons (p8) | ok |
| `.growth_fit` | R/growth-mixture.R:265 | ran 4 covariance configs + 6 edge configs (p2,p9); starts table inspected (p3) | ok |
| `.growth_standardize` | R/growth-mixture.R:377 | ran; data-unit vcov equals numDeriv Hessian in data units (p2, max rel SE err 1.7e-7) | ok |
| `.growth_apply_standard` | R/growth-mixture.R:420 | ran via predict(newdata) == fitted (p7, p9: diff 0) | ok |
| `.growth_random_map` | R/growth-mixture.R:452 | ran; covered by p2 vcov and BLUP check (p3, err 4e-16) | ok |
| `.growth_unstandardize_coefficients` | R/growth-mixture.R:473 | ran; trajectory table vs hand X beta incl. common, offset, factor, no-intercept random (p3,p9: err <= 9e-16) | ok |
| `.growth_unstandardize_params` | R/growth-mixture.R:510 | ran; LL vs nlme and hand integration (p1) | ok |
| `.growth_person_patterns` | R/growth-mixture.R:532 | ran via correlated predict (p7) | ok |
| `.growth_map_expectation` | R/growth-mixture.R:556 | ran; BLUP vs hand G Z' V^-1 r (p3, 4e-16) | ok |
| `.growth_unstandardize_fit` | R/growth-mixture.R:590 | ran; LL equals hand-integrated LL (p1, p7) | ok |
| `.growth_standardize_params` | R/growth-mixture.R:619 | read; covered by test-review-1026-growth-r2 | ok |
| `.growth_standard_params` | R/growth-mixture.R:658 | ran via inference (p2) | ok |
| `.growth_fitted_expectation` | R/growth-mixture.R:677 | ran; predict(newdata = fitted data) posteriors/RE identical (p7) | ok |
| `.growth_standard_view` | R/growth-mixture.R:692 | ran via vcov (p2) | ok |
| `.growth_select` | R/growth-mixture.R:701 | ran; starts table, warnings (p3, p7) | ok |
| `.growth_degenerate_covariances` | R/growth-mixture.R:799 | ran (AR fit flagged class_1 at floor, p7) | ok |
| `.growth_order_classes` | R/growth-mixture.R:822 | ran (proportional config, stopifnot invariant held, p2/p9) | ok |
| `.growth_count_parameters` | R/growth-mixture.R:846 | ran; equals length(coef) in 10 configs (p2,p9) | ok |
| `.trajectory_view` | R/trajectory-views.R:17 | ran (growth and LCGA) | ok |
| `.trajectory_inference` | R/trajectory-views.R:46 | ran | ok |
| `.trajectory_observed` | R/trajectory-views.R:53 | ran | ok |
| `.trajectory_band` | R/trajectory-views.R:72 | ran; trajectory table vs hand (p3,p9) | ok |
| `.trajectory_coefficients` | R/trajectory-views.R:100 | ran | ok |
| `.trajectory_is_fit` | R/trajectory-views.R:104 | read | ok |
| `.trajectory_check_three_step` | R/trajectory-views.R:120 | read; covered by test-growth-r3step | ok |
| `.trajectory_check_data` | R/trajectory-views.R:145 | read; covered by test-growth-parameter-inference | ok |
| `.trajectory_pieces` | R/trajectory-views.R:180 | read; weighted fits refused upstream (three-step.R:502) | ok |
| `.trajectory_uses` | R/trajectory-views.R:188 | read | ok |
| `.trajectory_three_step_inputs` | R/trajectory-views.R:203 | read; covered by test-growth-r3step, test-three-step-step-one | ok |
| `.trajectory_r3step_inputs` | R/trajectory-views.R:232 | read; covered by test-growth-r3step | ok |

## Fixed

1. **Classification plot reported unweighted average posteriors, including persons of
   weight zero** (`R/growth-mixture-plot.R:453-459`, `.growth_plot_classification`, used by
   growth and latent class growth fits). Weighted growth fit (150 persons, 10 at weight 0):
   `get_results(fit, "classes")$average_posterior` = 0.9424, 0.9129, but the plot subtitle said
   "Class 1 0.94, Class 2 0.92" (unweighted mean over all persons: 0.9392, 0.9179). For the
   latent class growth fit of the same data the table said 0.955 / 0.923 and the plot
   0.95 / 0.94. Cause: `mean(certainty[modal == k])` ignored the sampling weights the
   tables use (family convention: weighted AvePP over carrying units). Change: weighted mean
   over persons with positive weight; persons of weight zero keep their bar. Unweighted fits
   are unchanged. Test: "the classification plot reports the weighted average posterior of
   the tables" (hand-computed weighted mean). Broken on purpose: fails.
2. **Individuals view dropped whole classes** (`R/growth-mixture-plot.R:241-256`,
   `.growth_plot_individuals`). With 5 classes and the default `persons = 12`, 3 persons were
   picked per class in class order and `head(, 12)` cut all of class 5 (7 classes: same).
   Change: the per-class picks are taken in turns across classes before the cut, so every
   class is thinned evenly; when no cut happens the selected set is identical. Test: "the
   individuals view shows every class when persons must be cut". Broken on purpose: fails
   (class_5 missing).
3. **Plot arguments: unclassed errors and silent NaN bands** (`R/growth-mixture-plot.R:59-63`,
   `:134-137`, `:237-239`). `plot(fit, level = 1.5)` drew NaN bands with a base "NaNs produced"
   warning; `level = NA` drew silently; `spread = NA`, `facet = "yes"` gave unclassed
   `simpleError`; `max_persons = NA` gave "missing value where TRUE/FALSE needed";
   `persons = -3` gave "'length.out' must be a non-negative number". Change: `.latents_require()`
   (class `latents_bad_argument`) for spread/facet in the method, level and max_persons in
   `.growth_plot_trajectories()` and persons in `.growth_plot_individuals()` — the helpers so
   the latent class growth plot path (R/mixture-regression-predict.R) is covered too. No
   argument added. Test: "growth plot arguments outside their range are refused with a
   class". Broken on purpose: fails.

Visual output: `tmp/audit-1101-growth-mixture.html` (5-class individuals view, weighted
classification view, unchanged trajectories view).

## Logged, not fixed

- `.growth_bootstrap_inference` (and the profile/LTA bootstraps in other files) resamples persons of
  weight zero too, so the number of carrying persons varies across replicates; consistent across
  families, a design decision, not changed here.
- `.growth_match_classes` matches bootstrap classes using every original person, including zero-weight
  ones; harmless for matching.
- `.growth_correlated_refit` keeps the uncorrelated fit's `starts` table and `n_best_replicated`, and a
  `latents_random_boundary` warning raised by the uncorrelated stage is not re-judged after the process
  is added.
- `.growth_residual_time`: the Brownian negative-time check uses `na.rm = TRUE` on raw data; NA times on
  fitted rows are refused later (`latents_bad_data`), so no NA is hidden.
- `plot(what = "individuals", persons = 5)` reads a single number as a count, never as the identifier 5
  (documented: "a number or identifiers").
- `get_results(time = )` accepts any numeric column of the model data, including a numeric `id`.
- Bootstrap `covariance` attribute of growth tables has no dimnames (profile tables have them).
- `.growth_class_table` / `.growth_coefficient_table` use `sqrt(pmax(diag(vcov), 0))`: a negative
  diagonal (not seen on realistic data) would show SE 0 rather than NA.

## Needs another slice

- none found. (Weighted trajectory fits are refused by three_step()/r3step() at R/three-step.R:502,
  :1113, so `.trajectory_pieces()` not carrying `weights`/`kept` is not reachable with weights.)

## Test runs

```
test-audit-1101-growth-mixture.R: PASS 24 FAIL 0 SKIP 0 ERROR 0
test-growth-mixture.R: PASS 122 FAIL 0 SKIP 0 ERROR 0
test-growth-correlated.R: PASS 15 FAIL 0 SKIP 0 ERROR 0
test-growth-predict.R: PASS 19 FAIL 0 SKIP 0 ERROR 0
test-growth-parameter-inference.R: PASS 39 FAIL 0 SKIP 0 ERROR 0
test-growth-r3step.R: PASS 12 FAIL 0 SKIP 0 ERROR 0
test-multilevel-growth.R: PASS 33 FAIL 0 SKIP 0 ERROR 0
test-review-1026-growth.R: PASS 25 FAIL 0 SKIP 0 ERROR 0
test-review-1026-growth-r2.R: PASS 67 FAIL 0 SKIP 0 ERROR 0
test-review-1101-bootstrap-tables.R: PASS 67 FAIL 0 SKIP 0 ERROR 0
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 ERROR 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 ERROR 0
test-weights.R: PASS 76 FAIL 0 SKIP 0 ERROR 0
test-mixture-regression.R: PASS 73 FAIL 0 SKIP 0 ERROR 0
test-three-step-step-one.R: PASS 172 FAIL 0 SKIP 0 ERROR 0
test-audit-1101-growth-mixture.R (rerun after final edit): PASS 24 FAIL 0 SKIP 0 ERROR 0
```

devtools::document(): no NAMESPACE or man/ change.

equivalence/engine-golden/compare.R: identical 247 / 260, DIFFERENT 13. All 13 are bootstrap tables whose attribute `replicate_rows` is now `replicate_identity` (the R/bootstrap-inference.R change, outside this slice); all 17 growth cases are identical.
