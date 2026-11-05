# Audit 2026-11-01 — slice `mixture-regression`

Files: R/mixture-regression.R, R/mixture-regression-engine.R,
R/mixture-regression-enumerate.R, R/mixture-regression-inference.R,
R/mixture-regression-methods.R, R/mixture-regression-ordinal.R,
R/mixture-regression-predict.R, R/mixture-regression-starts.R
(134 functions and one constant, `.ordinal_min_log_gap`).

Checking scripts (synthetic data, run with `OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`, `devtools::load_all()`) are in the
session scratchpad
/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/mixreg/
(not kept in the repository): `c1.R` (Gaussian, Poisson, binomial M-steps vs
`lm()`/`glm()` with posterior weights, common terms vs stacked WLS, NB
dispersion vs `optimize()`), `c2.R` (one-class ordinal vs `MASS::polr`,
hand likelihood, category probabilities), `c3.R` (one-class random intercept
vs `nlme::lme` ML), `c4.R`/`c4b.R` (standard errors vs the inverse
`numDeriv` Hessian for every family and nesting; robust vs a hand sandwich),
`c5.R` (predict vs design times coefficients and prior-weighted mixtures),
`c6.R` (zero-weight invariance of every table in four nestings), `c7.R`,
`c15.R` (methods, refusals and their classes), `c8.R` (SE invariance to
shifted/scaled predictors), `c9.R` (enumeration and BLRT), `c10.R`,
`c17.R` (ordinal predictions, trajectory band vs numDeriv delta), `c11.R`
(two-level, permutation invariance, empty-table schemas), `c16.R` (robust
clustered on `id`, weighted), `c12.R` (growth bootstrap table attribute),
`c18.R` (ordinal with every slope shared: the crash), `c19.R`/`c20.R`
(start-stage census).

Slice-note verifications (all run, numbers from this session):

- Gaussian M-step vs `lm(..., weights = posterior)`: one step from zero
  differs by 7.3e-10 / 1.3e-9 (the 1e-10 ridge, see Logged); the weighted
  least-squares solution is the step's exact fixed point (1e-10 test);
  weighted residual variance 8.9e-16 (varying) and 8.9e-16 (equal).
  Poisson 3.6e-10, binomial counts 4.7e-15, shared-term block vs stacked WLS
  3.4e-10.
- Ordinal, one class, vs `MASS::polr`: log likelihood +7.2e-10 (latents
  higher); hand cumulative-logit likelihood at the estimates 6e-14 (c18.R,
  two classes, shared slope, offset); category probabilities match
  `predict(polr, type = "probs")` to 4e-7 (polr's own optimum).
- Random intercept, one class, vs `nlme::lme(method = "ML")`: 1.4e-13.
- predict: `class_response` = X beta (+ offset) exactly (0); priors = softmax
  of W gamma (1.1e-16); `response` = rowSums(prior * mean) (4.4e-16); ordinal
  probability table and expected scores vs hand (1.7e-16, 4.4e-16).
- Standard errors vs inverse numDeriv Hessian (max relative difference):
  Gaussian + membership covariate 1.6e-9, equal variance + common 1.0e-9,
  group 1.4e-9, two-level 3.4e-6, Poisson 1.0e-7, NB 6.2e-6, binomial counts
  2.4e-9, ordinal group 1.3e-8, ordinal shared slope + offset 7.5e-8,
  no-intercept 2.0e-9, two-level with both covariates 1.3e-8. Robust vs hand
  sandwich: group 2.5e-8, observation clustered on id 4.5e-9, weighted 4.3e-9.
- The growth bootstrap table (`parameter_inference(<growth>, method =
  "bootstrap")`) still carries `replicate_identity` (block, class, term,
  estimate); plain mixture regressions are refused by
  `parameter_inference()` with `latents_bad_argument` (by design).

## 1. Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.mixture_log1pexp` | R/mixture-regression-engine.R:34 | ran synthetic (binomial M-step vs glm, c1.R) | ok |
| `.mixture_linear_predictors` | R/mixture-regression-engine.R:41 | ran synthetic (predict vs X %*% beta, c5.R; shifted-design invariance c8.R) | ok |
| `.mixture_log_density` | R/mixture-regression-engine.R:76 | ran synthetic (LL by hand: ordinal c2/c18, RI c3; numDeriv SEs c4) | ok |
| `.mixture_negative_binomial_log_density` | R/mixture-regression-engine.R:108 | ran synthetic (expansion vs exact log1p sum at alpha 1e-4, 1e-6: 4e-16, 1e-15) | ok |
| `.mixture_expectation` | R/mixture-regression-engine.R:138 | ran synthetic (zero-weight invariance c6.R, all tables diff 0) | ok |
| `.mixture_unweighted_expectation` | R/mixture-regression-engine.R:147 | read; covered by c6.R | ok |
| `.mixture_family_derivatives` | R/mixture-regression-engine.R:161 | ran synthetic (Poisson/binomial M-step vs glm(weights) c1.R: <4e-10) | ok |
| `.mixture_block_diagonal` | R/mixture-regression-engine.R:189 | read; exercised by every Newton M-step (c1.R) | ok |
| `.mixture_regression_objective` | R/mixture-regression-engine.R:206 | read; exercised by c1.R | ok |
| `.mixture_update_regression` | R/mixture-regression-engine.R:226 | ran synthetic (vs lm/glm(weights = tau), common term vs stacked WLS, c1.R; fixed point test) | ok (ridge logged) |
| `.mixture_update_variance` | R/mixture-regression-engine.R:385 | ran synthetic (varying and equal vs weighted RSS, c1.R: <2e-15) | ok |
| `.mixture_update_dispersion` | R/mixture-regression-engine.R:415 | ran synthetic (vs optimize() of the weighted NB objective, c1.R: equal to 7 digits) | ok |
| `.mixture_maximization` | R/mixture-regression-engine.R:483 | read; covered by every fit | ok |
| `.mixture_em_settings` | R/mixture-regression-engine.R:497 | read | ok |
| `.mixture_observed` | R/mixture-regression-engine.R:507 | read; zero-weight invariance c6.R | ok |
| `.mixture_em` | R/mixture-regression-engine.R:520 | ran synthetic (max_iter 0/3, n_starts 0, c7.R) | ok (class-size rule logged) |
| `enumerate_regressions` | R/mixture-regression-enumerate.R:123 | ran synthetic (grid 1:3, failing cell, two-level grid, BLRT, c9.R) | logged (BIC n across nestings) |
| `.mixture_try_fit` | R/mixture-regression-enumerate.R:208 | ran synthetic (failing cell becomes note, c9.R) | ok |
| `.mixture_enumeration_row` | R/mixture-regression-enumerate.R:222 | ran synthetic (c9.R) | ok |
| `.mixture_blrt` | R/mixture-regression-enumerate.R:244 | ran synthetic (Gaussian and NB grids with bootstrap, c9.R); test-audit-1101-mixture-regression.R | fixed |
| `as.data.frame.latents_regression_enumeration` | R/mixture-regression-enumerate.R:346 | ran synthetic (c9.R) | ok |
| `get_results.latents_regression_enumeration` | R/mixture-regression-enumerate.R:353 | ran synthetic (fit/model/tree/all, missing model -> latents_bad_argument, c9.R) | ok |
| `summary.latents_regression_enumeration` | R/mixture-regression-enumerate.R:391 | ran synthetic (c9.R) | ok |
| `plot.latents_regression_enumeration` | R/mixture-regression-enumerate.R:401 | ran synthetic (pdf(NULL), c9.R) | ok |
| `print.latents_regression_enumeration` | R/mixture-regression-enumerate.R:416 | ran synthetic (c9.R) | ok |
| `.mixture_pack` | R/mixture-regression-inference.R:21 | ran synthetic (ordinal with all slopes common / o ~ 1: crashed); test-audit-1101-mixture-regression.R | fixed |
| `.mixture_unpack` | R/mixture-regression-inference.R:49 | ran synthetic (pack/unpack round trip inside numDeriv checks, c4.R) | ok |
| `.mixture_unit_scores` | R/mixture-regression-inference.R:87 | ran synthetic (robust SEs vs hand sandwich of numDeriv unit gradients: group 2.5e-8, obs+id 4.5e-9, weighted 4.3e-9; c4b/c16.R) | ok |
| `.mixture_total_score` | R/mixture-regression-inference.R:133 | ran synthetic (score_norm equals numDeriv gradient, c4.R) | ok |
| `.mixture_observed_information` | R/mixture-regression-inference.R:154 | ran synthetic (SE vs inverse numDeriv Hessian, 9 variants, c4/c4b/c11/c18.R) | ok |
| `.mixture_chart` | R/mixture-regression-inference.R:182 | ran synthetic (SEs invariant to shifting x by 50, u by -3000 x100, v by 40; ordinal and Gaussian, c8.R: <2e-13) | ok |
| `.mixture_inference` | R/mixture-regression-inference.R:288 | ran synthetic (observed/robust/opg, weighted, c4/c7/c16.R) | ok |
| `.mixture_boundary_coordinates` | R/mixture-regression-inference.R:351 | ran synthetic (NB boundary fit, c9.R; ordinal gap boundary, c4.R) | ok |
| `.mixture_invert` | R/mixture-regression-inference.R:374 | ran synthetic (singular two-level fit warns latents_singular_information, c5/c9.R) | ok |
| `.mixture_delta` | R/mixture-regression-inference.R:401 | ran synthetic (sigma/dispersion/prior rows; ordinal band SE vs numDeriv 4.6e-10, c17.R) | ok |
| `.mixture_tables` | R/mixture-regression-methods.R:5 | read | ok |
| `.mixture_resolve_inference` | R/mixture-regression-methods.R:13 | ran synthetic (stored, recomputed, weighted default, c7/c15.R) | ok |
| `.mixture_class_means` | R/mixture-regression-methods.R:31 | ran synthetic (ordinal expected score vs hand, c10.R; binomial/Poisson c5.R) | ok |
| `.mixture_coefficient_table` | R/mixture-regression-methods.R:43 | ran synthetic (vs polr c2.R; ordinal common c18.R) | ok |
| `.mixture_class_table` | R/mixture-regression-methods.R:80 | ran synthetic (zero-weight invariance c6.R) | ok |
| `.mixture_membership_table` | R/mixture-regression-methods.R:142 | ran synthetic (1-class empty table lacks `reference`, c11.R) | logged (layout frozen) |
| `.mixture_group_class_table` | R/mixture-regression-methods.R:170 | ran synthetic (single-level empty table lacks `n_groups_assigned`; SE delta order checked, c11.R) | logged (layout frozen) |
| `.mixture_relative_entropy` | R/mixture-regression-methods.R:212 | ran synthetic (zero-weight invariance c6.R) | ok |
| `.mixture_complete_entropy` | R/mixture-regression-methods.R:220 | ran synthetic (zero-weight invariance incl. two-level, c6.R) | ok |
| `.mixture_fit_table` | R/mixture-regression-methods.R:239 | ran synthetic (BIC/AIC vs stats::BIC/AIC, c7.R; c6.R) | ok (cross-nesting BIC logged) |
| `.mixture_assignment_table` | R/mixture-regression-methods.R:278 | ran synthetic (newdata posterior = fitted posterior, c5.R) | ok |
| `.mixture_group_table` | R/mixture-regression-methods.R:298 | ran synthetic (c10.R) | ok |
| `.mixture_fitted_table` | R/mixture-regression-methods.R:319 | ran synthetic (ordinal codes and expected scores, c10.R) | ok |
| `.mixture_classification_table` | R/mixture-regression-methods.R:337 | ran synthetic (zero-weight invariance c6.R) | ok |
| `.mixture_recovery_table` | R/mixture-regression-methods.R:372 | ran synthetic (row and group level, varying truth refused, c7/c10.R) | ok (weights logged) |
| `.mixture_table` | R/mixture-regression-methods.R:430 | ran synthetic (every table via "all", c7.R) | ok |
| `.mixture_present` | R/mixture-regression-methods.R:456 | ran synthetic (printed tables, c7.R) | ok |
| `get_results.latents_mixture_regression` | R/mixture-regression-methods.R:529 | ran synthetic (c7.R) | ok |
| `as.data.frame.latents_mixture_regression` | R/mixture-regression-methods.R:551 | ran synthetic (what =, extra argument refused, c15.R) | ok |
| `.mixture_require_level` | R/mixture-regression-methods.R:569 | ran synthetic (level = 2 -> latents_bad_argument, c7.R) | ok |
| `print.latents_mixture_regression` | R/mixture-regression-methods.R:577 | ran synthetic (with and without inference, c7.R) | ok |
| `summary.latents_mixture_regression` | R/mixture-regression-methods.R:614 | ran synthetic (c7/c15.R) | ok |
| `print.summary_latents_mixture_regression` | R/mixture-regression-methods.R:626 | ran synthetic (c7.R) | ok |
| `get_results.summary_latents_mixture_regression` | R/mixture-regression-methods.R:638 | ran synthetic (level refused, "all", c15.R) | ok |
| `as.data.frame.summary_latents_mixture_regression` | R/mixture-regression-methods.R:645 | ran synthetic (c15.R) | ok |
| `.mixture_summary_results` | R/mixture-regression-methods.R:661 | ran synthetic (c15.R) | ok |
| `coef.latents_mixture_regression` | R/mixture-regression-methods.R:675 | ran synthetic (c18.R) | ok |
| `vcov.latents_mixture_regression` | R/mixture-regression-methods.R:682 | ran synthetic (type =, vcov_type in ... warned, bad type refused, c15.R) | ok (message logged) |
| `confint.latents_mixture_regression` | R/mixture-regression-methods.R:693 | ran synthetic (positions, names, level, bad parm, c7.R) | ok |
| `.mixture_confint_rows` | R/mixture-regression-methods.R:716 | ran synthetic (c7.R) | ok |
| `logLik.latents_mixture_regression` | R/mixture-regression-methods.R:730 | ran synthetic (AIC()/BIC() agree with fit table, c7.R) | ok |
| `.mixture_weighted_sample` | R/mixture-regression-methods.R:747 | ran synthetic (zero-weight invariance c6.R) | ok |
| `nobs.latents_mixture_regression` | R/mixture-regression-methods.R:760 | ran synthetic (c7.R) | ok |
| `.mixture_new_spec` | R/mixture-regression-predict.R:16 | ran synthetic (newdata membership, group, two-level, missing id, c5/c15.R) | ok |
| `.mixture_prediction_regression` | R/mixture-regression-predict.R:68 | ran synthetic (class_response = X %*% beta + offset, c5.R) | ok |
| `.mixture_row_priors` | R/mixture-regression-predict.R:107 | ran synthetic (priors = softmax(W gamma), group and two-level, c5.R) | ok |
| `predict.latents_mixture_regression` | R/mixture-regression-predict.R:124 | ran synthetic (response/class_response/posterior/probabilities vs hand and polr, c2/c5/c10/c15.R) | ok |
| `.mixture_prediction_data` | R/mixture-regression-predict.R:171 | ran synthetic (missing column/NA -> latents_bad_data, c15.R) | ok |
| `.mixture_prediction_trials` | R/mixture-regression-predict.R:190 | ran synthetic (counts without trials refused, explicit trials, posterior mismatch refused, c15.R) | ok |
| `simulate.latents_mixture_regression` | R/mixture-regression-predict.R:229 | ran synthetic (seed restores caller stream, nsim = 0 refused, c7/c15.R) | ok |
| `.mixture_draw` | R/mixture-regression-predict.R:239 | ran synthetic (BLRT draws for Gaussian/NB, c9.R) | ok |
| `plot.latents_mixture_regression` | R/mixture-regression-predict.R:274 | ran synthetic (fitted/coefficients/posteriors, trajectory views for group fits, refused view, pdf(NULL), c7/c10.R) | ok |
| `.mixture_trajectory_views` | R/mixture-regression-predict.R:316 | read | ok |
| `.mixture_plot_views` | R/mixture-regression-predict.R:327 | read; c10.R | ok |
| `.mixture_plot_fitted` | R/mixture-regression-predict.R:333 | ran synthetic (c7/c10.R) | ok |
| `.mixture_model_data` | R/mixture-regression-predict.R:379 | read | ok |
| `.mixture_typical_rows` | R/mixture-regression-predict.R:387 | ran synthetic (through plot and trajectory band, c10/c17.R) | ok |
| `.mixture_plot_coefficients` | R/mixture-regression-predict.R:406 | ran synthetic (incl. ordinal thresholds, c10.R) | ok |
| `.mixture_plot_posteriors` | R/mixture-regression-predict.R:430 | ran synthetic (c7/c10.R) | ok |
| `.mixture_class_scales` | R/mixture-regression-predict.R:462 | read (Okabe-Ito, shapes and line types pair colour) | ok |
| `mixture_regression` | R/mixture-regression.R:327 | ran synthetic (every family and nesting, validation paths, c4-c19.R) | ok |
| `.mixture_as_formula` | R/mixture-regression.R:508 | ran synthetic (names, "intercept", formula, c7.R) | ok (logged) |
| `.mixture_weight_spec` | R/mixture-regression.R:540 | ran synthetic (NA/char/absent weights, varying within id; inline probe) | ok (logged) |
| `.mixture_gaussian_reference` | R/mixture-regression.R:573 | read; covered by test-review-1026-mixreg.R | ok |
| `.mixture_is_count` | R/mixture-regression.R:593 | ran synthetic (n_classes 0 refused, c7.R) | ok |
| `.mixture_with_seed` | R/mixture-regression.R:606 | ran synthetic (simulate seed restores stream, c15.R) | ok |
| `.mixture_spec` | R/mixture-regression.R:625 | ran synthetic (refusals: group w/o id, common, membership varying in group, too many classes, binary single level, c7.R) | ok (locale order logged) |
| `.mixture_row_support` | R/mixture-regression.R:852 | ran synthetic (zero weights, c6.R) | ok |
| `.mixture_outside_support` | R/mixture-regression.R:873 | read; covered by test-review-1020-zero-weights.R | ok |
| `.mixture_categorical_columns` | R/mixture-regression.R:896 | read; covered by .mixture_require_levels tests | ok |
| `.mixture_require_levels` | R/mixture-regression.R:925 | read; covered by test-review-1026-mixreg.R | ok |
| `.mixture_common_columns` | R/mixture-regression.R:961 | ran synthetic (common intercept, absent term, factor, c7/c5.R) | ok |
| `.mixture_membership_design` | R/mixture-regression.R:995 | ran synthetic (character membership covariate, c7.R) | ok |
| `.mixture_response` | R/mixture-regression.R:1015 | ran synthetic (cbind, factor, logical, counts, c1/c7.R) | ok |
| `.mixture_group_design` | R/mixture-regression.R:1066 | ran synthetic (varying membership in group refused, c7.R) | ok |
| `.mixture_carrying_spec` | R/mixture-regression-starts.R:17 | ran synthetic (zero-weight invariance across nestings, c6.R) | ok |
| `.mixture_null_params` | R/mixture-regression-starts.R:55 | ran synthetic (c1.R) | ok |
| `.mixture_soften` | R/mixture-regression-starts.R:82 | read | ok |
| `.mixture_residual_labels` | R/mixture-regression-starts.R:98 | ran synthetic (n_starts = 0 uses it alone, c7.R) | ok |
| `.mixture_seed_expectation` | R/mixture-regression-starts.R:129 | read; exercised by every start (c4/c11.R) | ok |
| `.mixture_run_starts` | R/mixture-regression-starts.R:160 | ran synthetic (starts tables, degenerate finalists census, c19/c20.R) | logged (degenerate finalists) |
| `.mixture_select` | R/mixture-regression-starts.R:219 | ran synthetic (unconverged/degenerate/boundary warnings, c7/c9/c20.R) | ok |
| `.mixture_count_parameters` | R/mixture-regression-starts.R:311 | ran synthetic (n_parameters vs packed length in every variant, c4.R) | ok (1 x H cells logged) |
| `.mixture_order_classes` | R/mixture-regression-starts.R:340 | ran synthetic (c11.R) | ok |
| `.mixture_class_orders` | R/mixture-regression-starts.R:355 | ran synthetic (c11.R) | ok |
| `.mixture_permute` | R/mixture-regression-starts.R:386 | ran synthetic (two-level with covariates permuted: LL, tau, rho unchanged, c11.R) | ok |
| `.ordinal_response` | R/mixture-regression-ordinal.R:50 | ran synthetic (ordered factor, integer codes, c2/c10.R) | ok |
| `.ordinal_check_identified` | R/mixture-regression-ordinal.R:107 | read; covered by test-ordinal-regression.R | ok |
| `.ordinal_check_patterns` | R/mixture-regression-ordinal.R:143 | read; covered by test-ordinal-regression.R | ok (zero-weight patterns logged) |
| `.ordinal_threshold_labels` | R/mixture-regression-ordinal.R:167 | ran synthetic (c2.R) | ok |
| `.ordinal_start_thresholds` | R/mixture-regression-ordinal.R:179 | read | ok |
| `.ordinal_log1mexp` | R/mixture-regression-ordinal.R:196 | read (Maechler 2012 split at log 2) | ok |
| `.ordinal_bounds` | R/mixture-regression-ordinal.R:214 | ran synthetic (via LL vs hand formula, c18.R: 6e-14) | ok |
| `.ordinal_log_density` | R/mixture-regression-ordinal.R:227 | ran synthetic (LL vs polr and hand formula, c2/c18.R) | ok |
| `.ordinal_derivatives` | R/mixture-regression-ordinal.R:246 | ran synthetic (SE vs numDeriv Hessian, ordinal group 1.3e-8, common+offset 7e-8) | ok |
| `.ordinal_indicators` | R/mixture-regression-ordinal.R:268 | read | ok |
| `.ordinal_ordered` | R/mixture-regression-ordinal.R:276 | read | ok |
| `.ordinal_update` | R/mixture-regression-ordinal.R:296 | ran synthetic (1-class vs MASS::polr LL +7e-10, c2.R) | ok |
| `.ordinal_min_log_gap` | R/mixture-regression-ordinal.R:402 | ran synthetic (boundary warning, c4.R) | ok |
| `.ordinal_threshold_jacobian` | R/mixture-regression-ordinal.R:411 | ran synthetic (threshold SEs vs polr, c2.R) | ok |
| `.ordinal_pack` | R/mixture-regression-ordinal.R:423 | ran synthetic (c2/c4b.R) | ok |
| `.ordinal_natural` | R/mixture-regression-ordinal.R:439 | ran synthetic (c2/c4b.R) | ok |
| `.ordinal_threshold_scores` | R/mixture-regression-ordinal.R:458 | ran synthetic (robust ordinal SE vs hand sandwich 2.5e-8, c4b.R) | ok |
| `.ordinal_probabilities` | R/mixture-regression-ordinal.R:472 | ran synthetic (predict probabilities vs polr, vs hand for 2 classes, c2/c10.R) | ok |
| `.ordinal_expected_scores` | R/mixture-regression-ordinal.R:492 | ran synthetic (class_response vs hand, c10.R) | ok |
| `.ordinal_categories` | R/mixture-regression-ordinal.R:508 | read; via simulate/BLRT draw types | ok |
| `.ordinal_draw` | R/mixture-regression-ordinal.R:520 | read; exercised by simulate on ordinal fits | ok |
| `.ordinal_probability_table` | R/mixture-regression-ordinal.R:535 | ran synthetic (row/class/category order vs hand, c10.R: 1.7e-16) | ok |
| `.ordinal_threshold_table` | R/mixture-regression-ordinal.R:562 | ran synthetic (vs polr SEs, c2.R) | ok |
| `.ordinal_trajectory_band` | R/mixture-regression-ordinal.R:598 | ran synthetic (estimate exact, SE vs numDeriv delta 4.6e-10, c17.R) | ok |

## 2. Fixed

**F1. Ordinal mixture with no class-specific slope crashed (unclassed error).**
`mixture_regression(oc ~ x + offset(e), data, 2, family = "ordinal", id = "id",
class_level = "group", common = "x")`, and equally `o ~ 1` with
`class_level = "group"` (both documented: `common` shares slopes, thresholds
always vary; a predictor-free ordinal mixture is identified through repeated
outcomes), stopped once inference ran:
`Error in names(object) <- nm: 'names' attribute [1] must be the same length
as the vector [0]` (a `simpleError`, not a latents condition; also for one
class). Cause: `.mixture_pack()` named the empty class-slope block with
`paste0("coefficient.", character(0), ".", NULL)`, which returns one string.
Change: R/mixture-regression-inference.R:25, an empty block when `ncol(spec$x)
== 0`. After the fix the fit's likelihood equals a hand cumulative-logit
mixture likelihood (5.7e-14) and its standard errors the inverse numDeriv
Hessian (7.5e-8). Guarded by test "an ordinal mixture with every slope shared
(or none) is fitted and inferred"; with the fix reverted it errors at
test-audit-1101-mixture-regression.R:47.

**F2. Bootstrap LRT replicates repeated boundary warnings.** The code's own
contract (comment at R/mixture-regression-enumerate.R:280) is that replicate
refits' routine warnings are not repeated on the console; a
negative-binomial grid with `bootstrap = 3` printed six `latents_boundary`
warnings from the replicate refits on top of the observed fits' own (c9.R).
With `bootstrap = 199` that is hundreds of identical warnings. Change:
R/mixture-regression-enumerate.R:301-302 muffle `latents_boundary` and
`latents_failed_starts` inside replicate refits only; such a replicate stays
valid (a Poisson-limit dispersion or a zero threshold gap is a legitimate
maximum, and a start that could not be built does not invalidate the refit);
the observed fits still warn. No number changes. Guarded by test "bootstrap
replicates do not repeat boundary warnings" (boundary warnings with
`bootstrap = 3` equal those with `bootstrap = 0`); fails (boundary count 6 vs observed-only) with the fix reverted.

## 3. Logged, not fixed

- **Empty membership / group-class tables have fewer columns (layout frozen;
  fix made and reverted on the coordinator's instruction).** One-class
  `get_results(fit, "membership")`: `model class term estimate std_error
  statistic p_value conf_low conf_high odds_ratio`; with classes:
  `model class term reference estimate ...`. Single-level `"group_classes"`:
  `group_class class probability group_share`; two-level: `group_class class
  probability group_share n_groups_assigned`. So tables of a class-count grid
  cannot be stacked with `rbind()`. Adding the columns changed the
  engine-golden fingerprint `group_classes@display` (length 4 vs 5) of all 21
  `reg_*` cases and `growth_lcga`, so both were reverted
  (R/mixture-regression-methods.R unchanged).

- **BIC across nestings in a crossed regression grid (owner decision; frozen
  criterion).** `enumerate_regressions(..., id = "id", n_group_classes = 1:2)`:
  cells with one group class are observation-nesting fits whose `bic` uses the
  number of rows; cells with two or more use the number of groups. On 120 rows
  in 24 groups with no group structure, the 2 x 2 cell has the same log
  likelihood as 2 x 1 (-159.3889) and two more parameters, yet BIC 347.38 vs
  352.29, and `best_bic` marks 2 x 2 (c9.R); likewise 1 x 2 (498.39) beats
  1 x 1 (500.04). `bic_rows` orders them correctly. multilpa uses groups for
  every cell and the growth mixture with `cluster` uses clusters at
  `n_group_classes = 1` ("comparable by BIC"). Fixing changes a reported
  number (the H = 1 `bic`) or the grid's marking, so it is left to the owner.
- **One class with several group classes is accepted (owner decision).**
  `n_classes = 1, n_group_classes = 2` fits a model whose group-class logits do
  not enter the likelihood, yet counts them (4 parameters vs 3) and gives a
  lower BIC through the group-count penalty (above). multilpa refuses the same
  cell (`latents_unidentified`, an NA row in its grid). Refusing it here would
  change the regression grid's rows, and the BLRT of 2 vs 1 class at H = 2
  uses that cell as its null.
- **Degenerate starts take finalist slots (estimation algorithm, frozen).**
  `.mixture_run_starts()` carries the best half of the screened starts by
  likelihood, degenerate ones included, but degenerate starts are never
  continued nor selectable. With collapsing Gaussian classes (high
  likelihood) usable starts are left at 50 screening iterations: K = 6, seed
  4, 4 of 6 finalist slots degenerate, 2 usable starts not continued (c20.R);
  on `y ~ 1` with pure noise the selected start is a screened one and warns
  "did not converge in 50 iterations ... Raise `max_iter`", which cannot help.
  Ranking finalists among non-degenerate starts changes no fit without
  degenerate starts, but changes which maximum is reported where they occur.
- The Gaussian Newton M-step adds a ridge `1e-10 (1 + max diag)`, so one step
  from zero differs from exact weighted least squares by ~1e-9 (relative);
  the WLS solution is its exact fixed point, so converged estimates are
  unaffected.
- `.mixture_em()` counts a degenerate class as one with fewer expected rows
  than `ncol(x) + 1`; the comment says "slopes and an intercept" but `x`
  already holds the intercept (the + 1 is the variance for Gaussian, an extra
  row otherwise). Conservative, starts only.
- `common = c("intercept", "x")` silently keeps the intercept class-specific
  (only `"intercept"` alone is refused, as documented).
- `group_levels <- sort(unique(id))` sorts character ids by the session's
  collation, so the "groups" row order and the random labels' unit order can
  differ between locales (C vs en_US) for mixed-case ids.
- A missing weight on a row already dropped by `missing = "omit"` is still
  refused (`latents_unsupported_weights`, "drop those rows first").
- `.ordinal_check_patterns()` counts covariate patterns over zero-weight rows
  too, so a zero-weight row can lift the identification capacity.
- Recovery tables count zero-weight units (as multilpa's do).
- `vcov(fit, type = "xx")` names `vcov_type` in its message (class correct).
- BLRT replicates that fail with an unexpected (non-data) error are counted as
  invalid rather than raised; they are warned about via
  `latents_failed_replicates`.
- A single-level ordinal fit on 300 rows (c4.R) reached a solution with two
  zero threshold gaps and warned `latents_boundary` as designed.

## 4. Needs another slice

- R/weights.R:87-93 (`.latents_weighted_vcov()`): the
  refusal of `vcov_type = "observed"` for a weighted mixture regression
  suggests `method = "bootstrap"`, an argument `mixture_regression()` does not
  have.

## 5. Test runs

All with `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`,
`devtools::load_all()`, `testthat::test_file()` (pass / fail / skip / error /
warning), after the three fixes (first run) and after reverting the
empty-table layout change (second run):

```
first run (fixes F1, F2 and the reverted schema change in place)
test-audit-1101-mixture-regression.R          30 0 0 0 0
test-mixture-regression.R                     73 0 0 0 0
test-mixture-regression-labels.R               7 0 0 0 0
test-mixture-regression-prediction-design.R   12 0 0 0 0
test-negative-binomial-regression.R           17 0 0 0 0
test-ordinal-regression.R                     38 0 0 0 0
test-release-audit-regression.R               58 0 0 0 0
test-release-regression.R                    100 0 0 0 0
test-review-1015-regression.R                 27 0 0 0 0
test-review-1026-mixreg-r2.R                  77 0 0 0 0
test-review-1026-mixreg.R                     20 0 0 0 0
test-review-regression.R                      50 0 0 0 0
test-api-compat-0921.R                        44 0 0 0 0
test-compat-0921.R                            10 0 0 0 0
test-grid-contract.R                         237 0 0 0 0
test-review-1101-bootstrap-tables.R           67 0 0 0 0

second run (final code: F1 and F2 only)
test-audit-1101-mixture-regression.R          26 0 0 0 0
test-mixture-regression.R                     73 0 0 0 0
test-api-compat-0921.R                        44 0 0 0 0
test-compat-0921.R                            10 0 0 0 0
test-grid-contract.R                         237 0 0 0 0
```

`git status tests/testthat/_snaps/` is clean after both runs
(tests/testthat/_snaps/mixture-regression.md unchanged).

Sabotage: with F1 and F2 reverted together, the new file reported an error
at the ordinal test (`names(object) <- nm`) and a failure at the bootstrap
warning test.

`Rscript equivalence/engine-golden/compare.R` (run once, with the
since-reverted schema change in place): DIFFERENT = the 13 known bootstrap
cases (`replicate_rows` -> `replicate_identity`) plus 22 cases differing only
at `probes/tables/value/group_classes@display` (length 4 vs 5), caused by
the empty group-class schema change, now reverted. No other difference;
every other reg_*, growth_* and ordinal case identical. Not rerun after the
revert (the revert restores the file byte for byte; `git diff` of
R/mixture-regression-methods.R is empty).

`devtools::document()` not needed: no roxygen changed (comments inside
function bodies only); NAMESPACE and man/ untouched by this slice.
