# Audit 2026-11-01 — slice `profile-fitting`

Files: R/multilpa.R, R/lpa.R, R/multilca.R, R/categorical.R,
R/extra-indicators.R, R/input-checks.R, R/utils.R (127 functions plus 5
constants / `%||%`). Scratch scripts behind the "ran synthetic" entries are in
the session scratchpad (`audit/num1.R`–`num8.R`, `probe1.R`, `probe2.R`,
`plumb1.R`).

Independent checks of the likelihood pieces (all on synthetic data, this session):

- Ordinal adjacent-category probabilities vs hand formula: max diff 1.1e-16, rows sum to 1.
- Ordinal Newton M-step vs `optim(BFGS)`: same objective (-47.11068), parameters within 3.2e-7.
- NB2 log-dispersion score and curvature vs 300-bit Rmpfr finite differences on a 120-case grid
  (y 0..400, mu 0.5..500, alpha 1e-7..5, covering the small-alpha series and the r >= 1000 asymptotic branch):
  max relative error 2.9e-13 (score), 5.2e-13 (curvature).
- NB Newton M-step vs `optim`: same objective (varying -904.7938, shared -909.7629), params within 1.5e-7;
  analytic Hessian vs `numDeriv::hessian` 2.7e-6 on a scale of 255.
- Poisson limit: binomial (under-dispersed) counts put both dispersions on the 1e-8 floor with means equal to the weighted means.
- Extra-block E-step density vs hand `dpois` + ordinal log P: 8.9e-16; centred/offset form reproduces the plain one to 1.8e-15.
- Scores (Fisher identity) vs `numDeriv::grad` of the observed-data log likelihood with NA cells:
  Poisson + ordinal 1.9e-7, NB 9.4e-8. Jacobian vs `numDeriv::jacobian` 3.4e-11.
- Count M-step with NA cells = hand posterior-weighted means (2.7e-15).
- Integer sampling weights vs repeated rows (Gaussian + categorical + ordinal + 2 counts, varying and equal variances):
  estimates agree to <= 8.3e-7, scaled log likelihoods to 1.9e-11.
- Parameter counts vs hand: 21, 19 (Poisson, varying/equal), 17 (NB varying), 15 (NB equal).
- Joined-pair MM M-step (`.multilpa_partial_maximize`) vs `optim` on the concave objective: 6e-8.
- 20000 draws from `.latents_draw_extra()`: ordinal shares and NB mean/variance match the model.

## 1. Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `lpa` | R/lpa.R:56 | ran synthetic (probe battery: id/n_group_classes refused, ordinal/count via ...); test-single-level.R | ok |
| `lca` | R/lpa.R:115 | ran synthetic (n_profiles refused; ordinal+count+NB fits; negative-df warning); test-single-level.R, audit test | ok |
| `enumerate_lpa` | R/lpa.R:175 | read; covered by test-single-level.R | ok |
| `enumerate_lca` | R/lpa.R:236 | read; shares .latents_lca_categorical (fix 2); test-single-level.R | ok |
| `.multilpa_single_level_call` | R/lpa.R:256 | ran synthetic (lca id=, lpa n_group_classes= -> latents_bad_argument) | ok |
| `.multilpa_quiet_single_level` | R/lpa.R:275 | read; covered by test-single-level.R | ok |
| `.multilpa_log_sum_exp` | R/multilpa.R:5 | ran synthetic (rows -1000/0 and -1001/log 3 -> 0, log 3) | ok |
| `.multilpa_expectation` | R/multilpa.R:30 | ran via weighted-vs-replicated fits (estimates agree 1e-7) | ok |
| `.multilpa_maximization` | R/multilpa.R:56 | ran via weighted-vs-replicated fits, varying and equal variances | logged (two unclassed internal stop()s) |
| `.multilpa_initialize` | R/multilpa.R:171 | read; ran via every fit (zero-weight row path in fit3) | ok |
| `.multilpa_group_class_seed` | R/multilpa.R:351 | read; covered by test-multilpa.R | ok |
| `.multilpa_initialization_center` | R/multilpa.R:374 | ran synthetic (zero-weight row ignored, NA skipped: 1, 52.5) | ok |
| `.multilpa_random_assignments` | R/multilpa.R:402 | read; covered by emEM tests in test-multilpa.R | ok |
| `.multilpa_ward_assignments` | R/multilpa.R:429 | read; covered by test-multilpa.R (start 1) | ok |
| `.multilpa_validate_start` | R/multilpa.R:470 | ran synthetic (labelled blocks reordered; bad sums -> latents_bad_start) | ok |
| `.multilpa_align_response_start` | R/multilpa.R:598 | ran synthetic (item and category labels reorder; unknown label -> latents_bad_start) | ok |
| `.multilpa_validate_response_start` | R/multilpa.R:645 | ran synthetic via validate_start | ok |
| `.multilpa_em` | R/multilpa.R:692 | read; covered by every fit | ok |
| `.multilpa_time_values` | R/multilpa.R:756 | ran synthetic (absent, NA, duplicate, T2/T10 -> classed) | ok |
| `.multilpa_single_level` | R/multilpa.R:810 | read; covered by test-single-level.R | logged (notice wording) |
| `multilpa` | R/multilpa.R:1513 | ran synthetic (75-case input-contract probe battery) | ok |
| `.multilpa_setup` | R/multilpa.R:1566 | ran synthetic (probe battery) | fixed (data label) |
| `.multilpa_warn_negative_df` | R/multilpa.R:2072 | ran synthetic (3 binary x 2 classes: none; x 3 classes and 2 ordinal x 3 classes: warns) | ok |
| `.latents_match_formal` | R/multilpa.R:2107 | ran synthetic (count_model = "nb", missing = "drop" -> latents_bad_argument) | ok |
| `.multilpa_check_columns` | R/multilpa.R:2123 | ran synthetic (dup vars, NA vars, id in vars, 1 row, absent id) | ok |
| `.multilpa_check_arguments` | R/multilpa.R:2160 | ran synthetic (n_profiles 0/1.5/"2"/NA, max_iter -1, tol 0, seed 1.5) | ok |
| `.multilpa_prepare_indicators` | R/multilpa.R:2218 | ran synthetic (character, Inf, NA, constant continuous) | ok |
| `.multilpa_contributing_rows` | R/multilpa.R:2290 | ran via zero-weight fit (LL equal to fit without the row, diff 4e-11) | ok |
| `.multilpa_prepare_groups` | R/multilpa.R:2311 | ran synthetic (factor ids first-occurrence order; NA/list ids refused) | ok |
| `.multilpa_label_parameters` | R/multilpa.R:2329 | read; covered by every fit's dimnames | ok |
| `.multilpa_count_parameters` | R/multilpa.R:2356 | ran synthetic (3 profiles, 2 classes, 4 full + cat (2,3) = 56 = hand) | ok |
| `.latents_profile_parameters` | R/multilpa.R:2394 | ran synthetic (VVI 4 -> 8, EEI -> 0) | ok |
| `.latents_warn_thin_profiles` | R/multilpa.R:2436 | read; covered by test-multilpa.R | ok |
| `.multilpa_covariate_model` | R/multilpa.R:2481 | read (pass-through) | ok |
| `.multilpa_start_kmeans` | R/multilpa.R:2525 | read (muffles two diagnosed kmeans warnings, recorded as failed starts) | ok |
| `.multilpa_categorical_levels` | R/categorical.R:17 | ran synthetic (logical, factor unused level, zero-weight-only value) | logged (NaN as missing) |
| `.multilpa_encode_categorical` | R/categorical.R:89 | ran synthetic (codes, levels, NA for excluded-only value) | ok |
| `.multilpa_attach_partial` | R/categorical.R:127 | read; covered by local-dependence tests | ok |
| `.multilpa_partial_of` | R/categorical.R:145 | read | ok |
| `.multilpa_categorical_observed` | R/categorical.R:159 | ran synthetic (c(1,NA),(2,NA) -> 2, 0) | ok |
| `.multilpa_partial_sets` | R/categorical.R:171 | ran synthetic via partial_maximize | ok |
| `.latents_warn_suspect_categorical` | R/categorical.R:189 | ran synthetic (continuous declared categorical -> latents_suspect_categorical) | ok |
| `.latents_report_once` | R/categorical.R:230 | read; covered by test-categorical-suspect.R | ok |
| `.multilpa_categorical_values` | R/categorical.R:252 | ran synthetic (typed values; latents_excluded attr for zero-weight value) | ok |
| `.multilpa_restore_excluded` | R/categorical.R:277 | read; covered by test-categorical.R | ok |
| `.multilpa_categorical_log_density` | R/categorical.R:294 | ran via weighted-vs-replicated mixed fit | ok |
| `.multilpa_categorical_maximize` | R/categorical.R:333 | ran via weighted-vs-replicated fit (response probs agree 2e-8) | ok |
| `.multilpa_partial_maximize` | R/categorical.R:388 | ran synthetic vs optim on the concave objective (max diff 6e-8) | ok |
| `.multilpa_bound_probabilities` | R/categorical.R:439 | ran synthetic (c(100,0,0) -> .998,.001,.001; zeros -> uniform) | ok |
| `.multilpa_categorical_initialize` | R/categorical.R:469 | read; covered by every categorical fit | ok |
| `.multilpa_categorical_design` | R/categorical.R:500 | ran synthetic (NA cell takes observed shares) | ok |
| `.multilpa_categorical_parameters` | R/categorical.R:520 | ran via n_parameters hand counts | ok |
| `.multilpa_categorical_thresholds` | R/categorical.R:534 | ran synthetic (= qlogis of cumulative sums) | ok |
| `multilca` | R/multilca.R:55 | ran synthetic (missing id, NA vars, matrix data) | ok |
| `.latents_lca_categorical` | R/multilca.R:93 | ran synthetic (count_dispersion value read as count names) | fixed |
| `.latents_model_classes` | R/input-checks.R:11 | read (constant) | ok |
| `.latents_enumeration_classes` | R/input-checks.R:15 | read (constant) | ok |
| `.latents_profile_fit_text` | R/input-checks.R:23 | read (constant) | ok |
| `.latents_require` | R/input-checks.R:37 | ran via probes (classed, sequential) | ok |
| `.latents_match_arg` | R/input-checks.R:67 | ran via probes | ok |
| `.latents_describe` | R/input-checks.R:84 | ran via probes (matrix, list) | ok |
| `.latents_argument_label` | R/input-checks.R:125 | ran via probes | ok |
| `.latents_name_list` | R/input-checks.R:135 | ran via probes | ok |
| `.latents_did_you_mean` | R/input-checks.R:149 | ran via probes (time = "zz" -> did you mean g) | ok |
| `.latents_check_data_frame` | R/input-checks.R:171 | ran via probes | ok |
| `.latents_require_fit` | R/input-checks.R:200 | read; covered by test-input-like checks in other verbs' tests | ok |
| `.latents_refuse_missing_values` | R/input-checks.R:244 | ran via probes (Inf, NA, error vs fiml) | ok |
| `.latents_str` | R/input-checks.R:285 | read; covered by test-api-compat-0921.R / str tests | ok |
| `str.multilpa` | R/input-checks.R:296 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.multilpa_covariates` | R/input-checks.R:298 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.multilpa_transitions` | R/input-checks.R:300 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.multilpa_lta` | R/input-checks.R:302 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.multilpa_additive` | R/input-checks.R:304 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.multilpa_cross_level` | R/input-checks.R:306 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.latents_growth_mixture` | R/input-checks.R:308 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.latents_mixture_regression` | R/input-checks.R:310 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.multilpa_enumeration` | R/input-checks.R:312 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.latents_transition_enumeration` | R/input-checks.R:314 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.latents_family_enumeration` | R/input-checks.R:316 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.latents_regression_enumeration` | R/input-checks.R:318 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `str.multilpa_diagnostics` | R/input-checks.R:320 | read; one-line delegation to .latents_str (no Rd: project convention for S3 methods) | ok |
| `.latents_prepare_extra` | R/extra-indicators.R:50 | ran synthetic (factor, numeric, NA under fiml, zero-weight-only value, 1e16 values -> ordinal_values recorded) | ok |
| `.latents_ordinal_labels` | R/extra-indicators.R:175 | ran synthetic (1e16 and 1e16+2 get 17-digit labels) | ok |
| `.latents_extra_names` | R/extra-indicators.R:185 | read | ok |
| `.latents_ordinal_log_probabilities` | R/extra-indicators.R:196 | ran synthetic vs hand adjacent-category logit (1e-16; rows sum to 1); audit test | ok |
| `.latents_extra_log_density` | R/extra-indicators.R:216 | ran synthetic (profile column = hand dpois + ordinal log P, 9e-16; offset form consistent 2e-15) | ok |
| `.latents_ordinal_counts` | R/extra-indicators.R:288 | ran via ordinal score / M-step checks | ok |
| `.latents_ordinal_maximize` | R/extra-indicators.R:308 | ran synthetic vs optim (same maximum, params 3e-7); audit test | ok |
| `.latents_extra_maximize` | R/extra-indicators.R:386 | ran synthetic (count means with NA = hand weighted means, 3e-15) | ok |
| `.latents_extra_initialize` | R/extra-indicators.R:440 | read; covered by every ordinal/count fit | ok |
| `.latents_extra_clustering` | R/extra-indicators.R:458 | read (NA mean imputation for clustering only) | ok |
| `.latents_extra_n_parameters` | R/extra-indicators.R:473 | ran (Poisson, NB varying, NB equal = hand counts 21/19/17/15) | ok |
| `.latents_extra_matrix` | R/extra-indicators.R:485 | read | ok |
| `.latents_check_extra_arguments` | R/extra-indicators.R:498 | ran synthetic (ordinal = 2 / NA -> unclassed simpleError) | fixed |
| `.latents_label_extra` | R/extra-indicators.R:527 | ran via fits' dimnames | ok |
| `.latents_measurement_model` | R/extra-indicators.R:553 | ran synthetic | ok |
| `.latents_extra_coordinates` | R/extra-indicators.R:574 | ran synthetic (names align with .latents_extra_labels row by row) | ok |
| `.latents_extra_width` | R/extra-indicators.R:620 | ran (= length of coordinates) | ok |
| `.latents_extra_decode` | R/extra-indicators.R:628 | ran synthetic (round trip in score check) | ok |
| `.latents_extra_scores` | R/extra-indicators.R:674 | ran synthetic vs numDeriv gradient of observed LL (Poisson+ordinal 2e-7, NB 9e-8) | fixed (na.rm removed) |
| `.latents_extra_jacobian` | R/extra-indicators.R:734 | ran synthetic vs numDeriv jacobian (3e-11) | ok |
| `.latents_ordinal_frame` | R/extra-indicators.R:750 | ran (get_results(fit, "ordinal"), order and SEs) | ok |
| `.latents_count_frame` | R/extra-indicators.R:791 | ran (get_results(fit, "count_means"), shared dispersion SE) | ok |
| `.latents_extra_catalogue` | R/extra-indicators.R:819 | ran via get_results | ok |
| `.latents_draw_extra` | R/extra-indicators.R:838 | ran synthetic (20000 draws: ordinal shares and NB mean/variance match model) | ok |
| `.latents_prepare_extra_like` | R/extra-indicators.R:873 | ran via predict() on 1e16 ordinal values | ok |
| `.latents_print_extra` | R/extra-indicators.R:925 | ran (print(fit)) | ok |
| `.latents_draw_extra_columns` | R/extra-indicators.R:947 | ran via get_results(fit, "data") | logged (needs-other-slice) |
| `.latents_extra_labels` | R/extra-indicators.R:969 | ran synthetic (aligned with coordinates) | ok |
| `.latents_extra_rows` | R/extra-indicators.R:1009 | read | ok |
| `.lta_extra_table` | R/extra-indicators.R:1028 | ran (lta() ordinal/count tables) | ok |
| `.latents_as_column_type` | R/extra-indicators.R:1049 | ran synthetic | ok |
| `.latents_negative_binomial` | R/extra-indicators.R:1068 | read | ok |
| `.latents_dispersion_rows` | R/extra-indicators.R:1075 | ran via n_parameters | ok |
| `.latents_set_count_model` | R/extra-indicators.R:1086 | ran (count_dispersion = "equal" reaches the engine: one dispersion) | logged (ignored under Poisson) |
| `.latents_min_dispersion` | R/extra-indicators.R:1096 | read (constant) | ok |
| `.latents_negative_binomial_scores` | R/extra-indicators.R:1110 | ran synthetic vs 300-bit Rmpfr finite differences (rel err 3e-13) | ok |
| `.latents_negative_binomial_dispersion_derivatives` | R/extra-indicators.R:1134 | ran synthetic, 120-case grid incl. both expansion regimes, vs Rmpfr (score 3e-13, curvature 5e-13) | ok |
| `.latents_moment_dispersion` | R/extra-indicators.R:1203 | read | ok |
| `.latents_negative_binomial_maximize` | R/extra-indicators.R:1226 | ran synthetic vs optim, varying and shared (same maximum, params 1e-7) | ok |
| `.latents_poisson_limit` | R/extra-indicators.R:1323 | ran synthetic (binomial counts -> floor, means = weighted means) | ok |
| `.latents_warn_poisson_limit` | R/extra-indicators.R:1351 | read; covered by test-extra-indicators.R | ok |
| `.latents_negative_binomial_hessian` | R/extra-indicators.R:1380 | ran synthetic vs numDeriv hessian (3e-6 on scale 255) | ok |
| ``%||%`` | R/utils.R:7 | read | ok |
| `.multilpa_any_fit` | R/utils.R:22 | read | ok |
| `.multilpa_is_seed` | R/utils.R:36 | ran synthetic (1,-1 TRUE; 1.5, NA, Inf, 2^31 FALSE) | ok |
| `.multilpa_check_seed` | R/utils.R:47 | ran via probes | ok |
| `.multilpa_select_start` | R/utils.R:85 | ran synthetic (tie window, converged rule) | ok |
| `.multilpa_row_max` | R/utils.R:133 | ran synthetic (NA propagates) | ok |
| `.multilpa_notice` | R/utils.R:148 | read; covered by test-single-level.R | ok |
| `.multilpa_single_level_notice` | R/utils.R:168 | read | ok |
| `.multilpa_missing_id` | R/utils.R:185 | ran via probes (multilca without id) | ok |

## 2. Fixed

**F1 — `ordinal`/`count` given as a non-character value raised an unclassed error**
(R/extra-indicators.R:498, `.latents_check_extra_arguments()`). Reproduction:
`lpa(d, c("x","o"), 2, ordinal = 2)` and `ordinal = NA_character_`, `count = 2`,
`lca(..., ordinal = NA_character_)` all gave
`simpleError | ``ordinal`` must be a character vector of indicator names`, while every
other argument contract of the fitting verbs is `latents_bad_argument`. Cause: a bare
`stopifnot()`. Change: `.latents_require()` (same messages, class `latents_bad_argument`),
@return of the @noRd block updated. Guard: test "a non-character ordinal or count argument
is a classed bad argument". Broken on purpose (back to `stopifnot`): test errors.

**F2 — `lca()`/`multilca()`/`enumerate_lca()` dropped an indicator from the categorical set
when `count_model`/`count_dispersion` was given without `count`** (R/multilca.R:93,
`.latents_lca_categorical()`). `extra$count` on `list(...)` matches partially, so
`count_dispersion = "equal"` was read as `count = "equal"`. Reproduction:
`lca(d, c("equal", "a"), 2, count_dispersion = "equal")` (character `equal`) failed with
"`equal` is not numeric ... A continuous indicator must be a numeric vector"; with a 0/1
numeric `equal` the indicator would silently have been fitted as Gaussian. Change:
`extra[["categorical"]]`, `extra[["ordinal"]]`, `extra[["count"]]` (exact). Guard: test
"lca() keeps an indicator categorical whatever count_* values are given" (asserts the
response blocks, no Gaussian columns, and log likelihood / n_parameters equal to the fit
without `count_dispersion`). Broken on purpose (back to `$`): 4 expectations fail.

**F3 — `multilpa()` named an internal symbol in its data refusal** (R/multilpa.R:1673,
`.multilpa_setup()`). Reproduction: `multilpa(as.matrix(d), c("x","y"), "g", 2)` →
"`data` must be a data frame, but `.data` is a double matrix [120 x 4]". Cause:
`multilpa()` binds argument values to generated symbols (`.latents_bind_values`) before
calling the setup, so `substitute(data)` there saw `.data`. Change: the label comes from the
caller's own expression, `.call[["data"]]` (class unchanged, `latents_bad_data`; a symbol
`m` is still named `m`). Guard: test "multilpa() names the caller's data, not an internal
symbol, when refusing it". Broken on purpose: 1 failure.

**F4 — `na.rm = TRUE` in the ordinal intercept scores** (R/extra-indicators.R:691,
`.latents_extra_scores()`). The sum cannot meet an NA on valid input (`code == k & observed`
is FALSE on a missing cell), so the `na.rm` could only hide a real fault (an NA posterior or
probability) and return a finite, wrong score. No number changes for valid input. Guard:
test "ordinal and count scores equal the numerical gradient with missing cells" (finite and
equal to `numDeriv::grad` of the observed-data log likelihood with NA ordinal and count
cells, tolerance 1e-6). Broken on purpose (`code == k` without `& observed`): 2 failures.

## 3. Logged, not fixed

- `.multilpa_maximization()` (R/multilpa.R:73, :152) raises two unclassed `stop()`s (zero effective membership; non-finite M-step parameters); they occur only inside start loops, are recorded per start and surface as `latents_all_starts_failed`. Same for the two k-means non-convergence `stop()`s in `.multilpa_initialize()`.
- `multilpa()`/`lpa()` without `n_profiles` give R's own "argument \"n_profiles\" is missing" (unclassed); `id` and `data` have dedicated messages, `n_profiles` does not.
- `count_dispersion = "equal"` with the default Poisson `count_model` is accepted and silently has no effect (documented as an NB-only option).
- NaN in a categorical or ordinal column is treated as missing (NA), while a NaN in a continuous or count column is refused as non-finite.
- `.multilpa_single_level()` notice says "latent profile" for `multilca(id = NULL)` whose indicators are partly ordinal/count (lca() muffles it).
- `lpa()`'s `...` documentation lists `categorical` but not `ordinal`/`count` (they work; multilpa() documents them).

## 4. Needs another slice

- R/accessors.R:234 (`.multilpa_model_frame`, via `.latents_draw_extra_columns()`): `get_results(fit, "data")` returns ordinal indicators as ordered factors of the labels even when the column was numeric, and shows `NA` for a value seen only on a zero-weight row (reproduced: row 5, `o = 9`, weight 0 → `<NA>`), whereas categorical columns keep their type and restore excluded-only values (`.multilpa_restore_excluded`). Changing it changes a printed table and the refit frame, so it was not touched here.

## 5. Test runs (verbatim summary lines)

test-audit-1101-profile-fitting.R: PASS 76 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-categorical.R: PASS 75 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-categorical-inference.R: PASS 27 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-categorical-suspect.R: PASS 15 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-extra-indicators.R: PASS 52 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-multilca.R: PASS 17 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-multilpa.R: PASS 156 FAIL 0 SKIP 0 WARN 0 ERROR 0
test-single-level.R: PASS 30 FAIL 0 SKIP 0 WARN 0 ERROR 0

Engine golden comparison (`Rscript equivalence/engine-golden/compare.R`, whole working tree incl. other slices' uncommitted edits):
`identical 247, DIFFERENT 13` of 260. All 13 differences are one attribute of the bootstrap
inference table (`replicate_rows` in the baseline vs `replicate_identity` now: svc_inference_*,
*_bootstrap, profiles_two_level_VEI / profiles_lca_single / profiles_count_two_level via their
`inference_bootstrap` service). That attribute is set in R/bootstrap-inference.R (modified outside
this slice; it was already in `git status` before this audit), not in any file of this slice. No
estimate, table or text of any case changed.

`devtools::document()` not run: the only roxygen edit is in an `@noRd` block
(`.latents_check_extra_arguments` @return), which generates no Rd; NAMESPACE and man/ unaffected.
