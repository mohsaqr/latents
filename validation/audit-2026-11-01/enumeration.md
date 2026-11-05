# Audit 2026-11-01 — slice `enumeration`

Files: R/enumeration.R, R/enumeration-type.R, R/selection.R, R/search.R,
R/candidates.R, R/lmr.R, R/grid-views.R, R/sensitivity.R, R/robust.R,
R/parallel.R, R/parallel-report.R (121 top-level functions).

Checking scripts (synthetic data, run with `OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`, `devtools::load_all()`, one R process
at a time, `cores = 1` throughout) are in the session scratchpad
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/enum/`:
`ic.R` (every information criterion of a two-level VVI candidate by hand,
with the complete-data entropy rebuilt from the parameters; relative
entropies and AvePP by hand), `lmr.R` (VLMR/aLMR against an independent W
eigen-decomposition and CompQuadForm), `sel.R` (deltas, knees, smallest
shares, classification), `blrt.R` (BLRT p-value, simulation moments),
`covnest.R`/`covmix.R` (covariate pairs), `scores.R`/`scores2.R` (robust
group scores against numDeriv), `sens*.R`, `grids.R`/`fam*.R`/`reg.R`
(transition, family and regression grids), `shares.R` (weights, zero-weight
group, empty rows), `routing.R`, `methods.R`, `dollar.R`, `plots*.R`,
`par.R` (serial mapper), `lmrprior.R`. No cluster or forked worker was
started; the parallel transport code was read only.

## 1. Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `enumerate_classes` | R/enumeration.R:189 | ran synthetic (`ic.R`: all 14 criteria of a 2x2 VVI candidate by hand; `routing.R`: 12 refusals, local dependence, one continuous indicator, mixed grid); test-enumeration.R | ok |
| `.latents_candidate_call` | R/enumeration.R:392 | ran (`dollar.R`: candidate call names the caller's data) | ok |
| `.multilpa_enumeration_models` | R/enumeration.R:411 | ran (`routing.R`: pattern + EEE refused, `model` + `variance_model` refused, one continuous indicator keeps EEI/VVI) | ok |
| `get_candidate` | R/enumeration.R:553 | ran (`methods.R`, `fam.R`, `grids.R`: name, counts, pair, model, transition and family paths) | fixed |
| `.latents_check_candidate_counts` | R/enumeration.R:669 | ran (`methods.R`: 2.5 gave an unclassed error) | fixed |
| `candidate_fit` | R/enumeration.R:687 | ran (`methods.R`: forwards to get_candidate) | ok |
| `summary.multilpa_enumeration` | R/enumeration.R:702 | ran (`methods.R`, LCA grid) | ok |
| `.multilpa_enumeration_minima` | R/enumeration.R:729 | ran (`ic.R`: best rows equal which.min over converged candidates by hand) | ok |
| `.multilpa_enumeration_grid` | R/enumeration.R:767 | ran (`ic.R`: AvePP min/max/mean/below at both levels by hand) | ok |
| `.multilpa_avepp_summary` | R/enumeration.R:814 | ran (`ic.R`, hand AvePP) | ok |
| `.multilpa_informative_classification` | R/enumeration.R:831 | ran (`sel.R`: classification table, one-class levels dropped) | ok |
| `.multilpa_enumeration_classification` | R/enumeration.R:855 | ran (`sel.R`) | ok |
| `print.summary_multilpa_enumeration` | R/enumeration.R:880 | ran (`methods.R`) | ok |
| `.multilpa_compact_candidates` | R/enumeration.R:924 | ran via print (`methods.R`: single-level BIC pair shown once) | ok |
| `as.data.frame.summary_multilpa_enumeration` | R/enumeration.R:950 | ran (`methods.R`) | ok |
| `.multilpa_enumeration_patterned` | R/enumeration.R:960 | ran (`routing.R`, `local_dependence = "y1:y3"`) | ok |
| `.multilpa_enumeration_arguments` | R/enumeration.R:972 | ran (`routing.R`: VVI pattern fitted with `variance_model = "varying"`) | ok |
| `.multilpa_simulate` | R/enumeration.R:983 | ran (`blrt.R`: 200 draws from a VVV fit, mixture mean and covariance match the model's to 2 decimals; group ids kept) | ok |
| `.multilpa_cov_simulate` | R/enumeration.R:1010 | ran via `bootstrap_lrt()` on covariate fits (`covnest.R`); read | ok |
| `.multilpa_carry_missingness` | R/enumeration.R:1041 | read; covered by test-bootstrap-lrt-extended.R (mocked call count) | ok |
| `.multilpa_draw_indicators` | R/enumeration.R:1054 | ran via `.multilpa_simulate()` (full covariance) | ok |
| `.multilpa_bootstrap_constraint` | R/enumeration.R:1133 | read; covered by test-bootstrap-constraints.R | ok |
| `.multilpa_check_covariate_nesting` | R/enumeration.R:1208 | ran (`covnest.R`: same names, other values; affine values; group-class step) | fixed |
| `.multilpa_cov_check_likelihood` | R/enumeration.R:1260 | ran (`covnest.R`: refuses a fit whose data differ) | ok |
| `bootstrap_lrt` | R/enumeration.R:1489 | ran (`blrt.R`: p = (1+b)/(V+1) = 0.5 by hand; `covnest.R`; `covmix.R`: plain/covariate pairs refused) | fixed |
| `print.multilpa_bootstrap_lrt` | R/enumeration.R:1782 | ran (`blrt.R`) | ok |
| `summary.multilpa_bootstrap_lrt` | R/enumeration.R:1824 | ran (`plots.R`) | ok |
| `.multilpa_bootstrap_test_frame` | R/enumeration.R:1847 | ran (`blrt.R`, as.data.frame) | ok |
| `print.summary_multilpa_bootstrap_lrt` | R/enumeration.R:1867 | ran (`plots.R`) | ok |
| `as.data.frame.multilpa_bootstrap_lrt` | R/enumeration.R:1887 | ran (`blrt.R`) | ok |
| `as.data.frame.summary_multilpa_bootstrap_lrt` | R/enumeration.R:1894 | read | ok |
| `plot.multilpa_bootstrap_lrt` | R/enumeration.R:1901 | ran (`plots.R`: returns a ggplot; visual output unchanged, not rendered) | ok |
| `.multilpa_enumeration_indices` | R/enumeration.R:1949 | ran (`ic.R`: AIC, KIC, BIC, SABIC, CAIC, AWE, ICL, CLC at both n conventions equal hand values, all.equal TRUE at 1e-10; EN computed independently from the parameters) | ok |
| `.multilpa_enumeration_criteria` | R/enumeration.R:1973 | read | ok |
| `.latents_refuse_shared_data_error` | R/enumeration.R:1996 | ran (`routing.R`: misspelled column -> latents_bad_data) | ok |
| `.latents_check_enumeration_type` | R/enumeration-type.R:7 | ran (`routing.R`) | logged |
| `.latents_enumeration_measurement_type` | R/enumeration-type.R:44 | ran (`routing.R`: mixed grid -> "mixed") | ok |
| `.latents_enumeration_type_label` | R/enumeration-type.R:55 | ran via print (`methods.R`: "(LCA)") | ok |
| `.latents_selection_views` | R/selection.R:15 | read | ok |
| `.latents_selection_frame` | R/selection.R:27 | ran via deltas/knees (`sel.R`) | ok |
| `.latents_smaller_neighbour` | R/selection.R:54 | ran (`sel.R`: deltas against hand differences) | ok |
| `.latents_deltas_table` | R/selection.R:77 | ran (`sel.R`: delta_profile and delta_group_class of profiles3_groups2 by hand) | ok |
| `.latents_kneedle` | R/selection.R:110 | ran (`sel.R`: perpendicular distances 0, 0.537, 0.251, 0 by hand; knee at 2 profiles) | ok |
| `.latents_knees_table` | R/selection.R:141 | ran (`sel.R`, `methods.R`) | ok |
| `.latents_selection_look` | R/selection.R:212 | ran via the views (`plots.R`) | ok |
| `.latents_class_word` | R/selection.R:226 | ran (LCA grid names "classes") | ok |
| `.gg_selection_deltas` | R/selection.R:239 | ran (`plots.R`: profile, transition and regression grids) | ok |
| `.gg_selection_heatmap` | R/selection.R:302 | ran (`plots.R`) | logged |
| `.gg_selection_smallest` | R/selection.R:332 | ran (`plots.R`) | ok |
| `.gg_selection_trust` | R/selection.R:383 | ran (`plots.R`) | ok |
| `.latents_trust_label` | R/selection.R:435 | read; covered by test-grid-contract.R | ok |
| `.gg_selection_view` | R/selection.R:451 | ran (`plots.R`) | ok |
| `.latents_search_defaults` | R/search.R:21 | read | ok |
| `.latents_check_search` | R/search.R:31 | read; covered by test-default-search.R | ok |
| `.latents_default_search` | R/search.R:58 | read; covered by test-default-search.R (frozen) | logged |
| `.latents_search_screen` | R/search.R:75 | read (pool, protect, `n_screen` refusals) | ok |
| `.latents_start_kind` | R/search.R:113 | read | ok |
| `.latents_start_kind_label` | R/search.R:122 | read | ok |
| `.latents_search_record` | R/search.R:140 | read | ok |
| `.latents_refit_controls` | R/search.R:161 | ran (`blrt.R`: refit_controls table repeats each fit's controls) | ok |
| `.latents_refit_controls_frame` | R/search.R:183 | ran (`blrt.R`) | ok |
| `.latents_search_columns` | R/search.R:209 | read | ok |
| `.latents_final_stage` | R/search.R:224 | read | ok |
| `.latents_candidate_names` | R/candidates.R:20 | ran (`dollar.R`, `routing.R`, `fam.R`) | logged |
| `.latents_candidate_data` | R/candidates.R:48 | ran (`dollar.R`) | logged |
| `.latents_is_grid` | R/candidates.R:67 | read | ok |
| `.latents_candidate_row` | R/candidates.R:81 | ran (`methods.R`: unknown name, unknown pair) | ok |
| `.latents_candidate_at` | R/candidates.R:132 | ran | ok |
| `.latents_named_table` | R/candidates.R:148 | ran | ok |
| `.latents_grid_dollar` | R/candidates.R:165 | ran (`dollar.R`: candidate, partial component, unknown name -> NULL) | ok |
| `$.multilpa_enumeration` | R/candidates.R:175 | ran (`dollar.R`) | ok |
| `$.latents_transition_enumeration` | R/candidates.R:179 | ran (`grids.R`) | ok |
| `$.latents_family_enumeration` | R/candidates.R:183 | read | ok |
| `.latents_grid_dollar_names` | R/candidates.R:190 | ran (`dollar.R`) | ok |
| `.DollarNames.multilpa_enumeration` | R/candidates.R:197 | ran (`dollar.R`) | ok |
| `.DollarNames.latents_transition_enumeration` | R/candidates.R:203 | read | ok |
| `.DollarNames.latents_family_enumeration` | R/candidates.R:209 | read | ok |
| `.latents_grid_result_tables` | R/candidates.R:226 | ran (`grids.R`, `fam.R`: minimum over converged) | ok |
| `.latents_vuong_parts` | R/lmr.R:24 | ran (`lmr.R`) | ok |
| `.latents_vuong_weights` | R/lmr.R:52 | ran (`lmr.R`: weights reproduce mean/sd of the eigenvalues of W built independently from A, B, B10 with solve(); p-values equal CompQuadForm::imhof/davies) | ok |
| `.latents_imhof_upper` | R/lmr.R:189 | ran (`lmr.R`: chi-square(3) = pchisq; mixed-sign weights, q = 0 and q < 0 equal CompQuadForm::imhof to 7 digits) | ok |
| `.latents_wynn` | R/lmr.R:279 | read; exercised by `.latents_imhof_upper()` | ok |
| `.latents_vuong_row` | R/lmr.R:306 | ran (`lmr.R`: adjusted LR = LR / (1 + 1/(d log n)), Lo, Mendell & Rubin 2001 eq. 15, all.equal TRUE) | ok |
| `.latents_lmr_table` | R/lmr.R:341 | ran (`lmr.R`, `methods.R`, `lmrprior.R`: boundary and prior notes, weights, noise) | logged |
| `.latents_grid_kind` | R/grid-views.R:20 | ran | ok |
| `.latents_grid_view_names` | R/grid-views.R:33 | read | ok |
| `.latents_grid_views` | R/grid-views.R:58 | ran (`fam.R`, `reg.R`) | ok |
| `.latents_grid_view_table` | R/grid-views.R:100 | ran (`fam.R`: plot_views) | ok |
| `.latents_refuse_non_grid` | R/grid-views.R:137 | read; covered by test-grid-contract.R | ok |
| `.latents_grid_offer_text` | R/grid-views.R:148 | ran (in a refusal) | ok |
| `.latents_require_grid_views` | R/grid-views.R:164 | ran (`fam2.R`: family deltas -> latents_unsupported_view) | ok |
| `.latents_grid_table` | R/grid-views.R:188 | ran (`reg.R`) | ok |
| `.latents_regression_posteriors` | R/grid-views.R:213 | ran (`reg.R`) | ok |
| `.latents_grid_quality_table` | R/grid-views.R:239 | ran (`grids.R`, `fam.R`, `reg.R`: regression AvePP minima by hand) | ok |
| `.latents_grid_criteria_titles` | R/grid-views.R:304 | ran (`reg.R`: bic_rows equal to bic drawn once) | ok |
| `.latents_grid_smallest_shares` | R/grid-views.R:329 | ran (`sel.R`, `shares.R` with weights, a zero-weight group and empty rows, `grids.R`, `reg.R`: equal min/sum of weighted posterior counts by hand) | ok |
| `.latents_plot_former_grid` | R/grid-views.R:378 | ran (`plots2.R`: `criterion = "bic_i"`, unknown argument warns) | ok |
| `sensitivity` | R/sensitivity.R:120 | ran (`sens.R`, `sens2.R`) | logged |
| `as.data.frame.latents_sensitivity` | R/sensitivity.R:215 | ran (`sens.R`) | ok |
| `plot.latents_sensitivity` | R/sensitivity.R:223 | ran (`plots.R`) | ok |
| `.sensitivity_profiles` | R/sensitivity.R:295 | ran (`sens3.R`: agreement equals the best label permutation by hand, 0.985) | ok |
| `.sensitivity_transitions` | R/sensitivity.R:331 | read; covered by test-sensitivity.R | ok |
| `.sensitivity_growth` | R/sensitivity.R:369 | read; covered by test-sensitivity.R | ok |
| `.multilpa_optimum_index` | R/sensitivity.R:411 | ran (`sens.R`: anchored grouping of a hand vector, NA kept, single value) | ok |
| `.multilpa_group_scores` | R/robust.R:13 | ran (`scores.R`: column sums equal numDeriv gradient of the log likelihood, VVI/EEI/VVV/EEE/VEI/EVE with a categorical item and FIML, max diff 5e-7; `scores2.R`: single group rows equal that group's numerical gradient, 1e-8) | ok |
| `.multilpa_full_measurement_group_score` | R/robust.R:82 | ran (`scores.R`, VVV/EEE incl. FIML) | ok |
| `.latents_check_cores` | R/parallel.R:12 | ran (`par.R`: 2.5 refused, 1 accepted, 1e10) | logged |
| `.latents_fork_safe` | R/parallel.R:46 | ran (`par.R`) | ok |
| `.latents_stop_cluster` | R/parallel.R:59 | read (no cluster started) | ok |
| `.latents_collect_conditions` | R/parallel.R:75 | ran via the serial mapper (`par.R`) | ok |
| `.latents_worker_context` | R/parallel.R:106 | read (no worker started) | ok |
| `.latents_task_collector` | R/parallel.R:137 | ran (`par.R`: a task's classed error reaches the caller) | ok |
| `.latents_signal_conditions` | R/parallel.R:153 | ran (`par.R`: warnings and messages in task order) | ok |
| `.latents_fork_pass` | R/parallel.R:168 | read (no fork started) | ok |
| `.latents_windows_alive` | R/parallel.R:220 | read | ok |
| `.latents_make_cluster` | R/parallel.R:242 | read (no cluster started) | ok |
| `.latents_abort_cluster` | R/parallel.R:267 | read | ok |
| `.latents_socket_pass` | R/parallel.R:318 | read (chunking, retry, cleanup) | ok |
| `.latents_warn_worker_lost` | R/parallel.R:374 | read; covered by test-parallel-report.R | ok |
| `.latents_map_candidates` | R/parallel.R:413 | ran with cores = 1 (`par.R`: costs reorder dispatch and results come back in index order; signal = FALSE; empty input); retry indexing read | ok |
| `.latents_worker_reporter` | R/parallel-report.R:7 | covered by test-parallel-report.R (run) | ok |

Counts: 121 functions checked; 4 marked fixed (two defects, each spanning two
functions); 8 logged; 0 needs-other-slice in the table (other-slice findings
are in section 4).

## 2. Fixed

**F1. `bootstrap_lrt()` on covariate fits that name the same covariate with
other values was refused as bad data, not as incomparable models**
(R/enumeration.R:1208-1253, `.multilpa_check_covariate_nesting()`).
Reproduction (`covnest.R`): a 2-profile and a 3-profile `profile_covariates =
"z"` fit, the second to a copy of the data with `z` redrawn. Before the fix:
`bootstrap_lrt(null, alt)` (no `data`) raised `latents_bad_inference_data`
"data do not reproduce the fitted model likelihood." — the caller supplied no
data. So the pair could NOT be accepted (the per-model likelihood check
always caught it; the affine case `2 z + 1` was caught too), but the refusal
named the wrong cause and class. Cause: `.multilpa_check_covariate_nesting()`
compared covariate names only. Change: after the name checks, the values are
compared with compare-lrt.R's `.latents_lrt_spans()` (called, not
duplicated) and a mismatch raises `latents_incomparable_models`, the class
`bootstrap_lrt()` already documents for incompatible data. The
`.multilpa_lrt_predictors_nested()` wrapper could not be called as is: it
stacks the profile designs, which have n x C rows, so a valid group-class
step (1 -> 2 group classes, same `z`) was reported as non-nested ("a profile
membership covariate of the first has values the second's covariates do not
reproduce", `covnest.R`). Group class 1's design block (n rows, holds the
covariates under either slope rule) is compared instead, plus the group
design when the null has group covariates. Valid pairs are unaffected (same
data: accepted; profile step and group-class step pass the check). Test:
test-audit-1101-enumeration.R "bootstrap_lrt() refuses a shared covariate
with other values as incomparable" and "the covariate value check passes
nested pairs at either level"; with the check disabled the first test fails
with the old `latents_bad_inference_data` message.

**F2. `get_candidate()` argument checks were unclassed, and a family grid
reported a fractional count truncated** (R/enumeration.R:573, 611, 669).
Reproduction (`methods.R`, `fam.R`): `get_candidate(lca_grid, 2.5)` raised a
`simpleError` from `stopifnot()`; `get_candidate(family_grid,
n_group_classes = 2.5, model = "additive_equal")` raised
`latents_unknown_candidate` "0 candidates have 2 group classes and model
additive_equal" although that candidate exists; a vector `model` on a family
grid built one message per element. Cause: `.latents_check_candidate_counts()`
used `stopifnot()`, the family path skipped it, and `model` was validated only
on the profile path. Change: the helper uses `.latents_require()`
(`latents_bad_argument`) with `n_profiles = 1L` as default so the family path
can check its count; the `model` check runs once for every grid kind; the
Conditions section of `?get_candidate` names the class. Same arguments,
same returns for valid input. Test: test-audit-1101-enumeration.R
"get_candidate() refuses a fractional count or a vector model ..." and "a
family grid checks its count ..."; both fail with the old code.

## 3. Logged, not fixed

- `.latents_check_cores(1e10)`: `as.integer()` coercion warning, then an unclassed "missing value" error (junk value; a realistic `cores` is refused with `latents_bad_argument`).
- `.latents_check_enumeration_type()`: `match.arg()` accepts a partial `type = "lp"` as `"lpa"`; the documentation lists full strings only. Harmless.
- `.latents_candidate_names()`: with `local_dependence` the model label enters the name, e.g. `"profiles2_VVI + y1:y3"`; `grid$<name>` then needs backticks. Works with `get_candidate()`.
- `.latents_candidate_data()`: a grid made by `enumerate_lpa()`/`enumerate_lca()` records `call = enumerate_classes(data = mydata, vars = vars, n_profiles = n_profiles, model = model, ...)` (wrapper formals, not values). Display only; each candidate's own call is correct.
- `.latents_lmr_table()`: on a `noise = TRUE` grid the whole `get_results(grid, "lmr")` is refused with `latents_unsupported_noise` (classed) instead of NA rows with a note, as prior and boundary candidates get.
- `.gg_selection_heatmap()`: tile labels use `"+%.0f"`, so a candidate within 0.5 of the lowest reads "+0".
- `sensitivity()`: with `n_starts = 1` under the kmeans search the only start is Ward's deterministic start, so every seed returns the identical refit and agreement 1 (`sens.R`: 6 seeds, identical log likelihood -444.3303). True of the search, not a code defect; the help could say so.
- `sensitivity()`: the comment "a covariate fit is class `multilpa`" is stale (covariate fits are class `multilpa_covariates` only); the refusal itself is correct.
- `.latents_default_search()`: the brief describes the default as kmeans with emEM opt-in; the code defaults EII/EEI/EEE to emEM (NEWS 1.1.11 and test-default-search.R pin this). Not touched (frozen).
- `bootstrap_lrt()`: direct-effect covariates with the same name but other values are still refused by the likelihood check as `latents_bad_inference_data` (F1 covers membership covariates). A pair whose covariate differs by an affine map passes the column-space check and is then refused by the likelihood check, correctly, since the two fits saw different data.
- `.multilpa_optimum_index()`: an all-`NA` logical vector fails its `is.numeric` guard; unreachable from `sensitivity()`, which passes `numeric` from `vapply()`.

## 4. Needs another slice

- R/additive-compare.R:286 and :288, `get_results.latents_family_enumeration()`: `match.arg(what)` / `match.arg(criterion)` raise an unclassed `simpleError` ("'arg' should be one of ...") for e.g. `what = "deltas"` or `"lmr"`; the transition-grid method (R/lta-compare.R:575) uses the classed `.latents_match_formal()`.

## 5. Test runs

One R process, `testthat::test_file(reporter = "silent")` per file, `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`, after the fixes. Not run: test-search.R and test-parallel*.R other than test-parallel-report.R (they start socket workers; brief says skip). The one skip in test-sensitivity.R is its parallel-seeds test (`skip_if(!.latents_fork_safe())`, Accelerate BLAS).

```
test-audit-1101-enumeration.R: PASS 13 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (16s) 
test-enumeration.R: PASS 315 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (53s) 
test-enumeration-type.R: PASS 122 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (4s) 
test-enumerate-verbs.R: PASS 30 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (4s) 
test-default-search.R: PASS 28 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (5s) 
test-lmr.R: PASS 103 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (13s) 
test-features-vlmr.R: PASS 176 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (26s) 
test-sensitivity.R: PASS 76 FAIL 0 SKIP 1 WARN 0 ERROR-tests 0 (26s) 
test-robust.R: PASS 35 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (2s) 
test-bootstrap-lrt-extended.R: PASS 24 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (58s) 
test-bootstrap-constraints.R: PASS 21 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (3s) 
test-grid-contract.R: PASS 237 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (49s) 
test-parallel-report.R: PASS 20 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (0s) 
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (5s) 
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (1s) 
test-review-1027-error-contract.R: PASS 216 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (3s) 
test-additive-compare.R: PASS 29 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (6s) 
test-audit-recent-api.R: PASS 44 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (9s) 
test-tidy-classes.R: PASS 61 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (4s) 
test-release-surface.R: PASS 88 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (6s) 
test-review-integration.R: PASS 39 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (10s) 
test-audit-reporting-repairs.R: PASS 25 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (8s) 
test-seeds.R: PASS 26 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (2s) 
test-weights.R: PASS 76 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (10s) 
test-direct-effects.R: PASS 55 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (7s) 
test-get-data.R: PASS 230 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (10s) 
test-tree.R: PASS 131 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (14s) 
test-plot-named.R: PASS 138 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (25s) 
test-covariance-structure.R: PASS 167 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (8s) 
test-informative-rows.R: PASS 57 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (37s) 
test-review-0970-estimation.R: PASS 27 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (3s) 
test-review-1015-inference.R: PASS 39 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (2s) 
test-lta-compare.R: PASS 31 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (53s) 
test-scaling.R: PASS 31 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (2s) 
test-features-joint-items.R: PASS 61 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (11s) 
test-audit-plots-predictions.R: PASS 86 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (4s) 
test-extra-indicators.R: PASS 52 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (21s) 
test-features-pair-covariance.R: PASS 69 FAIL 0 SKIP 0 WARN 0 ERROR-tests 0 (6s) 
```

Sabotage: with each fix disabled (value check forced FALSE; `stopifnot` restored in the count helper; family count check and shared `model` check removed) test-audit-1101-enumeration.R gave 3 errors (lines 34, 59, 78), then passed again with the fixes restored.

`devtools::document()`: NAMESPACE unchanged; man/get_candidate.Rd carries the new Conditions sentence (the only man/ change from this slice).

## 6. Engine golden comparison

`Rscript equivalence/engine-golden/compare.R` (all cases, one process, after
the fixes; it rewrote `equivalence/engine-golden/last-compare.csv`, as it
always does): 35 cases DIFFERENT, exit 1, none from this slice.

- 13: the known bootstrap attribute rename (`replicate_rows` ->
  `replicate_identity`): profiles_two_level_VEI, profiles_lca_single,
  profiles_count_two_level, lta_homogeneous_bootstrap, lta_general_bootstrap,
  lta_general_weights_bootstrap, svc_inference_lpa, svc_inference_lpa_full,
  svc_inference_two_level, svc_inference_bootstrap_lca,
  svc_inference_bootstrap_mixed_extra, svc_inference_bootstrap_vei,
  svc_inference_bootstrap_staged.
- 22 OTHER: every reg_* Gaussian/binomial/Poisson/NB/ordinal case of a
  single-level or group-level mixture regression plus growth_lcga, all with
  `first_diff_path = <case>/probes/tables/value/group_classes@display`,
  "length 4 vs 5": the printed `group_classes` table of a mixture regression
  now declares five display columns (R/mixture-regression-methods.R:495 adds
  "Groups assigned"). No file of this slice builds that table.
- Every case through this slice's code is identical: svc_enumerate_lpa_basic,
  svc_enumerate_lpa_all_structures, svc_enumerate_lca,
  svc_enumerate_classes_two_level, svc_candidate_fit,
  svc_candidate_fit_ambiguous, all seven svc_bootstrap_lrt_* (including
  covariates and fixed), family_enumeration, family_bootstrap_lrt(_refused),
  lta_enumeration*, reg_enumerate*, svc_sensitivity_two_level,
  svc_sensitivity_options.
