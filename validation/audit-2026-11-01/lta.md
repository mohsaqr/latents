# Audit 2026-11-01 — slice `lta`

Files: R/lta-engine.R, R/lta-starts.R, R/lta-methods.R, R/lta-inference.R,
R/lta-bootstrap.R, R/lta-compare.R, R/lta-three-step.R, R/sequences.R, R/tna.R.
149 functions checked (147 definitions + 2 method aliases), 6 fixed, 4 logged,
1 needs another slice. No reported number changed; no engine, default, signature
or table layout changed; no roxygen changed (devtools::document() not needed).

Scratch evidence scripts (not part of the package):
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/lta/` (`sim.R`, `brute.R`, `scores.R`).

## Independent checks of the engine (no defect found)

- **Forward-backward vs brute-force path sum.** An independent likelihood that
  enumerates every latent path (K^T paths per person, initial logits against
  the last profile, transition logits against staying, covariates read at the
  destination occasion, emissions by dnorm) was compared with
  `fit$log_likelihood`. Differences (engine − brute):
  homogeneous 0, 2 group classes 1.1e-13, occasion transitions 0,
  transition + initial covariates with 2 classes 0, occasion measurement −1.1e-13,
  second order 0, second order + covariate 1.1e-13, stayer 0; on T = 4 with
  20 of 60 persons dropping the last occasion and rows shuffled: occasion +
  covariates −1.1e-13, second order + covariate + 2 classes −1.1e-13,
  stayer + covariate + occasion measurement 0. This also shows the covariate
  designs line up with occasions (destination row) after shuffling and dropout.
- **Row-stochastic transitions.** `get_results(fit, "transitions")` rows sum to
  1 per class × occasion × origin in every structure above.
- **Analytic scores.** `colSums(.lta_group_scores())` vs central differences of
  the exact log-likelihood at the fitted point: max abs diff 2.1e-8 (occasion),
  2.5e-8 (categorical + equal variance), 2.1e-8 (covariates, 2 classes),
  3.2e-8 (second order), 2.2e-8 (stayer + occasion measurement), 2.7e-8 (full
  covariance), 2.0e-8 (sampling weights), 2.1e-8 (FIML). `.lta_pack(.lta_unpack())`
  round trip ≤ 2.2e-16; `length(.lta_names()) == length(.lta_theta())` in all.
- **Covariance.** `vcov()` SEs vs `optimHess` of the exact likelihood: relative
  difference 0 for measurement and transition logits, ≤ 2.6e-4 for an initial
  covariate held at mean 50 (raw-scale differencing in the check, not the package).
- **n_parameters** hand-counted: 21 (2 classes, initial + transition covariate),
  39 (stayer, 4 occasion-specific blocks, transition covariate), 17
  (categorical + equal variance, occasion transitions) — all match.
- **Simulation.** `.lta_simulate()` from a second-order fit with transition and
  initial covariates and a categorical indicator (n = 1500), refitted: every
  coefficient within sampling error of the generating fit.
- **Three-step.** `.lta_occasion_posteriors()` vs a hand computation
  w_k f_k(y_it): max diff 2.2e-16 (homogeneous engine) and 3.3e-16 (general
  engine, indicator at +50). `.lta_three_step_evaluate()` scores vs central
  differences: 6.8e-8.
- **`.lta_nested()` follow-up.** Two fits sharing covariate `z` with different
  values (null `z`, alternative `z, w` fitted to data with `z` changed):
  `.lta_nested()` is TRUE, `.lta_lrt_designs_nested()` is FALSE, and
  `bootstrap_lrt()` refuses with `latents_bad_inference_data` whichever data is
  supplied, because `.lta_check_data()` rebuilds both fits' designs from the one
  `data` and compares them exactly (tolerance 0) before any simulation. The path
  therefore cannot accept such a pair, so no call to
  `.lta_lrt_designs_nested()` was added; the refusal is pinned by a regression
  test.

## Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.lta_covariate_matrix` | R/lta-engine.R:15 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_designs` | R/lta-engine.R:70 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_n_intercepts` | R/lta-engine.R:106 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_design_chart` | R/lta-engine.R:130 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_chart_designs` | R/lta-engine.R:154 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_chart_matrix` | R/lta-engine.R:171 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_map_coefficients` | R/lta-engine.R:204 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_chart_jacobian` | R/lta-engine.R:237 | ran synthetic: vcov() SEs vs optimHess of the exact likelihood with analytic gradient (covariate held at mean 50), rel. diff <= 2.6e-4 (raw-scale differencing), others 0; test-transition-inference | ok |
| `.lta_check_transition_rank` | R/lta-engine.R:282 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_check_design_rank` | R/lta-engine.R:332 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_check_covariate_gaps` | R/lta-engine.R:374 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_log_transition` | R/lta-engine.R:405 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_log_initial` | R/lta-engine.R:421 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_moments` | R/lta-engine.R:430 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_reduce_moments` | R/lta-engine.R:463 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_log_density` | R/lta-engine.R:495 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_augment` | R/lta-engine.R:528 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_identity_transition` | R/lta-engine.R:557 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_identity_transition2` | R/lta-engine.R:565 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_log_transition2` | R/lta-engine.R:576 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_expectation` | R/lta-engine.R:592 | ran synthetic: path-sum likelihood over all latent paths vs engine, 11 structures (homogeneous, 2 group classes, occasion, covariates, occasion measurement, 2nd order +/- covariate, stayer, dropout + shuffled rows), max diff 1.1e-13; test-audit-1101-lta | ok |
| `.lta_maximization` | R/lta-engine.R:667 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_row_weights` | R/lta-engine.R:773 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_em` | R/lta-engine.R:786 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_initialize` | R/lta-engine.R:828 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_fit_general` | R/lta-engine.R:889 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_setup_general` | R/lta-engine.R:912 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_run_general_start` | R/lta-engine.R:1064 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_finalise_general` | R/lta-engine.R:1103 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_used_rows` | R/lta-engine.R:1249 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_occupied_occasions` | R/lta-engine.R:1275 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_identified_columns` | R/lta-engine.R:1327 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_empty_origins` | R/lta-engine.R:1362 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_zero_empty_origins` | R/lta-engine.R:1389 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_identified_blocks` | R/lta-engine.R:1417 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_count_parameters` | R/lta-engine.R:1423 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_result` | R/lta-engine.R:1441 | ran synthetic: fits of the 11 structures above; n_parameters hand-checked (21 cov+2 classes, 39 stayer+occasion meas., 17 categorical/equal); test-lta-general, test-review-1026-lta* | ok |
| `.lta_bind_context` | R/lta-starts.R:7 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_plan` | R/lta-starts.R:33 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_run_plan` | R/lta-starts.R:50 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_restore_structure_log` | R/lta-starts.R:65 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_report_failed_plan` | R/lta-starts.R:82 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_finish_plan` | R/lta-starts.R:96 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_execute_start` | R/lta-starts.R:110 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_start_tasks` | R/lta-starts.R:118 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_setup_task` | R/lta-starts.R:127 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_pool_phase` | R/lta-starts.R:139 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_candidate_error` | R/lta-starts.R:147 | read; exercised by every lta() fit and enumerate_classes(time =) run here (cores = 1); test-lta-parallel | ok |
| `.lta_tables` | R/lta-methods.R:5 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_coefficient_table` | R/lta-methods.R:20 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_probability_table` | R/lta-methods.R:87 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_second_order_table` | R/lta-methods.R:148 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_occasion_weights` | R/lta-methods.R:201 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_terms_used` | R/lta-methods.R:214 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_logits_identified` | R/lta-methods.R:230 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_group_mean` | R/lta-methods.R:243 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_responses_table` | R/lta-methods.R:259 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_assignment_view` | R/lta-methods.R:289 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_assignments` | R/lta-methods.R:314 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `.lta_table` | R/lta-methods.R:347 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `get_results.multilpa_lta` | R/lta-methods.R:396 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `as.data.frame.multilpa_lta` | R/lta-methods.R:424 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `print.multilpa_lta` | R/lta-methods.R:431 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `summary.multilpa_lta` | R/lta-methods.R:468 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `print.summary_multilpa_lta` | R/lta-methods.R:489 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `get_results.summary_multilpa_lta` | R/lta-methods.R:513 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `as.data.frame.summary_multilpa_lta` | R/lta-methods.R:525 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `coef.multilpa_lta` | R/lta-methods.R:531 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `vcov.multilpa_lta` | R/lta-methods.R:548 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `logLik.multilpa_lta` | R/lta-methods.R:559 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `nobs.multilpa_lta` | R/lta-methods.R:567 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `parameter_inference.multilpa_lta` | R/lta-methods.R:571 | ran synthetic: every get_results(fit, what) table on covariate, occasion and categorical-only fits; transition rows sum to 1 (max err < 1e-12); test-audit-1101-lta, test-lta-general | ok |
| `plot.multilpa_lta` | R/lta-methods.R:624 | read; covered by test-lta-general / plot snapshot tests (no visual change made) | ok |
| `.lta_transition_columns` | R/lta-inference.R:15 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_uses_covariance_chart` | R/lta-inference.R:29 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_covariance_view` | R/lta-inference.R:40 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_covariance_coordinates` | R/lta-inference.R:55 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_covariance_names` | R/lta-inference.R:67 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_names` | R/lta-inference.R:77 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_pack` | R/lta-inference.R:147 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_unpack` | R/lta-inference.R:195 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_group_scores` | R/lta-inference.R:281 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_chart_scores` | R/lta-inference.R:412 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_rowsum` | R/lta-inference.R:445 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_inference_scales` | R/lta-inference.R:462 | ran synthetic: vcov() SEs vs optimHess of the exact likelihood with analytic gradient (covariate held at mean 50), rel. diff <= 2.6e-4 (raw-scale differencing), others 0; test-transition-inference | ok |
| `.lta_inference` | R/lta-inference.R:485 | ran synthetic: vcov() SEs vs optimHess of the exact likelihood with analytic gradient (covariate held at mean 50), rel. diff <= 2.6e-4 (raw-scale differencing), others 0; test-transition-inference | ok |
| `.lta_identified_coordinates` | R/lta-inference.R:609 | read + covered by test-review-1026-lta-open, test-lta-thin-profiles; estimated column checked in synthetic tables | ok |
| `.lta_empty_origin_names` | R/lta-inference.R:640 | read + covered by test-review-1026-lta-open, test-lta-thin-profiles; estimated column checked in synthetic tables | ok |
| `.lta_theta` | R/lta-inference.R:672 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_parameters` | R/lta-inference.R:678 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_quasi_newton` | R/lta-inference.R:716 | ran synthetic: every general fit above (finish converged, max|score| <= 2.3e-4 on full-cov EM-only path, <= 1.2e-5 otherwise) | ok |
| `.lta_probability_boundary` | R/lta-inference.R:770 | ran synthetic: boundary warnings on 2nd-order/stayer fits with zero moves; read | ok |
| `.lta_extra_names` | R/lta-inference.R:806 | ran synthetic: analytic scores vs central differences of the log-likelihood, 8 structures (occasion, categorical+equal var, covariates+2 classes, 2nd order, stayer+occasion measurement, full covariance, weights, FIML), max diff 3.2e-8; pack/unpack round trip <= 2.2e-16 | ok |
| `.lta_bootstrap_inference` | R/lta-bootstrap.R:23 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_bootstrap_frame` | R/lta-bootstrap.R:93 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_bootstrap_estimates` | R/lta-bootstrap.R:136 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_bootstrap_replicates` | R/lta-bootstrap.R:152 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_align_estimates` | R/lta-bootstrap.R:183 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_profile_signature` | R/lta-bootstrap.R:213 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_class_signature` | R/lta-bootstrap.R:235 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_relabel_profiles` | R/lta-bootstrap.R:272 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_relabel_classes` | R/lta-bootstrap.R:311 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_permute_transitions` | R/lta-bootstrap.R:325 | read; label alignment checked against LEARNINGS (stay-referenced logits rename, initial rebased); covered by test-lta-bootstrap, test-review-1101-bootstrap-tables | ok |
| `.lta_bootstrap_dispatch` | R/lta-bootstrap.R:359 | ran synthetic: iter = 1 gave an unclassed simpleError; test-audit-1101-lta | fixed |
| `.lta_arguments` | R/lta-compare.R:7 | ran synthetic: bootstrap_lrt and simulation on homogeneous and general fits | ok |
| `.lta_general_view` | R/lta-compare.R:31 | ran synthetic: bootstrap_lrt and simulation on homogeneous and general fits | ok |
| `.lta_simulate` | R/lta-compare.R:79 | ran synthetic: simulate from 2nd-order + covariates + categorical fit (n = 1500) and refit; estimates recovered within sampling error | ok |
| `.lta_draw_continuous` | R/lta-compare.R:167 | read; ran through .lta_simulate (diagonal) and test-lta-covariance-boundary (full) | ok |
| `.lta_nested` | R/lta-compare.R:197 | ran synthetic: shared covariate z with different values in the two fits -> .lta_nested TRUE, .lta_lrt_designs_nested FALSE, but bootstrap_lrt refuses via .lta_check_data (latents_bad_inference_data) for either data; test-audit-1101-lta | logged |
| `.lta_refit` | R/lta-compare.R:240 | ran synthetic: bootstrap_lrt refusals; covered by test-lta-compare, test-bootstrap-lrt-extended | ok |
| `.lta_check_data` | R/lta-compare.R:267 | ran synthetic: bootstrap_lrt refusals; covered by test-lta-compare, test-bootstrap-lrt-extended | ok |
| `.lta_bootstrap_lrt` | R/lta-compare.R:318 | ran synthetic: bootstrap_lrt refusals; covered by test-lta-compare, test-bootstrap-lrt-extended | logged |
| `.lta_enumerate` | R/lta-compare.R:414 | ran synthetic: enumerate_classes(time =, n_profiles = 2:3) rows equal direct lta() fits (logLik -380.5727, BIC 834.923); n_profiles = c(2, NA) / "a" / integer(0) crashed unclassed; test-audit-1101-lta | fixed |
| `get_results.latents_transition_enumeration` | R/lta-compare.R:568 | ran synthetic on a 2-candidate grid; test-audit-1101-lta | fixed |
| `as.data.frame.latents_transition_enumeration` | R/lta-compare.R:585 | ran synthetic on a 2-candidate grid; test-audit-1101-lta | ok |
| `print.latents_transition_enumeration` | R/lta-compare.R:592 | ran synthetic on a 2-candidate grid; test-audit-1101-lta | ok |
| `summary.latents_transition_enumeration` | R/lta-compare.R:613 | ran synthetic on a 2-candidate grid; test-audit-1101-lta | ok |
| `plot.latents_transition_enumeration` | R/lta-compare.R:621 | read; covered by existing enumeration plot tests | ok |
| `.lta_three_step_inputs` | R/lta-three-step.R:20 | ran synthetic: occasion posteriors vs hand computation w_k f_k(y) (max diff 3.3e-16, both engines); step-three scores vs central differences (6.8e-8); r3step with transition + initial covariates; test-lta-three-step | ok |
| `.lta_occasion_posteriors` | R/lta-three-step.R:103 | ran synthetic: occasion posteriors vs hand computation w_k f_k(y) (max diff 3.3e-16, both engines); step-three scores vs central differences (6.8e-8); r3step with transition + initial covariates; test-lta-three-step | ok |
| `.lta_three_step_basis` | R/lta-three-step.R:137 | ran synthetic: occasion posteriors vs hand computation w_k f_k(y) (max diff 3.3e-16, both engines); step-three scores vs central differences (6.8e-8); r3step with transition + initial covariates; test-lta-three-step | ok |
| `.lta_three_step_evaluate` | R/lta-three-step.R:171 | ran synthetic: occasion posteriors vs hand computation w_k f_k(y) (max diff 3.3e-16, both engines); step-three scores vs central differences (6.8e-8); r3step with transition + initial covariates; test-lta-three-step | ok |
| `.lta_three_step_fit` | R/lta-three-step.R:213 | ran synthetic: occasion posteriors vs hand computation w_k f_k(y) (max diff 3.3e-16, both engines); step-three scores vs central differences (6.8e-8); r3step with transition + initial covariates; test-lta-three-step | logged |
| `.lta_three_step` | R/lta-three-step.R:302 | ran synthetic: occasion posteriors vs hand computation w_k f_k(y) (max diff 3.3e-16, both engines); step-three scores vs central differences (6.8e-8); r3step with transition + initial covariates; test-lta-three-step | ok |
| `.multilpa_require_time` | R/sequences.R:6 | ran synthetic: get_results(fit, "sequences") long/wide and "sequence_summary" on a multilpa fit with time; test-sequences, test-tidy-sequences | ok |
| `.multilpa_time_name` | R/sequences.R:27 | ran synthetic: get_results(fit, "sequences") long/wide and "sequence_summary" on a multilpa fit with time; test-sequences, test-tidy-sequences | ok |
| `.multilpa_sequence_matrix` | R/sequences.R:46 | ran synthetic: get_results(fit, "sequences") long/wide and "sequence_summary" on a multilpa fit with time; test-sequences, test-tidy-sequences | ok |
| `.multilpa_occasion_names` | R/sequences.R:67 | ran synthetic: get_results(fit, "sequences") long/wide and "sequence_summary" on a multilpa fit with time; test-sequences, test-tidy-sequences | ok |
| `.multilpa_sequences` | R/sequences.R:81 | ran synthetic: get_results(fit, "sequences") long/wide and "sequence_summary" on a multilpa fit with time; test-sequences, test-tidy-sequences | ok |
| `.multilpa_sequence_summary` | R/sequences.R:127 | ran synthetic: get_results(fit, "sequences") long/wide and "sequence_summary" on a multilpa fit with time; test-sequences, test-tidy-sequences | ok |
| `get_tna` | R/tna.R:73 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_tna.default` | R/tna.R:79 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_tna.multilpa_transitions` | R/tna.R:85 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_group_tna` | R/tna.R:141 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_group_tna.default` | R/tna.R:147 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_group_tna.multilpa_transitions` | R/tna.R:153 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | fixed |
| `get_tna.multilpa_lta` | R/tna.R:174 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_group_tna.multilpa_lta` | R/tna.R:178 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `.lta_refuse_tna` | R/tna.R:184 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `.multilpa_warn_tna_empty_rows` | R/tna.R:196 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `.multilpa_require_tna` | R/tna.R:208 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `.multilpa_state_labels` | R/tna.R:222 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `.multilpa_aggregate_transitions` | R/tna.R:233 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `.multilpa_as_group_tna` | R/tna.R:273 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_tna.multilpa` | R/tna.R:288 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_group_tna.multilpa` | R/tna.R:299 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | fixed |
| `.multilpa_modal_sequences` | R/tna.R:332 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | logged |
| `get_tna.multilpa_covariates` | R/tna.R:295 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | ok |
| `get_group_tna.multilpa_covariates` | R/tna.R:311 | ran synthetic: get_tna()/get_group_tna() on homogeneous lta (2 classes) and multilpa(time =) fits, tna installed; test-tna, test-tna-fits | fixed |

## Fixed

1. **`.lta_bootstrap_dispatch()` — unclassed refusal of bootstrap controls.**
   `parameter_inference(fit, data = d, method = "bootstrap", iter = 1)` on an
   `lta()` fit (either engine; the homogeneous engine reaches it from
   R/transition-inference.R:355) gave
   `class: "simpleError" "error" "condition"`, where the profile families give
   `latents_bad_argument` for the same controls (test-release-inference.R:197).
   Cause: `stopifnot()` in the validator. Change: R/lta-bootstrap.R:365, the same
   named checks through `.latents_require()` (same messages, now
   `latents_bad_argument`). Test: test-audit-1101-lta.R "bootstrap controls of a
   transition fit are refused with a class" (iter = 1, n_starts = 0, tol = -1).
   Sabotage (restoring `stopifnot`) made it fail.
2. **`get_group_tna()` — unclassed refusal of `label`** (methods for
   `multilpa_transitions`, `multilpa`, and the `multilpa_covariates` alias).
   `get_group_tna(fit, label = 3)` gave `simpleError`. Change: R/tna.R:155 and
   R/tna.R:301, `.latents_require()` with the same message (existing tests that
   match "single string" still pass). Test: "get_group_tna() refuses a bad label
   with a class". Sabotage made it fail.
3. **`get_results()` on a transition enumeration — unclassed refusal.**
   `get_results(grid, "bset")` and `get_results(grid, criterion = "bci")` gave
   `simpleError` from `match.arg()`; every fit's table lookup is classed
   (TRANS-5). Change: R/lta-compare.R:573 and :575, `.latents_match_formal()`
   (same partial matching and defaults, `latents_bad_argument`, with a
   did-you-mean hint). Test: "a transition enumeration refuses unknown tables
   and bad counts with a class". Sabotage made it fail.
4. **`.lta_enumerate()` — crash on a malformed `n_profiles` / `n_group_classes`.**
   `enumerate_classes(d, vars, "id", time = "time", n_profiles = c(2, NA))` gave
   `simpleError: missing value where TRUE/FALSE needed`; `n_profiles = "a"` gave
   `non-numeric argument to binary operator`; `integer(0)` gave
   `'names' attribute [1] must be the same length as the vector [0]`.
   Cause: `any(n_profiles < 2)` evaluated before any type check. Change:
   R/lta-compare.R:440-445, a `.latents_require()` contract (non-empty finite
   whole numbers; `n_group_classes` ≥ 1) before the existing
   `latents_bad_transition` check for fewer than two profiles. Test: same block
   (c(2, NA), "2", integer(0), n_group_classes = 0). Sabotage (guard disabled)
   made it fail.

(The table counts 6 fixed functions: the four fixes touch
`.lta_bootstrap_dispatch`, `get_group_tna.multilpa_transitions`,
`get_group_tna.multilpa` (+ its `multilpa_covariates` alias),
`get_results.latents_transition_enumeration` and `.lta_enumerate`.)

The forward-backward test in the same file (path-sum likelihood with
covariates, two group classes, and shuffled rows) was sabotaged by reading
transition covariates at the origin instead of the destination occasion
(R/lta-engine.R `.lta_designs`): both expectations failed; restored.

## Logged, not fixed

- `.lta_nested()` (R/lta-compare.R:197) checks transition/initial covariates by
  name only; the bootstrap path is protected by `.lta_check_data()`
  (refuses with `latents_bad_inference_data`, see above). A refusal there under
  `latents_bad_nesting` would be more precise but needs no code to be correct.
- `.lta_bootstrap_lrt()` (R/lta-compare.R:349) refuses a non-nested pair with
  `latents_bad_argument`, where `compare_models(test = "lrt")` uses
  `latents_bad_nesting`. Existing tests pin `latents_bad_argument`
  (test-lta-compare.R:74-106); left as is (frozen API). Adding
  `latents_bad_nesting` as a second class would be compatible if the owner wants it.
- `.lta_three_step_fit()` (R/lta-three-step.R:213): a 3-profile step one on
  2-profile data (one profile with D row 0.65/0.16/0.19) makes step three stop
  at optim's 2000 iterations and refuse with `latents_no_converge` (largest score
  0.00121). Classed refusal of an over-extracted model; not a wrong number.
- `.multilpa_modal_sequences()` (R/tna.R:332): a fit with no classified row
  would hit `max(lengths(list()))` (-Inf, warning) before tna. Junk input only.

## Needs another slice

- R/additive-compare.R:286 and :288 — `get_results.latents_family_enumeration()`
  uses `match.arg()` for `what` and `criterion`, so an unknown table name is an
  unclassed `simpleError` (same defect as fix 3; `.latents_match_formal()` is the
  classed equivalent).

## Test runs

Run: `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`, `testthat::test_file()` per file after `devtools::load_all()`:

```
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 ERROR 0
test-audit-1101-lta.R: PASS 20 FAIL 0 SKIP 0 ERROR 0
test-bootstrap-lrt-extended.R: PASS 24 FAIL 0 SKIP 0 ERROR 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 ERROR 0
test-lta-acceleration.R: PASS 39 FAIL 0 SKIP 0 ERROR 0
test-lta-bootstrap.R: PASS 102 FAIL 0 SKIP 0 ERROR 0
test-lta-compare.R: PASS 31 FAIL 0 SKIP 0 ERROR 0
test-lta-covariance-boundary.R: PASS 52 FAIL 0 SKIP 0 ERROR 0
test-lta-general.R: PASS 35 FAIL 0 SKIP 0 ERROR 0
test-lta-parallel.R: PASS 99 FAIL 0 SKIP 2 ERROR 0
test-lta-thin-profiles.R: PASS 52 FAIL 0 SKIP 0 ERROR 0
test-lta-three-step.R: PASS 16 FAIL 0 SKIP 0 ERROR 0
test-release-audit-lta.R: PASS 41 FAIL 0 SKIP 0 ERROR 0
test-release-transitions.R: PASS 83 FAIL 0 SKIP 0 ERROR 0
test-review-1015-transitions.R: PASS 38 FAIL 0 SKIP 0 ERROR 0
test-review-1026-lta-open.R: PASS 45 FAIL 0 SKIP 0 ERROR 0
test-review-1026-lta-r2.R: PASS 47 FAIL 0 SKIP 0 ERROR 0
test-review-1026-lta.R: PASS 11 FAIL 0 SKIP 0 ERROR 0
test-review-1027-error-contract.R: PASS 216 FAIL 0 SKIP 0 ERROR 0
test-review-lta.R: PASS 305 FAIL 0 SKIP 0 ERROR 0
test-sequences.R: PASS 52 FAIL 0 SKIP 0 ERROR 0
test-tidy-sequences.R: PASS 26 FAIL 0 SKIP 0 ERROR 0
test-tna-fits.R: PASS 11 FAIL 0 SKIP 0 ERROR 0
test-tna.R: PASS 62 FAIL 0 SKIP 0 ERROR 0
test-transition-inference.R: PASS 23 FAIL 0 SKIP 0 ERROR 0
test-transitions-numerics.R: PASS 46 FAIL 0 SKIP 0 ERROR 0
test-transitions.R: PASS 118 FAIL 0 SKIP 0 ERROR 0
```

Engine golden fingerprint (`Rscript equivalence/engine-golden/compare.R`, repo root):

```
identical: 247; within_tolerance: 0; DIFFERENT: 13; ERROR: 0
```

All 13 differences are bootstrap tables (`profiles_two_level_VEI`,
`profiles_lca_single`, `profiles_count_two_level`, `lta_homogeneous_bootstrap`,
`lta_general_bootstrap`, `lta_general_weights_bootstrap`, seven
`svc_inference_*`), and each first difference is the attribute list only:
`replicate_rows` in the baseline vs `replicate_identity` now. That attribute is
set in R/bootstrap-inference.R:281, which was already modified in the working
tree before this audit (not a file of this slice). None of this slice's edits
touches a bootstrap attribute or any estimate.
