# Audit 2026-11-01 — slice `transitions-three-step`

Files: R/transitions.R, R/transition-inference.R, R/three-step.R,
R/three-step-categorical.R, R/three-step-ml.R, R/align-labels.R,
R/alignment-service.R. 110 functions checked; 5 defects fixed (7 functions marked fixed), 9 items logged,
2 need another slice. No reported number changed for valid input; no engine,
default, signature, table layout or roxygen changed (devtools::document() not
needed). The F7/F8/F10 changes made earlier today in three-step.R,
three-step-categorical.R and align-labels.R were left as they are.

Scratch evidence scripts (not part of the package):
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/tts/`
(`ts1.R`, `ts2.R`, `so1.R`, `so2.R`, `tr1.R`, `tr2.R`, `align1.R`, `contract1.R`,
`lta_contract.R`, `paths.R`, `evidence.R`).

## Independent checks (no defect found)

- **BCH weights = D^-1.** On a 2-class multilpa fit, D computed by hand as
  `sum_i p_it 1(W_i = s) / sum_i p_it` equals `get_results(fit, "classification_errors")`
  (max diff 0) and `get_results(fit, "bch_weights")` equals `solve(D)[W, ]` (max diff 0).
- **Weighted means.** BCH means = closed-form weighted least squares
  `solve(X'WX, X'Wy)` on the unit-by-class long format (diff 4.4e-16; `lm()` refuses
  negative weights); proportional = `lm(y ~ 0 + class, weights = posterior)`
  (2.695926, 0.3772016 both); modal = `tapply()`.
- **Sandwich.** CR1 cluster SEs by hand `G/(G-1) * crossprod(rowsum(w (y - mu)/T, g))`
  = package (0.1400523, 0.1188245); independent n/(n-1) the same; pair difference
  -3.176058 with SE `sqrt(V11 + V22 - 2 V12)` = 0.1895952 both; overall Wald 280.6222 both.
- **ML numeric step three.** An independent likelihood
  `sum log sum_t pi_t D[t, W_i] dnorm(y; mu_t, s_t)` with `pi = colMeans(posteriors)`
  maximised by BFGS from the package estimate moves nothing (3.068059, -0.03313691;
  sd 2.046513/2.046514, 1.070726); equal variance likewise. Sandwich from a numDeriv
  Hessian of independent scores: SEs 0.1271254, 0.08286887 and pair SE 0.1427282, identical.
- **Categorical.** BCH probabilities = `crossprod(solve(D)[W,], Y) / colSums` (diff 0);
  ML log likelihood at the estimate equals an independent BFGS optimum (parameter
  shift 1.5e-14); SEs from numDeriv Hessian + softmax Jacobian agree to 3.8e-12.
- **BH families.** three_step pairs: `p_adjusted == p.adjust(p_value, "BH")` over all
  rows (categorical 9 rows); r3step: BH over covariate rows only, intercept rows NA.
- **Step-one correction.** The closure's base point reproduces the stored posteriors
  (< 1.5e-15) for 13 fit kinds (FIML, mixed, ordinal + count, VVV/EEE/VEI, local
  dependence, direct effects, two group classes at both levels, person centring and
  scaling). The SE^2 increment equals `diag(G V1 G')` with G from numDeriv
  (3.177765e-04, 3.545946e-04 both; pair 8.42817e-05 both). It is >= 0 and vanishes
  with separation: BCH SE fixed/estimated 0.104/0.186 (sep 1.2), 0.0806/0.0847 (2),
  0.0672/0.0673 (3), 0.0673/0.0673 (5); ML the same pattern.
- **Transitions.** Rows of every transition table sum to 1. Re-evaluating the decoded
  likelihood reproduces `fit$log_likelihood` (gap 0 to 2.3e-13). `vcov(scale =
  "unconstrained")` vs the inverse of a numDeriv Hessian of that likelihood: max SE
  ratio deviation 9.9e-9 (VVI), 6.0e-8 (VVV), 1.2e-8 (EEE), 7.3e-9 (mixed), 1.1e-8
  (grid), unbalanced data. A numDeriv Jacobian of an independent softmax (reference
  last) times that covariance reproduces every initial/transition probability SE
  exactly (observed 0.04216585, 0.0254102, 0.03792118; robust 0.04247928, 0.02519267,
  0.03933089).
- **align_labels().** Beyond the existing permuted-start refits (lpa VVV, lca,
  two-level multilpa), fits with ordinal + negative-binomial counts, continuous and
  categorical local dependence, direct effects, a noise class and a two-level
  covariate fit were permuted (3,1,2): every get_results() table matched the original
  up to relabelling (ordinal locations differ only by the designed rebasing; the
  adjacent-category intercept shift `+ c * reference` was derived and checked), and
  parameter_inference() accepted each permuted fit's likelihood.

## Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.multilpa_sequence_layout` | R/transitions.R:27 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_sequence_emission` | R/transitions.R:70 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_forward_backward` | R/transitions.R:107 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_sequence_moments` | R/transitions.R:133 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_transition_validate` | R/transitions.R:190 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_transition_prevalence` | R/transitions.R:220 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_transition_expectation` | R/transitions.R:247 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_transition_maximization` | R/transitions.R:311 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_transition_em` | R/transitions.R:376 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.multilpa_transition_initialize` | R/transitions.R:430 | ran synthetic: lta() on simulated 2-state data (balanced/unbalanced, grid, VVI/VVV/EEE, categorical); log likelihood re-evaluated through decode+expectation (gap 0 to 2e-13); covered by test-transitions.R, test-transitions-numerics.R | ok |
| `.lta_check_time_order` | R/transitions.R:475 | ran synthetic (text time "1".."15" -> latents_bad_time); covered by test-transitions.R, test-review-1015-transitions.R | ok |
| `.multilpa_count_transition_parameters` | R/transitions.R:514 | read (count = cross-sectional + G K(K-1)); length(theta) == n_parameters asserted in every vcov run; test-transitions.R | ok |
| `lta` | R/transitions.R:748 | ran synthetic: argument contracts (17 bad inputs), course_engagement fits; test-audit-1101-transitions-three-step.R | fixed |
| `.lta_setup` | R/transitions.R:788 | ran synthetic (as lta()); unclassed stopifnot replaced | fixed |
| `.lta_run_start` | R/transitions.R:983 | ran synthetic via lta(); fields read (counts, BIC, empty rows, small_classes); test-transitions.R | ok (one item logged) |
| `.lta_finalise` | R/transitions.R:998 | ran synthetic via lta(); fields read (counts, BIC, empty rows, small_classes); test-transitions.R | ok (one item logged) |
| `.multilpa_transitions_table` | R/transitions.R:1176 | ran synthetic (estimated = NA, stable = c(T,F)) | fixed |
| `.multilpa_optional_flag` | R/transitions.R:1194 | ran synthetic: rows per origin sum to 1 (all structures), stable/estimated restriction; test-transitions.R | ok |
| `.multilpa_transition_frame` | R/transitions.R:1207 | ran synthetic: rows per origin sum to 1 (all structures), stable/estimated restriction; test-transitions.R | ok |
| `.multilpa_restrict_transitions` | R/transitions.R:1230 | ran synthetic: rows per origin sum to 1 (all structures), stable/estimated restriction; test-transitions.R | ok |
| `.multilpa_initial_frame` | R/transitions.R:1243 | ran synthetic: rows per origin sum to 1 (all structures), stable/estimated restriction; test-transitions.R | ok |
| `as.data.frame.multilpa_transitions` | R/transitions.R:1258 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `print.multilpa_transitions` | R/transitions.R:1265 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `logLik.multilpa_transitions` | R/transitions.R:1285 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `nobs.multilpa_transitions` | R/transitions.R:1293 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `summary.multilpa_transitions` | R/transitions.R:1300 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `as.data.frame.summary_multilpa_transitions` | R/transitions.R:1343 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `print.summary_multilpa_transitions` | R/transitions.R:1350 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `coef.multilpa_transitions` | R/transitions.R:1379 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `plot.multilpa_transitions` | R/transitions.R:1391 | ran synthetic (print/summary/logLik/nobs/coef/plot on course_engagement fit); @noRd per owner rule; test-transitions.R | ok |
| `.multilpa_transition_view` | R/transition-inference.R:28 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `.multilpa_transition_shared_blocks` | R/transition-inference.R:52 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `.multilpa_transition_encode` | R/transition-inference.R:73 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `.multilpa_transition_decode` | R/transition-inference.R:107 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `.multilpa_transition_group_scores` | R/transition-inference.R:149 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `.multilpa_transition_group_counts` | R/transition-inference.R:221 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `.multilpa_transition_jacobian` | R/transition-inference.R:265 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `.multilpa_transition_labels` | R/transition-inference.R:300 | ran synthetic: numDeriv Hessian of the decoded likelihood vs vcov(scale = "unconstrained") (SE ratio within 6e-8, 5 structures); numDeriv softmax Jacobian x V reproduces every initial/transition probability SE exactly; test-transition-inference.R | ok |
| `parameter_inference.multilpa_transitions` | R/transition-inference.R:335 | ran synthetic (observed and robust SEs = J V J' by hand; level = 2) | fixed |
| `.multilpa_transition_covariance` | R/transition-inference.R:404 | ran synthetic (step = 0, -1e-4, "a" via vcov()) | fixed |
| `.multilpa_check_transition_regularity` | R/transition-inference.R:487 | ran synthetic (boundary 2-class fit refused latents_boundary_fit; wrong data latents_bad_inference_data); test-transition-inference.R | ok |
| `.multilpa_check_transition_data` | R/transition-inference.R:536 | ran synthetic (boundary 2-class fit refused latents_boundary_fit; wrong data latents_bad_inference_data); test-transition-inference.R | ok |
| `vcov.multilpa_transitions` | R/transition-inference.R:562 | ran synthetic (type = warns latents_unused_argument; step now checked in the shared covariance) | fixed |
| `confint.multilpa_transitions` | R/transition-inference.R:587 | ran synthetic (parm name/position, bad parm, level = 1.5) | ok |
| `.multilpa_level_assignments` | R/three-step.R:20 | ran synthetic: D by hand from posteriors and modal = get_results(fit, "classification_errors") (max diff 0); empty-row FIML fit (1362 of 1422 kept); test-three-step.R | ok |
| `.multilpa_check_three_step_fit` | R/three-step.R:56 | ran synthetic (data.frame, list -> latents_not_a_fit); test-three-step.R, test-lta-three-step.R | ok |
| `.multilpa_classification_errors` | R/three-step.R:126 | ran synthetic: BCH table = solve(D)[modal, ] exactly (2-class toy); test-three-step.R | ok |
| `.multilpa_error_frame` | R/three-step.R:141 | ran synthetic: D by hand from posteriors and modal = get_results(fit, "classification_errors") (max diff 0); empty-row FIML fit (1362 of 1422 kept); test-three-step.R | ok |
| `.multilpa_error_matrix` | R/three-step.R:164 | ran synthetic: D by hand from posteriors and modal = get_results(fit, "classification_errors") (max diff 0); empty-row FIML fit (1362 of 1422 kept); test-three-step.R | ok |
| `.multilpa_bch_weights` | R/three-step.R:185 | ran synthetic: BCH table = solve(D)[modal, ] exactly (2-class toy); test-three-step.R | ok |
| `.multilpa_weight_frame` | R/three-step.R:199 | ran synthetic: BCH table = solve(D)[modal, ] exactly (2-class toy); test-three-step.R | ok |
| `.multilpa_invert_errors` | R/three-step.R:219 | ran synthetic: BCH table = solve(D)[modal, ] exactly (2-class toy); test-three-step.R | ok |
| `three_step` | R/three-step.R:458 | ran synthetic: BCH means = closed-form WLS on the long format (4e-16), CR1 and independent SEs by hand exact, proportional = lm(weights = posteriors), modal = tapply, overall Wald by hand; lta/FIML/group-level paths; 20 contract cases | ok |
| `.multilpa_step_one_correction` | R/three-step.R:665 | ran synthetic: base point reproduces stored posteriors (<1.5e-15) for 13 fit kinds; SE^2 increment = diag(G V1 G') with numDeriv G (exact to print); increment >= 0 and -> 0 as separation grows (sep 1.2/2/3/5) | ok |
| `.multilpa_covariance_root` | R/three-step.R:725 | ran synthetic: base point reproduces stored posteriors (<1.5e-15) for 13 fit kinds; SE^2 increment = diag(G V1 G') with numDeriv G (exact to print); increment >= 0 and -> 0 as separation grows (sep 1.2/2/3/5) | ok |
| `.multilpa_step_clusters` | R/three-step.R:748 | ran synthetic (CR1 by hand); test-three-step-clusters.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_thin_classes` | R/three-step.R:787 | ran synthetic (CR1 by hand); test-three-step-clusters.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_thin` | R/three-step.R:799 | ran synthetic (CR1 by hand); test-three-step-clusters.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_pairs` | R/three-step.R:831 | ran synthetic: difference and SE sqrt(V11+V22-2V12) by hand exact (BCH and ML); BH over valid pairs | ok |
| `.multilpa_outcome_values` | R/three-step.R:874 | ran synthetic (group-level outcome, non-constant refused latents_bad_outcome) | ok |
| `.multilpa_step_weights` | R/three-step.R:891 | ran synthetic (bch/proportional/modal vs hand) | ok |
| `.multilpa_logit_prior` | R/three-step.R:908 | read; exercised by r3step runs | ok |
| `r3step` | R/three-step.R:1096 | ran synthetic: BH family = covariate terms only (p.adjust by hand equal, intercepts NA); lta, groups, by_group_class robust/observed, transition covariates; contract cases | ok |
| `.r3step_require_classes` | R/three-step.R:1279 | ran synthetic (one-class level, missing column, initial covariates on lpa) | ok |
| `.r3step_check_columns` | R/three-step.R:1294 | ran synthetic (one-class level, missing column, initial covariates on lpa) | ok |
| `.r3step_transitions` | R/three-step.R:1306 | ran synthetic (one-class level, missing column, initial covariates on lpa) | ok |
| `.multilpa_r3step_design` | R/three-step.R:1343 | ran synthetic via r3step(); test-three-step.R | ok |
| `.multilpa_r3step_check_rank` | R/three-step.R:1365 | ran synthetic via r3step(); test-three-step.R | ok |
| `.multilpa_r3step_standard_errors` | R/three-step.R:1385 | ran synthetic via r3step(); test-three-step.R | ok |
| `.multilpa_r3step_frame` | R/three-step.R:1403 | ran synthetic via r3step(); test-three-step.R | ok |
| `.multilpa_r3step_by_group_class` | R/three-step.R:1456 | ran synthetic (FIML 2-group-class multilpa, robust and observed); test-three-step.R | ok |
| `.multilpa_r3step_group_model` | R/three-step.R:1500 | ran synthetic (FIML 2-group-class multilpa, robust and observed); test-three-step.R | ok |
| `.multilpa_r3step_group_frame` | R/three-step.R:1597 | ran synthetic (FIML 2-group-class multilpa, robust and observed); test-three-step.R | ok |
| `.multilpa_measured_items` | R/three-step.R:1652 | ran synthetic (indicator as outcome -> multilpa_indicator_reused) | ok |
| `.multilpa_three_step_categorical` | R/three-step-categorical.R:26 | ran synthetic: BCH probabilities = crossprod(solve(D)[W,], Y)/colSums exactly; pairs BH = p.adjust over rows; overall on logits; FIML/lta paths; F7/F10 tests | ok |
| `.multilpa_step_categorical_fit` | R/three-step-categorical.R:167 | ran synthetic: log likelihood at estimate = BFGS optimum of an independent likelihood (par shift 1.5e-14); sandwich SEs from numDeriv Hessian + softmax Jacobian match (3.8e-12) | ok |
| `.multilpa_step_categorical_derivative` | R/three-step-categorical.R:225 | ran synthetic: log likelihood at estimate = BFGS optimum of an independent likelihood (par shift 1.5e-14); sandwich SEs from numDeriv Hessian + softmax Jacobian match (3.8e-12) | ok |
| `.multilpa_step_categorical_scores` | R/three-step-categorical.R:241 | ran synthetic: log likelihood at estimate = BFGS optimum of an independent likelihood (par shift 1.5e-14); sandwich SEs from numDeriv Hessian + softmax Jacobian match (3.8e-12) | ok |
| `.multilpa_step_categorical_information` | R/three-step-categorical.R:258 | ran synthetic: log likelihood at estimate = BFGS optimum of an independent likelihood (par shift 1.5e-14); sandwich SEs from numDeriv Hessian + softmax Jacobian match (3.8e-12) | ok |
| `.latents_block_diagonal` | R/three-step-categorical.R:277 | read; exercised by categorical ML information | ok |
| `.multilpa_step_categorical_influence` | R/three-step-categorical.R:296 | ran synthetic: log likelihood at estimate = BFGS optimum of an independent likelihood (par shift 1.5e-14); sandwich SEs from numDeriv Hessian + softmax Jacobian match (3.8e-12) | ok |
| `.multilpa_step_categorical_pairs` | R/three-step-categorical.R:330 | ran synthetic (9 rows, BH/holm by hand); test-review-1101-three-step-categorical.R | ok |
| `.multilpa_step_overall` | R/three-step-categorical.R:388 | ran synthetic: Wald = d^2 / var(d) by hand (2 classes); categorical df (K-1)(C-1) | ok |
| `.multilpa_step_ml_fit` | R/three-step-ml.R:26 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_ml_polish` | R/three-step-ml.R:261 | ran synthetic (mocked singular and malformed Hessian) | fixed |
| `.multilpa_step_ml_responsibilities` | R/three-step-ml.R:311 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_ml_scores` | R/three-step-ml.R:330 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_ml_coordinates` | R/three-step-ml.R:353 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_ml_scale` | R/three-step-ml.R:369 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_ml_hessian` | R/three-step-ml.R:388 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_ml_log_likelihood` | R/three-step-ml.R:403 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_step_ml_influence` | R/three-step-ml.R:425 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `.multilpa_three_step_ml` | R/three-step-ml.R:465 | ran synthetic: estimates = BFGS optimum of an independent likelihood with pi = colMeans(posteriors), D fixed (varying and equal); sandwich SEs (numDeriv) and pair SE exact; test-three-step-ml.R, test-features-ml-three-step.R | ok |
| `align_labels` | R/align-labels.R:48 | ran synthetic: refits from permuted starts (existing tests) plus permuted ordinal+NB count, continuous and categorical local dependence, direct effects, noise and 2-level covariate fits: every get_results() table equal up to relabelling (ordinal locations rebased by design) | ok |
| `.multilpa_permutations` | R/alignment-service.R:29 | read + ran via align probes; test-alignment-service.R | ok |
| `.multilpa_match_order` | R/alignment-service.R:56 | read + ran via align probes; test-alignment-service.R | ok |
| `.latents_measurement_signature` | R/alignment-service.R:126 | ran via align probes; test-alignment-service.R | ok (one item logged) |
| `.multilpa_extra_signature` | R/alignment-service.R:167 | ran via align probes; test-alignment-service.R | ok (one item logged) |
| `.multilpa_profile_signature` | R/alignment-service.R:209 | ran via align probes; test-alignment-service.R | ok (one item logged) |
| `.multilpa_cov_profile_signature` | R/alignment-service.R:223 | ran via align probes; test-alignment-service.R | ok (one item logged) |
| `.multilpa_group_signature` | R/alignment-service.R:236 | ran via align probes; test-alignment-service.R | ok (one item logged) |
| `.multilpa_cov_group_signature` | R/alignment-service.R:248 | ran via align probes; test-alignment-service.R | ok (one item logged) |
| `.latents_rebase_logits` | R/alignment-service.R:275 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |
| `.multilpa_permute_extra` | R/alignment-service.R:289 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |
| `.latents_permute_gaussian_diagnostics` | R/alignment-service.R:323 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |
| `.multilpa_permute_profiles` | R/alignment-service.R:340 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |
| `.multilpa_permute_group_classes` | R/alignment-service.R:396 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |
| `.multilpa_cov_permute` | R/alignment-service.R:434 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |
| `.multilpa_pool_align` | R/alignment-service.R:516 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |
| `.multilpa_align_labels` | R/alignment-service.R:555 | ran synthetic: permuted fits reproduce the likelihood (parameter_inference accepted) and every table up to relabelling; adjacent-category intercept rebasing derived and checked; test-align-labels.R, test-alignment-service.R | ok |

## Fixed

**1. lta() raised unclassed errors for argument mistakes multilpa() classes**
(R/transitions.R:834-848, `.lta_setup`; also R/transitions.R:823 `mover_stayer`).
Reproduction before the fix (`evidence.R`, old code sourced from HEAD):
`lta(course_engagement, v, NULL, 2, "sequence")` -> `simpleError | is.character(id) is not TRUE`;
id omitted -> `missingArgError | argument "id" is missing, with no default`;
time omitted -> `missingArgError | argument "time" is missing, with no default`;
`vars = 1:2` -> `simpleError | is.character(vars) is not TRUE`;
`mover_stayer = "yes"`, `min_probability = 2`, `categorical = NA` -> `simpleError`.
multilpa() raises `latents_bad_argument` for the same mistakes. Cause: a bare
`stopifnot()` in `.lta_setup()` (shared by both transition engines). Change: the
checks go through `.latents_require()` (class `latents_bad_argument`), with
`!missing(id)` / `!missing(time)` so an omitted required argument gets the same
classed message; `id` must be a single non-missing string. After: all of the above
raise `latents_bad_argument`. Guard: test-audit-1101-transitions-three-step.R,
"lta() classes its argument mistakes as multilpa() does".

**2. parameter_inference() on an lta() fit left `level`/`step` mistakes unclassed**
(R/transition-inference.R:345-349). Before: `parameter_inference(fit, level = 2)` ->
`simpleError`, while the profile-fit method classes it (`latents_bad_argument`,
review 2026-10-27 INF-6; test-release-inference.R pins it there). Change:
`.latents_require()` for `level`. Guard: same test file, "the transition table's
flags and the Wald arguments are classed".

**3. vcov() on an lta() fit took any `step`** (R/transition-inference.R:404-409,
`.multilpa_transition_covariance`). Before: `vcov(fit, step = 0)` ->
`simpleError | infinite or missing values in 'x'` (division by zero in the
Hessian); `vcov(fit, step = -1e-4)` silently returned the default covariance;
`parameter_inference()` refused both. Change: the `step` check moved into the
covariance helper both methods call, classed. Guard: same test (step 0, -1e-4,
"a", and vcov() still equals parameter_inference()'s covariance).

**4. get_results(fit, "transitions", estimated = NA) raised an unclassed error**
(R/transitions.R:1178-1185, `.multilpa_transitions_table`). Before:
`simpleError | estimated must be NULL, TRUE or FALSE`. Change: `.latents_require()`.
Guard: same test.

**5. The ML step-three Newton polish swallowed every error of its solve**
(R/three-step-ml.R:276-287, `.multilpa_step_ml_polish`). `tryCatch(solve(...),
error = function(condition) NULL)` treated any error (including a mistake in the
call) as "no Newton direction" and returned the EM estimate silently. Change: the
one intended rejection is tested explicitly (non-finite or `rcond() <
.Machine$double.eps`, solve()'s own tolerance), and solve() is then called
unguarded. Results are unchanged on every path (all ML tests and both Latent GOLD
kits pass). Guard: same test file, mocked singular Hessian (estimate returned as
it came) and malformed Hessian (error surfaces; with the old code it was swallowed).

Each guard was broken on purpose by restoring the HEAD versions of the three files:
tests 1 and 3 errored (wrong condition class) and test 4 failed its malformed-Hessian
expectation; the files were then restored and the tests pass.

## Logged, not fixed

- `three_step()`, `r3step()`, `lta()` use plain `match.arg()` for `method`, `level`,
  `contrast`, `adjust`, `vcov_type` etc., so an unknown value raises base R's
  unclassed error (`three_step(fit, d, "y", method = "foo")`). Package-wide pattern
  (dozens of verbs); not changed in one slice.
- Categorical `contrast = "pairs"` puts every pair x category in the BH family, so a
  binary outcome has each pair twice (differences d and -d). BH is unaffected (the
  adjusted values equal BH over the distinct p-values), but Holm/Bonferroni count
  each test twice. The family is as documented (one row per pair and category).
- Categorical step three requires G > K independent groups, but its covariance has
  K(C-1) parameters; with K < G <= K(C-1) some category contrasts get a singular,
  near-zero variance. Only with a handful of clusters.
- `.lta_finalise()` sets `small_classes` from all rows (including wholly empty ones
  under `drop_empty = TRUE`), while `effective_profile_counts` uses the counted rows;
  the general engine does the same. Differs only for a near-empty profile with empty rows.
- `print.summary_multilpa_transitions()` says "at least one variance is at
  min_variance" also for a full-covariance fit, where the boundary is an eigenvalue.
- `.multilpa_transition_expectation()` / `.multilpa_transition_maximization()` raise
  unclassed `stop()` on non-finite counts; caught per start, can only reach the user
  from vcov() at a displaced point of a pathological fit.
- The profile mixing column in the alignment signature is the unweighted
  `colMeans(profile_probabilities)` over group classes, not the group-weighted
  marginal; it is only a tie-breaker.
- `.multilpa_step_ml_fit()` evaluates `if (stats::var(values) <= 0)`, NA for a
  single unit (not a realistic input).
- `r3step(by_group_class = TRUE)` on a one-group-class fit fits the pooled model
  under group-class labels rather than refusing; the docs do not say either way.

## Needs another slice

- R/covariate-inference.R:394 (`.multilpa_cov_covariance`): `step` is checked by a
  bare `stopifnot()`; `parameter_inference(covariate_fit, step = -1)` raises
  `simpleError` (verified), while INF-6 classes it `latents_bad_argument` for
  plain profile fits.
- R/additive-inference.R:189 (`.additive_inference`): the same bare `stopifnot()`
  for `step` (read only; whether a classed check upstream catches it first was not
  verified).

## Equivalence and golden checks

- `latentgold-categorical-distal`: 6 expectations passed, 0 failed (tmp/equivalence-kits/latentgold-categorical-distal/test-results.csv).
- `latentgold-step-one`: 8 expectations passed, 0 failed (tmp/equivalence-kits/latentgold-step-one/test-results.csv).
  Both runs then stopped in `matrix_record_evidence()` with "runtime source changed
  during the comparison": other slices were editing R/ concurrently, so the
  provenance record was refused; the comparisons themselves completed.
- `Rscript equivalence/engine-golden/compare.R` (full): 242 identical, 18 DIFFERENT.
  13 are the known bootstrap attribute rename (`replicate_rows` -> `replicate_identity`).
  Mine: `lta_refusal_no_time` (`simpleError` -> `latents_bad_argument`, fix 1, intended).
  Not this slice's files: `svc_fit_staged_one_class` (`simpleError` ->
  `latents_bad_argument`, fit_staged), `cov_weights`, `family_additive_weights`,
  `reg_gaussian_weights_observed_refused` (refusal message text changed).

## Test runs (verbatim)

```
test-audit-1101-transitions-three-step.R: PASS 25 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-align-labels.R: PASS 102 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-alignment-service.R: PASS 64 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-assignment-alignment.R: PASS 26 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-review-1101-align-ordinal.R: PASS 34 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-review-1101-three-step-categorical.R: PASS 727 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-three-step-categorical.R: PASS 35 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-three-step-clusters.R: PASS 41 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-three-step-ml.R: PASS 61 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-three-step-step-one.R: PASS 172 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-three-step.R: PASS 74 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-tidy-three-step.R: PASS 74 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-features-ml-three-step.R: PASS 78 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-lta-three-step.R: PASS 16 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-transition-inference.R: PASS 23 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-transitions-numerics.R: PASS 46 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-transitions.R: PASS 118 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-release-transitions.R: PASS 83 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-review-1015-transitions.R: PASS 38 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-growth-r3step.R: PASS 12 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 ERROR 0 WARN 0
```
(NOT_CRAN=true, one R process, testthat::test_file per file; the same files gave
identical counts before the fixes, minus the new file.)
