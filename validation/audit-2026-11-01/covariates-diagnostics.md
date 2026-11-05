# Audit 2026-11-01 — slice `covariates-diagnostics`

Files: R/covariates.R, R/covariate-inference.R, R/direct-effects.R,
R/cross-level.R, R/local-dependence.R, R/bivariate-residuals.R,
R/gaussian-diagnostics.R, R/diagnostics.R, R/diagnostics-checks.R,
R/report.R. 160 functions (plus two constants in diagnostics-checks.R)
checked; 4 code defects fixed, 1 documentation defect fixed, 9 items logged,
2 items for another slice. No estimate, default, signature or table layout
changed for valid input; the one printed value that changes is the
restricted cross-level replication count, which was wrong.

Regression tests: `tests/testthat/test-audit-1101-covariates-diagnostics.R`.
Scratch evidence scripts (not part of the package):
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/covdiag/`
(`mstep.R` M-step vs nnet, `hess.R` SEs vs numDeriv, `thin.R` thin-profile
warnings by structure, `bvr.R` bivariate residuals by hand, `diag.R`/`fam.R`
diagnostics across 13 fits, `cl.R`/`clget.R`/`clarg.R` cross-level,
`direct.R` direct effects vs lm, `ld.R` local dependence, `gd.R` Gaussian
diagnostics by hand, `ic.R` criteria/entropy/classification by hand,
`center.R` centred and covariate residual scale, `stored.R` stored data and
covariate methods, `rep.R` report() across 13 fits, `sabotage.py` the
break-each-fix run).

## Independent checks (no defect found)

- **Membership-logit M-step** (`.multilpa_weighted_logits`) against
  `nnet::multinom` with posterior weights, 3 categories, 60 rows: max
  coefficient difference 2.8e-7 (row posteriors), 3.4e-8 (group-class joint
  counts, rows summing above 1), 9.5e-10 (5 zero-count rows); the residual is
  nnet's own stopping tolerance. In the package test (nnet is not in
  Suggests): two categories vs `glm(quasibinomial)` to 1e-8, three
  categories vs a direct `optim()` maximization to 1e-5 with the package's
  objective no worse.
- **Covariate SEs vs inverse numerical Hessian** (`numDeriv::hessian` of the
  log likelihood in the encode coordinates): shared slopes with a group
  covariate, relative differences ≤ 5.2e-8 on all 5 coefficients; group-class
  slopes ≤ 8.6e-9; all unconstrained coordinates (incl. log variances) ≤
  5.4e-6. Log likelihood reproduced exactly (difference 0).
- **Bivariate residuals by hand**: profile-conditional 2x2 Pearson
  chi-square 0.006849948 = package; pooled 0.2209987 = package; Cramér's V
  observed 0.006039675 = package; Gaussian Fisher z on posterior-weighted
  residuals r = 0.4254912, z = 5.585561 = package.
- **BH**: `diagnostics()` table `p_adjusted` = `p.adjust(p_value, "BH")`;
  the check's `value` = count of BH-adjusted p < 0.05; `get_results(fit,
  "residuals")` stays unadjusted by default (LEARNINGS 2026-10-31).
- **report()/diagnostics() claims vs tables**: smallest class, entropy,
  lowest AvePP, lowest OCC, BIC, parameter count, replication read from the
  tables on 13 fits (lpa, 1-profile lpa, lca, weighted, noise, multilpa,
  multilca, covariates, transitions, additive, restricted cross-level,
  regression, growth); `report()` runs with plots on all 9 profile-family
  fits and refuses the 4 other families with `latents_not_a_fit`.
- **Thresholds**: the only verdict threshold is `.latents_residual_alpha`
  (0.05, BH), documented on `?diagnostics` and printed as the criterion;
  "fail" for non-convergence and effective class size below one, "check"
  for failed/unfinished starts, bounds, empty classes and flagged pairs, as
  documented; every other check is `info`.
- **Weighted criteria and classification** (weights with 2 zero-weight
  groups): `information_criteria` wide = fit `aic`/`bic`/`bic_individual`;
  ICL by hand from the complete-data entropy; relative entropy by hand
  0.888985; `n_modal` excludes zero-weight rows (93/81/106 vs 97/89/114 over
  all rows).
- **Direct effects**: one-profile fit = `lm()` (slopes 0.3130264,
  0.2909154; means = lm intercepts at z = 0; SEs = ML-variance OLS SEs
  0.01094176, 0.1126164), covariate mean 40 (uncentring path).
- **Local dependence**: half-observed pair rows restore to the original
  items exactly; joint table probabilities sum to 1; continuous pair
  covariance/correlation = `covariances[, , k]`; VVI + a:b counts 15
  parameters (6 means + 6 variances + 2 covariances + 1).
- **Gaussian diagnostics** with FIML, VVV and 6 zero-weight rows: minimum
  standardized eigenvalue, Gaussian mass, variance ratio, Kish n and
  positive rows equal hand computations.
- **Centring**: residual tables and Gaussian diagnostics agree between plain
  and covariate fits under `centering = "none"/"grand"/"person"`.
- **MoEClust equivalence** (`equivalence/test-moeclust-covariate-structures.R`):
  see the test-run section.

## Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.multilpa_softmax` | R/covariates.R:6 | ran synthetic (mstep.R, via the M-step) | ok |
| `.multilpa_log_softmax` | R/covariates.R:15 | ran synthetic (mstep.R) | ok |
| `.multilpa_weighted_logits` | R/covariates.R:58 | ran vs nnet::multinom (mstep.R); test audit-1101 (glm, optim) | ok |
| `.multilpa_cov_expectation` | R/covariates.R:193 | ran (hess.R: log likelihood reproduced, diff 0) | ok |
| `.multilpa_fit_covariates` | R/covariates.R:238 | ran (thin.R, hess.R, ic.R); MoEClust kit | fixed (F1) |
| `print.multilpa_covariates` | R/covariates.R:576 | ran (stored.R) | ok |
| `summary.multilpa_covariates` | R/covariates.R:591 | ran (stored.R) | ok |
| `print.summary_multilpa_covariates` | R/covariates.R:617 | ran (stored.R) | ok |
| `as.data.frame.summary_multilpa_covariates` | R/covariates.R:646 | read; covered by test-covariates.R | ok |
| `.multilpa_cov_stored_data` | R/covariates.R:664 | ran (stored.R: equals data, shared and group-class slopes) | ok |
| `logLik.multilpa_covariates` | R/covariates.R:717 | ran (stored.R: df = n_parameters, AIC/BIC equal) | ok |
| `nobs.multilpa_covariates` | R/covariates.R:725 | ran (stored.R) | ok |
| `.multilpa_cov_check_covariates` | R/covariates.R:737 | read; covered by test-covariates.R | ok |
| `.multilpa_design_basis` | R/covariates.R:778 | read; ran through mstep.R/hess.R | ok |
| `.multilpa_cov_extreme` | R/covariates.R:817 | read; covered by test-covariate-separation.R | ok |
| `.multilpa_cov_designs` | R/covariates.R:860 | ran (stored.R, center.R, ic.R with zero weights) | ok |
| `.multilpa_cov_start` | R/covariates.R:939 | ran through every covariate fit | ok |
| `.multilpa_cov_initial` | R/covariates.R:952 | ran through every covariate fit; MoEClust kit | ok |
| `.multilpa_cov_run` | R/covariates.R:1007 | ran through every covariate fit; MoEClust kit | ok |
| `.multilpa_cov_assemble` | R/covariates.R:1066 | ran (thin.R parameter counts; stored.R) | ok |
| `.multilpa_cov_encode` | R/covariate-inference.R:12 | ran (hess.R) | ok |
| `.multilpa_cov_is_full` | R/covariate-inference.R:42 | read; covered by test-covariate-measurement.R | ok |
| `.multilpa_cov_chart_width` | R/covariate-inference.R:51 | read; covered by test-covariate-structure-inference.R | ok |
| `.multilpa_cov_spread_rows` | R/covariate-inference.R:60 | read; covered by test-covariate-measurement.R | ok |
| `.multilpa_cov_cholesky_coordinates` | R/covariate-inference.R:74 | read; covered by test-covariate-measurement.R | ok |
| `.multilpa_cov_cholesky_decode` | R/covariate-inference.R:91 | read; covered by test-covariate-measurement.R | ok |
| `.multilpa_cov_decode` | R/covariate-inference.R:109 | ran (hess.R) | ok |
| `.multilpa_cov_group_scores` | R/covariate-inference.R:176 | ran (hess.R: analytic-gradient Hessian = numDeriv) | ok |
| `.multilpa_cov_full_scores` | R/covariate-inference.R:275 | read; covered by test-covariate-measurement.R, test-covariate-missing.R | ok |
| `parameter_inference.multilpa_covariates` | R/covariate-inference.R:330 | ran (hess.R: SEs = numDeriv to 5e-8) | ok |
| `.multilpa_cov_covariance` | R/covariate-inference.R:383 | ran (hess.R) | ok |
| `.multilpa_cov_parameter_names` | R/covariate-inference.R:618 | ran (stored.R) | ok |
| `vcov.multilpa_covariates` | R/covariate-inference.R:629 | ran (hess.R unconstrained SEs to 5.4e-6); test-review-1101-confint-fixed.R | ok |
| `.multilpa_cov_drop_held` | R/covariate-inference.R:682 | read; covered by test-review-1101-confint-fixed.R | ok |
| `.multilpa_cov_estimation_labels` | R/covariate-inference.R:700 | ran (hess.R via vcov(scale = "unconstrained")) | ok |
| `.multilpa_cov_inference_frame` | R/covariate-inference.R:722 | ran (hess.R) | ok |
| `.multilpa_cov_held_parameters` | R/covariate-inference.R:783 | read; covered by test-review-1101-confint-fixed.R | ok |
| `.multilpa_cov_front_estimate` | R/covariate-inference.R:804 | ran (hess.R/stored.R via coef) | ok |
| `.multilpa_cov_front_jacobian` | R/covariate-inference.R:836 | ran (hess.R) | ok |
| `.multilpa_cov_chart_front` | R/covariate-inference.R:867 | read; covered by test-covariate-structure-inference.R | ok |
| `.multilpa_cov_blocks` | R/covariate-inference.R:897 | ran (hess.R) | ok |
| `.multilpa_cov_natural_estimate` | R/covariate-inference.R:920 | ran (hess.R, stored.R) | ok |
| `.multilpa_cov_natural_jacobian` | R/covariate-inference.R:955 | ran (hess.R) | ok |
| `.multilpa_cov_cholesky_natural` | R/covariate-inference.R:1008 | read; covered by test-covariate-measurement.R | ok |
| `.multilpa_cov_cholesky_jacobian` | R/covariate-inference.R:1028 | read; covered by test-covariate-measurement.R | ok |
| `.multilpa_cov_labels` | R/covariate-inference.R:1052 | ran (hess.R: labels align with SE rows) | ok |
| `coef.multilpa_covariates` | R/covariate-inference.R:1098 | ran (stored.R) | ok |
| `confint.multilpa_covariates` | R/covariate-inference.R:1119 | ran (stored.R); test-review-1101-confint-fixed.R | ok |
| `.latents_parse_direct_effects` | R/direct-effects.R:27 | read; covered by test-direct-effects.R | ok |
| `.latents_check_direct_conflicts` | R/direct-effects.R:66 | read; covered by test-direct-effects.R | ok |
| `.latents_add_direct` | R/direct-effects.R:83 | read; ran through direct.R | ok |
| `.latents_center_direct` | R/direct-effects.R:105 | ran (direct.R, covariate mean 40) | ok; logged (unweighted scales) |
| `.latents_uncenter_direct` | R/direct-effects.R:156 | ran (direct.R: means = lm intercepts) | ok |
| `.latents_direct_key` | R/direct-effects.R:169 | read; covered by test-features-direct-effects.R | ok |
| `.latents_direct_terms` | R/direct-effects.R:179 | read; covered by test-direct-effects.R | ok |
| `.latents_prepare_direct` | R/direct-effects.R:197 | read; covered by test-direct-effects.R | ok |
| `.latents_direct_covariates` | R/direct-effects.R:242 | read; covered by test-direct-effects.R | ok |
| `.latents_direct_matrix` | R/direct-effects.R:272 | ran (direct.R) | ok |
| `.latents_direct_shift` | R/direct-effects.R:285 | ran (direct.R, bvr residuals) | ok |
| `.latents_direct_maximize` | R/direct-effects.R:307 | ran (direct.R: = lm slopes) | ok |
| `.latents_direct_scores` | R/direct-effects.R:344 | ran (direct.R: SE = ML OLS SE) | ok |
| `.latents_direct_like` | R/direct-effects.R:367 | read; covered by test-features-direct-effects.R | ok |
| `.latents_direct_frame` | R/direct-effects.R:379 | ran (direct.R) | ok |
| `.cross_level_arguments` | R/cross-level.R:20 | read | ok |
| `.cross_level_group_density` | R/cross-level.R:28 | read; ran through cl.R | ok |
| `.cross_level_between_update` | R/cross-level.R:40 | read (divisors checked by hand); ran cl.R | ok |
| `.cross_level_initialize` | R/cross-level.R:63 | read; ran cl.R | ok |
| `.cross_level_em` | R/cross-level.R:128 | ran (cl.R, fam.R) | ok |
| `.cross_level_count` | R/cross-level.R:160 | read; hand count | ok |
| `.cross_level_fit` | R/cross-level.R:177 | ran (clarg.R) | fixed (F4) |
| `.cross_level_observed` | R/cross-level.R:289 | read | ok |
| `.cross_level_select` | R/cross-level.R:301 | ran (cl.R) | ok |
| `.cross_level_restricted` | R/cross-level.R:334 | ran (cl.R) | ok |
| `.cross_level_result` | R/cross-level.R:401 | ran (cl.R, fam.R) | ok |
| `.cross_level_tables` | R/cross-level.R:493 | read | ok |
| `get_results.multilpa_cross_level` | R/cross-level.R:500 | ran (clget.R) | fixed (F3) |
| `.cross_level_table` | R/cross-level.R:519 | ran (clget.R, fam.R) | ok |
| `as.data.frame.multilpa_cross_level` | R/cross-level.R:587 | ran (clget.R) | fixed (F3); row.names logged |
| `print.multilpa_cross_level` | R/cross-level.R:594 | ran (fam.R) | ok |
| `summary.multilpa_cross_level` | R/cross-level.R:612 | read; covered by test-cross-level.R | ok |
| `logLik.multilpa_cross_level` | R/cross-level.R:620 | read; covered by test-cross-level.R | ok |
| `nobs.multilpa_cross_level` | R/cross-level.R:627 | read; covered by test-cross-level.R | ok |
| `parameter_inference.multilpa_cross_level` | R/cross-level.R:633 | read | ok |
| `plot.multilpa_cross_level` | R/cross-level.R:642 | ran (clget.R) | fixed (F3) |
| `.latents_parse_pairs` | R/local-dependence.R:17 | read; covered by test-local-dependence.R | ok |
| `.latents_classify_pairs` | R/local-dependence.R:53 | read; covered by test-local-dependence.R | ok |
| `.latents_covariance_pattern` | R/local-dependence.R:77 | ran (ld.R) | ok |
| `.latents_pattern_pairs` | R/local-dependence.R:107 | ran (ld.R) | ok |
| `.latents_pattern_label` | R/local-dependence.R:121 | ran (ld.R) | ok |
| `.latents_mask_covariances` | R/local-dependence.R:131 | read; ran through ld.R | ok |
| `.latents_pattern_zeros` | R/local-dependence.R:143 | ran (ld.R, parameter count 15) | ok |
| `.latents_joined_names` | R/local-dependence.R:155 | ran (ld.R) | ok |
| `.latents_join_pairs` | R/local-dependence.R:183 | ran (ld.R, half-observed rows) | ok |
| `.latents_partial_label` | R/local-dependence.R:281 | ran (ld.R) | ok |
| `.latents_partial_levels` | R/local-dependence.R:297 | ran (ld.R) | ok |
| `.latents_partial_distinct` | R/local-dependence.R:305 | ran (ld.R) | ok |
| `.latents_partial_map` | R/local-dependence.R:330 | ran (ld.R) | ok |
| `.latents_restore_partial` | R/local-dependence.R:366 | ran (ld.R: labels restored exactly) | ok |
| `.latents_fit_local_dependence` | R/local-dependence.R:392 | ran (ld.R) | ok |
| `.latents_categorical_labels` | R/local-dependence.R:479 | ran (ld.R) | ok |
| `.latents_local_dependence_table` | R/local-dependence.R:516 | ran (ld.R: probabilities sum to 1; covariance/correlation by hand) | ok |
| `.multilpa_bivariate_residuals` | R/bivariate-residuals.R:13 | ran (bvr.R, ld.R, center.R) | ok; logged (BH family = all rows) |
| `.multilpa_cramers_v` | R/bivariate-residuals.R:87 | ran (bvr.R: V by hand) | ok |
| `.multilpa_gaussian_residuals` | R/bivariate-residuals.R:108 | ran (bvr.R: Fisher z by hand) | ok |
| `.multilpa_discrete_items` | R/bivariate-residuals.R:173 | ran (bvr.R, ld.R) | ok |
| `.multilpa_categorical_residuals` | R/bivariate-residuals.R:211 | ran (bvr.R: 2x2 chi-square by hand); test audit-1101 | ok; logged (V of expected ~1e-16) |
| `.multilpa_pair_values` | R/bivariate-residuals.R:275 | read; covered by test-bivariate-residuals.R | ok |
| `.multilpa_mixed_residuals` | R/bivariate-residuals.R:308 | read; covered by test-bivariate-residuals.R | ok |
| `.latents_gaussian_diagnostics` | R/gaussian-diagnostics.R:20 | ran (gd.R: every column by hand) | ok; logged (positive_rows, Kish) |
| `.latents_print_gaussian_diagnostics` | R/gaussian-diagnostics.R:90 | ran (stored.R print) | ok |
| `.latents_gaussian_minimum` | R/gaussian-diagnostics.R:112 | ran (gd.R) | ok |
| `.multilpa_entropy_sum` | R/diagnostics.R:7 | ran (diag.R, ic.R) | ok |
| `.multilpa_relative_entropy` | R/diagnostics.R:30 | ran (diag.R: weighted by hand) | ok |
| `.multilpa_criterion_definitions` | R/diagnostics.R:51 | ran (ic.R); formulas read against references | ok |
| `.multilpa_information_criteria` | R/diagnostics.R:88 | ran (ic.R: = fit aic/bic; ICL by hand) | ok |
| `.multilpa_criteria_wide` | R/diagnostics.R:172 | ran (ic.R) | ok |
| `.multilpa_criterion_labels` | R/diagnostics.R:196 | ran (ic.R) | ok |
| `.multilpa_posterior_levels` | R/diagnostics.R:208 | ran (diag.R) | ok |
| `.multilpa_classification_levels` | R/diagnostics.R:224 | ran (diag.R) | ok |
| `.multilpa_classification_table` | R/diagnostics.R:244 | ran (diag.R, ic.R) | ok |
| `.multilpa_level_weights` | R/diagnostics.R:268 | ran (diag.R: zero-weight rows excluded) | ok |
| `.multilpa_average_posteriors` | R/diagnostics.R:292 | ran (fam.R) | ok |
| `.multilpa_average_posterior_matrix` | R/diagnostics.R:315 | ran (fam.R) | ok |
| `.multilpa_classification_summary` | R/diagnostics.R:344 | ran (diag.R, ic.R) | ok |
| `.multilpa_average_posterior_rows` | R/diagnostics.R:372 | ran (fam.R) | ok |
| `.multilpa_class_labels` | R/diagnostics.R:394 | ran (fam.R noise fit: class 0) | ok |
| `.multilpa_odds_correct` | R/diagnostics.R:410 | read; ran (diag.R) | ok |
| `.multilpa_entropy_table` | R/diagnostics.R:433 | ran (diag.R, ic.R) | ok; logged (single-level groups row) |
| `.latents_class_label` | R/diagnostics-checks.R:7 | ran (fam.R) | ok |
| `.latents_check_row` | R/diagnostics-checks.R:27 | ran (fam.R) | ok |
| `.latents_estimation_checks` | R/diagnostics-checks.R:37 | ran (fam.R, 13 fits; cl.R) | ok |
| `.latents_classification_inputs` | R/diagnostics-checks.R:97 | ran (fam.R, 13 fits) | ok |
| `.latents_level_classification` | R/diagnostics-checks.R:139 | ran (fam.R) | ok |
| `.latents_level_posteriors` | R/diagnostics-checks.R:146 | ran (fam.R) | ok |
| `.latents_posterior_weights` | R/diagnostics-checks.R:160 | ran (fam.R) | ok |
| `.latents_posterior_classification` | R/diagnostics-checks.R:173 | ran (fam.R) | ok |
| `.latents_require_classification_weights` | R/diagnostics-checks.R:192 | read; covered by test-astra-diagnostics.R | ok |
| `.latents_classifying_levels` | R/diagnostics-checks.R:204 | ran (fam.R) | ok |
| `.latents_family_diagnostic_tables` | R/diagnostics-checks.R:237 | ran (fam.R) | ok |
| `.latents_level_checks` | R/diagnostics-checks.R:265 | ran (fam.R, diag.R) | ok; logged (noise class 0) |
| `.latents_family_smallest_share` | R/diagnostics-checks.R:338 | ran (fam.R regression/growth/additive) | ok |
| `.latents_level_entropy` | R/diagnostics-checks.R:351 | ran (fam.R) | ok |
| `.latents_residual_checks` | R/diagnostics-checks.R:359 | ran (bvr.R: count = BH by hand; diag.R weighted) | ok |
| `.latents_information_check` | R/diagnostics-checks.R:408 | ran (fam.R, 13 fits) | ok |
| `.latents_diagnostic_checks` | R/diagnostics-checks.R:426 | ran (fam.R) | ok |
| `.latents_print_checks` | R/diagnostics-checks.R:449 | ran (rep.R) | ok |
| `diagnostics` | R/report.R:139 | ran (fam.R, bvr.R, diag.R) | fixed (F5, docs) |
| `.multilpa_has_group_classes` | R/report.R:203 | ran (fam.R) | ok |
| `as.data.frame.multilpa_diagnostics` | R/report.R:210 | read; covered by test-diagnostics.R | ok |
| `summary.multilpa_diagnostics` | R/report.R:223 | read; covered by test-review-1027-diagnostics.R | ok |
| `print.multilpa_diagnostics` | R/report.R:236 | ran (rep.R) | ok |
| `.multilpa_replication_table` | R/report.R:320 | ran (fam.R, 13 fits; cl.R) | fixed (F2) |
| `.multilpa_replication_phrase` | R/report.R:369 | ran (stored.R/rep.R prints) | ok |
| `.multilpa_reject_extra_arguments` | R/report.R:407 | read; covered by test-review-1027-error-contract.R | ok |
| `.latents_warn_unused_arguments` | R/report.R:429 | ran (clget.R) | ok |
| `.latents_former_argument` | R/report.R:452 | read; covered by test-api-compat-0921.R | ok |
| `plot.multilpa_diagnostics` | R/report.R:465 | read; covered by test-diagnostics.R | ok |
| `report` | R/report.R:522 | ran (rep.R, 13 fits, plots on) | ok |
| `.multilpa_supported_views` | R/report.R:578 | ran (rep.R, other.R) | ok; needs-other-slice (probabilities view, 1 profile) |
| `.multilpa_available_views_text` | R/report.R:631 | read; covered by test-report-surface.R | ok |
| `.latents_residual_alpha` (constant) | R/diagnostics-checks.R:16 | read; documented on ?diagnostics | ok |
| `.latents_check_references` (constant) | R/diagnostics-checks.R:18 | read | ok |

## Fixed

**F1. Covariate fits warned of thin profiles that are not thin, with a wrong
parameter count.** `lpa(model = "VEI", profile_covariates = "z")` on 150 +
7 units warned "Profile 2: 7.0 Gaussian membership mass for 8 free means and
variances" while the same model without covariates did not (VEI has 4 means
+ 1 volume = 5 own parameters); likewise EVI and VEE, and a covariate fit
with `fixed = "variances"` (variances held, so nothing can collapse).
Reproduction (`thin.R`, before): `VEI plain: | covariate:
latents_thin_profiles`, `EVI ... covariate: latents_thin_profiles`, `VEE
... covariate: latents_thin_profiles`, `fixed variances plain: | covariate:
latents_thin_profiles`. Cause: `.multilpa_fit_covariates()` called
`.latents_profile_parameters()` without `structure` and `fixed`, so it
counted the VVI/VVV spread. Change: R/covariates.R:555-561 passes
`structure, fixed = fixed`. After: no warning for VEI/EVI/VEE/held, VVI and
VVV still warn. Test: "a covariate fit counts a constrained structure's own
parameters as multilpa() does", "a covariate fit holding its variances has
no thin profile to report".

**F2. Restricted cross-level diagnostics reported the factors as starts.**
A fit with `n_starts = 5` printed "best log likelihood reached by 2 of 2
starts" and the replication table `n_starts = 2, n_converged = 2, share =
1`: its `starts` table has one row per factor (`part`), not per start.
Cause: `.multilpa_replication_table()` counted rows of any table with
`log_likelihood` and `converged`. Change: R/report.R:320-331 treats a
per-factor record as no per-start record: `n_starts` comes from the fit
(5), share = 2/5 = 0.4, converged/failed counts are NA (not recorded), so
the "starts failed or unfinished" row is omitted rather than invented.
Test: "the restricted cross-level replication counts starts, not its two
factors".

**F3. Cross-level accessors raised an unclassed error or ignored
arguments.** `get_results(fit, "bogus")`, `as.data.frame(fit, what =
"bogus")` and `plot(fit, what = "bogus")` raised `simpleError` from
`match.arg()`; `get_results(fit, "groups", level = "groups")` returned the
table without a word about `level`. Change: R/cross-level.R:500-509 uses
`.latents_match_arg()` (same partial matching, class
`latents_bad_argument`) and warns `latents_unused_argument` for extras (the
0.9.21 policy for methods that ignored arguments); `as.data.frame()` passes
`...` on (587-590); `plot()` uses `.latents_match_formal()` (649). Test:
"cross-level tables refuse a misspelt table and name an ignored argument".

**F4. An invalid cross-level variance mode ended as "all starts failed".**
`multilpa(..., family = "full_cross_level", between_variance = "bogus")`
raised `latents_all_starts_failed`; `variance_model = "bogus"` likewise.
Cause: multilpa() passes `.latents_match_formal(...)` as an unevaluated
promise, first forced inside a start's error handler. Change:
R/cross-level.R:181-192 forces both and validates them with
`latents_bad_argument` before any start. Test: "an invalid cross-level
variance mode is refused, not counted as failed starts".

**F5 (documentation). `?diagnostics` named the wrong condition.** It said a
family without bivariate residuals raises `latents_no_group_classes` when the
table is asked for by name; the table is absent from that object's
catalogue and `get_results()` raises `latents_bad_argument` (checked on
additive, cross-level, regression and growth fits). Change:
R/report.R:93-98. Test: the last expectation of the F3 test.

Every fix was broken once on purpose (`sabotage.py`): F1 failed 2 tests, F2
failed 5 expectations, F3 and F4 errored their tests; all restored and
`diff -q` against the backups clean.

## Logged, not fixed

- `.multilpa_categorical_residuals`: the Cramér's V of the model-implied
  table is ~1e-16 instead of exactly 0 (rounding in the product of margins);
  no test or verdict reads it.
- `.multilpa_entropy_table`: a single-level fit carries a `groups` row (one
  class, relative entropy NA) in `get_results(x, "entropy")` and
  `diagnostics()`; frozen table layout, the checks already skip it.
- `.latents_center_direct`: `scales` is the unweighted SD of the effect's
  covariate; it only scales finite-difference steps, no reported number.
- `.multilpa_bivariate_residuals`: the BH family is every row (all pairs in
  all profiles under `by = "profile"`); a documented design choice.
- `diagnostics()`: the `latents_no_group_classes` handler around the
  residuals is unreachable (the residual builder never raises it); harmless.
- Noise fits: "smallest class" and "classes with no assigned unit" count the
  noise component (class 0); an empty noise component gives status "check".
- `as.data.frame.multilpa_cross_level` ignores `row.names` (now warns on any
  other extra argument).
- `.latents_gaussian_diagnostics`: `positive_rows` is the same for every
  component by design (comment GD-1).
- Kish `weight_effective_n` exceeds `gaussian_mass` for unweighted fits
  (posteriors ≤ 1); this is the definition, not an error.

## Needs another slice

- R/plot-gg.R:859 `.gg_view_probabilities()`: on a one-profile fit the
  "probabilities" view draws a one-point line per group class, and ggplot2
  prints "`geom_line()`: Each group consists of only one observation" during
  `report(lpa(..., n_profiles = 1))`. Either the view refuses a single
  profile (`latents_nothing_to_plot`, which `report()` already lists as not
  drawn) or `.multilpa_supported_views()` (R/report.R:578) drops it; not
  changed here because it changes plotted output.
- R/multilpa.R:1624-1627 hands `.latents_match_formal(variance_model)` and
  `.latents_match_formal(between_variance)` to `.cross_level_fit()` as lazy
  promises; F4 forces them at the receiving end. The additive path refuses
  `between_variance = "bogus"` correctly (`latents_bad_argument`, checked);
  any other dispatch that passes a match as an argument expression and first
  forces it inside a per-start `tryCatch` would turn a refusal into failed
  starts the same way.

## Test runs

One R process, `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`,
`testthat::test_file()` per file after all fixes:

```
test-audit-1101-covariates-diagnostics.R: PASS 38 FAIL 0 ERROR 0 SKIP 0
test-covariate-fixed.R: PASS 16 FAIL 0 ERROR 0 SKIP 0
test-covariate-fold.R: PASS 30 FAIL 0 ERROR 0 SKIP 0
test-covariate-inference.R: PASS 40 FAIL 0 ERROR 0 SKIP 0
test-covariate-measurement.R: PASS 33 FAIL 0 ERROR 0 SKIP 0
test-covariate-missing.R: PASS 16 FAIL 0 ERROR 0 SKIP 0
test-covariate-separation.R: PASS 12 FAIL 0 ERROR 0 SKIP 0
test-covariate-slopes.R: PASS 14 FAIL 0 ERROR 0 SKIP 0
test-covariate-structure-inference.R: PASS 65 FAIL 0 ERROR 0 SKIP 0
test-covariate-surface.R: PASS 40 FAIL 0 ERROR 0 SKIP 0
test-covariates.R: PASS 23 FAIL 0 ERROR 0 SKIP 0
test-direct-effects.R: PASS 55 FAIL 0 ERROR 0 SKIP 0
test-features-direct-effects.R: PASS 78 FAIL 0 ERROR 0 SKIP 0
test-cross-level.R: PASS 14 FAIL 0 ERROR 0 SKIP 0
test-cross-level-thin-profiles.R: PASS 48 FAIL 0 ERROR 0 SKIP 0
test-local-dependence.R: PASS 50 FAIL 0 ERROR 0 SKIP 0
test-bivariate-residuals.R: PASS 41 FAIL 0 ERROR 0 SKIP 0
test-diagnostics.R: PASS 133 FAIL 0 ERROR 0 SKIP 0
test-review-1027-diagnostics.R: PASS 158 FAIL 0 ERROR 0 SKIP 0
test-report-surface.R: PASS 40 FAIL 0 ERROR 0 SKIP 0
test-review-1101-confint-fixed.R: PASS 77 FAIL 0 ERROR 0 SKIP 0
test-api-compat-0921.R: PASS 44 FAIL 0 ERROR 0 SKIP 0
test-compat-0921.R: PASS 10 FAIL 0 ERROR 0 SKIP 0
test-astra-diagnostics.R: PASS 80 FAIL 0 ERROR 0 SKIP 0
test-tidy-diagnostics.R: PASS 51 FAIL 0 ERROR 0 SKIP 0
test-gaussian-diagnostics.R: PASS 41 FAIL 0 ERROR 0 SKIP 0
test-gaussian-diagnostic-tables.R: PASS 77 FAIL 0 ERROR 0 SKIP 0
test-audit-reporting-repairs.R: PASS 25 FAIL 0 ERROR 0 SKIP 0
test-review-1027-replication.R: PASS 29 FAIL 0 ERROR 0 SKIP 0
test-thin-profiles.R: PASS 39 FAIL 0 ERROR 0 SKIP 0
test-review-1015-covariates.R: PASS 32 FAIL 0 ERROR 0 SKIP 0
test-gaussian-warning-support.R: PASS 28 FAIL 0 ERROR 0 SKIP 0
test-review-1027-error-contract.R: PASS 216 FAIL 0 ERROR 0 SKIP 0
test-features-pair-covariance.R: PASS 69 FAIL 0 ERROR 0 SKIP 0
```

MoEClust kit (`testthat::test_dir("equivalence", filter =
"^moeclust-covariate-structures$")`, as run-package-kit-tests.R runs it,
without recording evidence), after F1:

```
PASS 55 FAIL 0 ERROR 0 SKIP 0
```

`devtools::document()`: NAMESPACE unchanged; man/diagnostics.Rd changed (F5);
it also rewrote man/compare_models.Rd from another slice's roxygen.

`Rscript equivalence/engine-golden/compare.R` (8.8 min): 242 identical, 18
DIFFERENT: 13 the bootstrap attribute rename (`replicate_rows` ->
`replicate_identity`); 5 others, none from this slice's files:
`cov_weights`, `family_additive_weights`,
`reg_gaussian_weights_observed_refused` (refusal message text from
R/weights.R), `lta_refusal_no_time` and `svc_fit_staged_one_class`
(`simpleError` -> `latents_bad_argument`, other slices). The 22
reg_*/growth_lcga `group_classes@display` differences did not appear in
this run.
