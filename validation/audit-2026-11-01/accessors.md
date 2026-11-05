# Audit 2026-11-01 — slice `accessors`

Files: R/accessors.R, R/methods.R, R/methods-documentation.R, R/get-data.R,
R/latents-table.R, R/predict.R, R/predict-docs.R, R/descriptives.R, R/data.R,
R/pool-imputations.R, R/simulation-service.R, R/conditions.R, R/condition-calls.R.
150 functions plus 9 documentation objects (6 datasets, 3 topics) checked;
4 code defects fixed, 4 documentation defects fixed, 8 items logged,
3 items need another slice. No estimate, default, signature, table layout or
printed fit table changed for valid input.

Scratch evidence scripts (not part of the package):
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/accessors/`
(`bench.R` table/column dump, `checks.R` table consistency, `pool.R` Rubin by
hand, `pred.R` predict, `rec.R` recovery/ICC, `ord.R` zero weights,
`data.R` dataset claims, `cond.R`/`inherit.R` condition registry, `misc.R`
printing and services).

## Independent checks (no defect found)

- **Table consistency** on five fits (two-level LPA, single-level LPA, mixed
  LCA with a continuous item, membership-covariate fit, LTA): `logLik()` ==
  fit, `df` == `n_parameters`, `AIC()`/`BIC()` == stored `aic`/`bic`,
  `nobs()` == groups; wide posteriors sum to 1, long == wide, assignments'
  `profile` == wide modal, `uncertainty` == 1 − max; individual counts ==
  `colSums(posteriors)`; profile means == `fit$means`, sd == sqrt(var),
  covariance diagonal == variances; standardized means and SEs recomputed by
  hand from `get_results(x, "data")`; `summary(x)$tables` identical to
  `get_results(x, "all")`; information criteria aic/bic/sabic/caic/kic
  recomputed from L and q; replication share. All 24 checks pass on all 5.
- **Columns vs roxygen**: every table of `get_results(x, "all")` dumped for the
  five fits, a starting-value set, a diagnostics object, an enumeration grid
  and a comparison. Mismatches found: `"replication"` (3 undocumented columns)
  and the diagnostics `"checks"` table (unlisted) — fixed. `"criteria"`
  duplicated beside `"best"` in a grid's `"all"` is frozen by
  test-api-compat-0921.R (deliberate).
- **Rubin's rules by hand** (3 imputations, covariate fit): Q, U, B, T,
  df = (m − 1)(1 + 1/r)², fmi = (r + 2/(df + 3))/(r + 1), t intervals,
  p-values and the full pooled covariance U + (1 + 1/m)·cov(Q) all equal the
  package's values. The df is Rubin's (1987) large-sample df, as documented;
  the package does not use Barnard–Rubin, and says so.
- **Follow-up 2** (pooling covariate fits with `fixed = "measurement"`): works;
  `vcov()` 2×2 over the two estimated logits, sqrt(diag) == pooled SEs,
  `coef()` names == `vcov()` names. Regression test added.
- **predict()**: posteriors on the fitted rows == stored posteriors, newdata
  = the fitting data gives identical output, group class == group posteriors,
  sum of single-level log densities == logLik (3 fits).
- **Recovery** (`truth =`): reproduced with `table()`; group-level truth
  counts each group once; proportions sum to 1 per truth value.
- **ICC(1)**: 0.1789578 from `descriptives()` and from `anova(lm())` mean
  squares with n0.
- **Datasets**: dims, types, ranges, factor levels and every numerical claim in
  R/data.R (per-student counts, 49.7% `worried == 1`, srl correlations,
  12.2% / 7.3% omitted occasions, 15 students per school, constant
  motivation/programme) match the shipped data.
- **Condition registry**: 98 `latents_` classes raised (errorCondition,
  warningCondition, classed `message(structure())`, `.multilpa_notice()`) ==
  the classes documented in R/conditions.R, apart from `latents_partial_pair`,
  which is documented as "Not raised". Every "also inherits from" claim
  matches the class vectors at the raising sites.
- **Services**: `.latents_draw_rows()` reproduces (0.2, 0.3, 0.5) over 1e5
  draws; `.latents_blrt_p_value()` = (1 + #T* ≥ T)/(V + 1) by hand (0.75) and
  warns `latents_failed_replicates`; replicate errors and warnings recorded;
  `.latents_bounded_call()` turns a 100000-row data frame into
  `<data.frame[100000 x 2]>`; `.latents_bind_values()` avoids a caller's `.x`.

## Per-function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `as.data.frame.multilpa` | R/accessors.R:3 | ran synthetic: equals get_results(x) on bench objects; test-get-data.R | ok |
| `.multilpa_assignments` | R/accessors.R:29 | ran synthetic: profile == wide modal, uncertainty = 1 - max posterior (all 5 bench fits); zero-weight fit; test-audit-1101 (zero weight) | ok |
| `.multilpa_recovery_frame` | R/accessors.R:78 | ran synthetic: course_engagement recovery cross-checked against table() of assignments (both levels, groups counted once); proportions sum to 1 per truth value | ok |
| `.multilpa_recovery_rows` | R/accessors.R:109 | ran synthetic: course_engagement recovery cross-checked against table() of assignments (both levels, groups counted once); proportions sum to 1 per truth value | ok |
| `.multilpa_constant_within` | R/accessors.R:161 | ran synthetic: course_engagement recovery cross-checked against table() of assignments (both levels, groups counted once); proportions sum to 1 per truth value | ok |
| `.multilpa_model_frame` | R/accessors.R:186 | ran synthetic: data rows = n_obs (5 fits); weighted fit with zero-weight ordinal/categorical rows | needs-other-slice |
| `.multilpa_resolve_data` | R/accessors.R:254 | ran synthetic: assignments(data = d) with weighted/categorical fit; covered by test-get-data.R | ok |
| `.multilpa_same_values` | R/accessors.R:275 | read + ran synthetic (excluded categorical value, zero-weight ordinal row accepted); covered by test-review-1015-accessors.R, test-get-data.R | ok |
| `.multilpa_alignment_values` | R/accessors.R:297 | read + ran synthetic (excluded categorical value, zero-weight ordinal row accepted); covered by test-review-1015-accessors.R, test-get-data.R | ok |
| `.multilpa_reconstructed_csv_values` | R/accessors.R:332 | read + ran synthetic (excluded categorical value, zero-weight ordinal row accepted); covered by test-review-1015-accessors.R, test-get-data.R | ok |
| `.multilpa_indicator_reproduced` | R/accessors.R:369 | read + ran synthetic (excluded categorical value, zero-weight ordinal row accepted); covered by test-review-1015-accessors.R, test-get-data.R | ok |
| `.multilpa_check_row_count` | R/accessors.R:392 | ran synthetic: assignments(data = d) with weighted/categorical fit; covered by test-get-data.R | ok |
| `.multilpa_check_alignment` | R/accessors.R:418 | read + ran synthetic (excluded categorical value, zero-weight ordinal row accepted); covered by test-review-1015-accessors.R, test-get-data.R | ok |
| `.multilpa_stage_frame` | R/accessors.R:513 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA); covariance diagonal == variances | ok |
| `.multilpa_continuous_names` | R/accessors.R:539 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA); covariance diagonal == variances | ok |
| `.multilpa_profile_frame` | R/accessors.R:547 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA): means == fit$means, sd = sqrt(var), SEs present | ok |
| `.multilpa_standardization` | R/accessors.R:613 | ran synthetic: standardized mean and SE recomputed by hand from get_results(x, "data") (5 fits) | logged |
| `.multilpa_measurement_errors` | R/accessors.R:644 | read + ran synthetic (SE columns filled on bench fits); refusal path covered by test-review-1015-accessors.R | ok |
| `.multilpa_table_errors` | R/accessors.R:662 | read + ran synthetic (SE columns filled on bench fits); refusal path covered by test-review-1015-accessors.R | ok |
| `.multilpa_match_error` | R/accessors.R:692 | read + ran synthetic (SE columns filled on bench fits); refusal path covered by test-review-1015-accessors.R | ok |
| `.multilpa_variance_error` | R/accessors.R:714 | read + ran synthetic (SE columns filled on bench fits); refusal path covered by test-review-1015-accessors.R | ok |
| `.multilpa_response_frame` | R/accessors.R:729 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA) (LCA: 16 rows, thresholds, SEs) | ok |
| `.multilpa_probability_frame` | R/accessors.R:773 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA); covariance diagonal == variances | ok |
| `.multilpa_posterior_frame` | R/accessors.R:794 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA): long == wide, rows sum to 1, group class == predict() | ok |
| `.multilpa_group_posterior_frame` | R/accessors.R:827 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA): long == wide, rows sum to 1, group class == predict() | ok |
| `as.data.frame.multilpa_enumeration` | R/accessors.R:858 | ran synthetic: equals get_results(x) on bench objects; test-get-data.R | ok |
| `.multilpa_is_probability` | R/accessors.R:868 | read | ok |
| `print.multilpa_enumeration` | R/accessors.R:875 | ran synthetic (printed on bench objects; printed means == profiles table) | ok |
| `starting_values` | R/accessors.R:987 | ran synthetic: get_results(starting_values(fit), "all"); measurement-only start used for held-measurement refits | ok |
| `print.multilpa_start` | R/accessors.R:1069 | ran synthetic (printed on bench objects; printed means == profiles table) | ok |
| `summary.multilpa_start` | R/accessors.R:1090 | ran synthetic (printed on bench objects; printed means == profiles table) | ok |
| `as.data.frame.multilpa_start` | R/accessors.R:1102 | ran synthetic: equals get_results(x) on bench objects; test-get-data.R | ok |
| `.multilpa_start_profile_frame` | R/accessors.R:1111 | ran synthetic: get_results(starting_values(fit), "all"); measurement-only start used for held-measurement refits | ok |
| `.multilpa_start_variances` | R/accessors.R:1134 | ran synthetic: get_results(starting_values(fit), "all"); measurement-only start used for held-measurement refits | ok |
| `.multilpa_start_response_frame` | R/accessors.R:1148 | ran synthetic: get_results(starting_values(fit), "all"); measurement-only start used for held-measurement refits | ok |
| `.multilpa_start_probability_frame` | R/accessors.R:1179 | ran synthetic: get_results(starting_values(fit), "all"); measurement-only start used for held-measurement refits | ok |
| `.multilpa_start_covariance_frame` | R/accessors.R:1201 | ran synthetic: get_results(starting_values(fit), "all"); measurement-only start used for held-measurement refits | ok |
| `as.data.frame.multilpa_covariates` | R/accessors.R:1224 | ran synthetic: equals get_results(x) on bench objects; test-get-data.R | ok |
| `.multilpa_coefficient_frame` | R/accessors.R:1233 | ran synthetic (covariate fit: 2 rows, matches coef()) | ok |
| `print.multilpa` | R/methods.R:3 | ran synthetic (printed on bench objects; printed means == profiles table) | ok |
| `summary.multilpa` | R/methods.R:94 | ran synthetic: summary tables identical to get_results(x, "all") (5 fits) | ok |
| `.multilpa_summary_advice` | R/methods.R:171 | read | ok |
| `as.data.frame.summary_multilpa` | R/methods.R:179 | ran synthetic: equals get_results(x) on bench objects; test-get-data.R | ok |
| `.multilpa_count_frame` | R/methods.R:194 | ran synthetic: individual counts == colSums(posteriors); zero-weight rows excluded (118 of 120) | logged |
| `.multilpa_covariance_frame` | R/methods.R:229 | ran synthetic (bench: 2-level LPA, single-level LPA, mixed LCA, covariate fit, LTA); covariance diagonal == variances | ok |
| `.multilpa_fit_frame` | R/methods.R:262 | ran synthetic: ll/aic/bic_groups/bic_individual == stored fit (5 fits) | ok |
| `.multilpa_one` | R/methods.R:309 | read | ok |
| `.multilpa_covariate_fit_frame` | R/methods.R:317 | ran synthetic: ll/aic/bic_groups/bic_individual == stored fit (5 fits) | ok |
| `print.summary_multilpa` | R/methods.R:344 | ran synthetic (printed on bench objects; printed means == profiles table) | ok |
| `logLik.multilpa` | R/methods.R:391 | ran synthetic: logLik == fit, df == n_parameters, AIC/BIC == stored aic/bic, nobs == groups (zero weights excluded); test-audit-1101 | ok |
| `nobs.multilpa` | R/methods.R:400 | ran synthetic: logLik == fit, df == n_parameters, AIC/BIC == stored aic/bic, nobs == groups (zero weights excluded); test-audit-1101 | ok |
| `get_results` | R/get-data.R:738 | ran synthetic: every catalogue table compared with the roxygen catalogue | fixed (doc) |
| `get_results.default` | R/get-data.R:744 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.multilpa` | R/get-data.R:753 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.multilpa_covariates` | R/get-data.R:759 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.multilpa_transitions` | R/get-data.R:765 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.multilpa_enumeration` | R/get-data.R:772 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.latents_comparison` | R/get-data.R:778 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.multilpa_bootstrap_lrt` | R/get-data.R:784 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.multilpa_start` | R/get-data.R:790 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.multilpa_diagnostics` | R/get-data.R:796 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.summary_multilpa` | R/get-data.R:802 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.summary_multilpa_covariates` | R/get-data.R:808 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.summary_multilpa_transitions` | R/get-data.R:815 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.summary_multilpa_enumeration` | R/get-data.R:821 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `get_results.summary_multilpa_bootstrap_lrt` | R/get-data.R:827 | ran synthetic: every table of "all" dumped (columns vs roxygen) on fits, start, diagnostics, grid, comparison; test-get-data.R | ok |
| `.multilpa_coerce` | R/get-data.R:841 | ran synthetic: equals get_results(x) on bench objects; test-get-data.R | ok |
| `.multilpa_dispatch_table` | R/get-data.R:861 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_match_table` | R/get-data.R:882 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_check_table_arguments` | R/get-data.R:908 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_absent_table_conditions` | R/get-data.R:932 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_all_tables` | R/get-data.R:957 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_catalogue` | R/get-data.R:989 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_table` | R/get-data.R:1074 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_any_summary` | R/get-data.R:1082 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_stored_catalogue` | R/get-data.R:1099 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_measurement_catalogue` | R/get-data.R:1116 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_mixing_catalogue` | R/get-data.R:1139 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_unit_catalogue` | R/get-data.R:1162 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_diagnostic_catalogue` | R/get-data.R:1187 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_sequence_catalogue` | R/get-data.R:1222 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_provenance_catalogue` | R/get-data.R:1239 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_hide_fabricated_id` | R/get-data.R:1260 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_default_level` | R/get-data.R:1270 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_diagnostic_residuals` | R/get-data.R:1281 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_sequence_length_frame` | R/get-data.R:1298 | ran synthetic via get_results(x, "all") on 9 objects; test-get-data.R, test-api-compat-0921.R | ok |
| `.multilpa_print_tables` | R/get-data.R:1321 | ran synthetic (print(fit), print(summary(fit))); wide means == long table | ok |
| `.multilpa_print_table_footer` | R/get-data.R:1351 | ran synthetic (print(fit), print(summary(fit))); wide means == long table | ok |
| `.multilpa_check_print_arguments` | R/get-data.R:1364 | ran synthetic (print(fit), print(summary(fit))); wide means == long table | ok |
| `.multilpa_wide_means` | R/get-data.R:1387 | ran synthetic (print(fit), print(summary(fit))); wide means == long table | ok |
| `.multilpa_wide_responses` | R/get-data.R:1401 | ran synthetic (print(fit), print(summary(fit))); wide means == long table | ok |
| `.multilpa_print_class_sizes` | R/get-data.R:1420 | ran synthetic (print(fit), print(summary(fit))); wide means == long table | ok |
| `.multilpa_print_primary` | R/get-data.R:1453 | ran synthetic (print(fit), print(summary(fit))); wide means == long table | ok |
| `.latents_table` | R/latents-table.R:18 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_column` | R/latents-table.R:37 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_column_decimals` | R/latents-table.R:44 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_label` | R/latents-table.R:54 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_format_p` | R/latents-table.R:76 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_cells` | R/latents-table.R:82 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_default_display` | R/latents-table.R:116 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_render_card` | R/latents-table.R:131 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `.latents_render` | R/latents-table.R:149 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `print.latents_table` | R/latents-table.R:183 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `as.data.frame.latents_table` | R/latents-table.R:194 | ran synthetic latents_table (label, number, ci, p formats; as.data.frame strips attributes) | ok |
| `predict.multilpa` | R/predict.R:6 | ran synthetic: posterior == stored posteriors, newdata = fitting data identical, sum(log density) == logLik (single level); misnamed argument and bad type reproduced | fixed |
| `predict.multilpa_covariates` | R/predict.R:55 | read (classed refusal; covered by test-predict.R) | ok |
| `predict.multilpa_transitions` | R/predict.R:66 | read (classed refusal; covered by test-predict.R) | ok |
| `.multilpa_training_rows` | R/predict.R:76 | ran synthetic through predict() on 3 fits + zero-weight fit; test-predict.R, test-predict-empty-rows.R | ok |
| `.multilpa_fit_data` | R/predict.R:86 | ran synthetic through predict() on 3 fits + zero-weight fit; test-predict.R, test-predict-empty-rows.R | ok |
| `.multilpa_prediction_rows` | R/predict.R:104 | ran synthetic through predict() on 3 fits + zero-weight fit; test-predict.R, test-predict-empty-rows.R | ok |
| `.multilpa_row_log_density` | R/predict.R:214 | ran synthetic through predict() on 3 fits + zero-weight fit; test-predict.R, test-predict-empty-rows.R | ok |
| `.multilpa_prediction_table` | R/predict.R:228 | ran synthetic through predict() on 3 fits + zero-weight fit; test-predict.R, test-predict-empty-rows.R | ok |
| `descriptives` | R/descriptives.R:95 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `descriptives.default` | R/descriptives.R:101 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `descriptives.data.frame` | R/descriptives.R:111 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `descriptives.multilpa` | R/descriptives.R:145 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `descriptives.multilpa_transitions` | R/descriptives.R:179 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `descriptives.multilpa_covariates` | R/descriptives.R:183 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `.multilpa_assignment_values` | R/descriptives.R:194 | ran synthetic (zero-weight fit: weight-0 rows go to the NA stratum while assignments give them a profile) | logged |
| `.multilpa_describe_by` | R/descriptives.R:229 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `.multilpa_describe` | R/descriptives.R:284 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `.multilpa_check_descriptive_columns` | R/descriptives.R:316 | ran synthetic: by profile / group_class on 2-level fit, zero-weight fit; test-descriptives-safety.R | ok |
| `.multilpa_icc` | R/descriptives.R:347 | ran synthetic: ICC(1) recomputed from anova(lm()) mean squares and n0 (0.1789578 both) | ok |
| `pool_imputations` | R/pool-imputations.R:133 | ran synthetic: 3 imputations, covariate fit, with and without fixed = "measurement" (follow-up 2); test-audit-1101 | ok |
| `.multilpa_completed_data` | R/pool-imputations.R:189 | ran synthetic through pool_imputations(); test-pool-imputations.R | ok |
| `.multilpa_pool_step` | R/pool-imputations.R:234 | ran synthetic through pool_imputations(); test-pool-imputations.R | ok |
| `.multilpa_pool_keys` | R/pool-imputations.R:247 | ran synthetic through pool_imputations(); test-pool-imputations.R | ok |
| `.multilpa_rubin` | R/pool-imputations.R:264 | ran synthetic: Q, U, B, T, df, fmi, CI, p and the full pooled covariance recomputed by hand (test-audit-1101) | ok |
| `.multilpa_rubin_covariance` | R/pool-imputations.R:307 | ran synthetic: Q, U, B, T, df, fmi, CI, p and the full pooled covariance recomputed by hand (test-audit-1101) | ok |
| `.multilpa_pool_long` | R/pool-imputations.R:330 | ran synthetic through pool_imputations(); test-pool-imputations.R | ok |
| `.multilpa_pool_fit_frame` | R/pool-imputations.R:344 | ran synthetic through pool_imputations(); test-pool-imputations.R | ok |
| `get_results.latents_pooled` | R/pool-imputations.R:362 | ran synthetic: "all" refused, extra argument dropped silently | fixed |
| `as.data.frame.latents_pooled` | R/pool-imputations.R:377 | ran synthetic | fixed (guards what = "all") |
| `print.latents_pooled` | R/pool-imputations.R:386 | ran synthetic | logged |
| `summary.latents_pooled` | R/pool-imputations.R:405 | ran synthetic | logged |
| `vcov.latents_pooled` | R/pool-imputations.R:426 | ran synthetic: vcov(pooled, type = "robust") silently ignored | fixed |
| `coef.latents_pooled` | R/pool-imputations.R:444 | ran synthetic: confint at level 0.95 == conf_low/conf_high; names == vcov names | ok |
| `confint.latents_pooled` | R/pool-imputations.R:466 | ran synthetic: confint at level 0.95 == conf_low/conf_high; names == vcov names | ok |
| `.multilpa_pooled_names` | R/pool-imputations.R:505 | ran synthetic: confint at level 0.95 == conf_low/conf_high; names == vcov names | ok |
| `plot.latents_pooled` | R/pool-imputations.R:514 | read (Okabe-Ito, colour paired with shape); not changed, no render | ok |
| `.latents_local_seed` | R/simulation-service.R:46 | read; covered by seed and bootstrap tests (test-tidy-core.R) | ok |
| `.latents_draw_rows` | R/simulation-service.R:80 | ran synthetic: 1e5 draws reproduce (0.2, 0.3, 0.5); one-column and one-row inputs | logged |
| `.multilpa_resample_groups` | R/simulation-service.R:100 | read; covered by seed and bootstrap tests (test-tidy-core.R) | ok |
| `.latents_cluster_bootstrap` | R/simulation-service.R:139 | ran synthetic (stub refit/replicate: failure rows, warnings recorded, p = (1 + #T* >= T)/(V + 1) by hand) | ok |
| `.latents_blrt_replicates` | R/simulation-service.R:177 | ran synthetic (stub refit/replicate: failure rows, warnings recorded, p = (1 + #T* >= T)/(V + 1) by hand) | ok |
| `.latents_blrt_p_value` | R/simulation-service.R:227 | ran synthetic (stub refit/replicate: failure rows, warnings recorded, p = (1 + #T* >= T)/(V + 1) by hand) | ok |
| `.latents_bind_values` | R/condition-calls.R:26 | ran synthetic (100000-row data frame bounded to a placeholder; symbol collision with caller `.x` resolved) | ok |
| `.latents_bounded_call` | R/condition-calls.R:70 | ran synthetic (100000-row data frame bounded to a placeholder; symbol collision with caller `.x` resolved) | ok |
| `.latents_value_placeholder` | R/condition-calls.R:84 | ran synthetic (100000-row data frame bounded to a placeholder; symbol collision with caller `.x` resolved) | ok |
| `.latents_bounded_condition` | R/condition-calls.R:94 | ran synthetic (100000-row data frame bounded to a placeholder; symbol collision with caller `.x` resolved) | ok |
| `course_engagement`, `student_esm`, `study_hours`, `srl`, `growth_scores`, `growth_schools` (dataset docs) | R/data.R:100-392 | ran: dims, column types, ranges, levels, per-student counts, ICC-free claims (worried = 1 share 0.497, srl correlations 0.78-0.86 / -0.30 to -0.10, 12.2% and 7.3% omitted occasions, 15 students per school, constant motivation/programme) all match | ok |
| `latents-conditions` (topic) | R/conditions.R:456 | grep + regex audit: 98 raised `latents_` classes == documented set (less `latents_partial_pair`, documented "Not raised"); inheritance claims checked against class vectors; new working registry test | ok (logged: partial_pair, multilpa_indicator_reused) |
| `latents-model-methods` (topic) | R/methods-documentation.R:92 | ran coef/vcov/confint on held-measurement fits (follow-up 1) | fixed (doc) |
| `latents-predict` (topic) | R/predict-docs.R:131 | ran predict on zero-weight fit | fixed (doc) |

## Fixed

1. **`predict()` of a profile/class model dropped a misnamed argument and gave
   an unclassed error for a bad `type`.** `predict(fit, data = srl)` returned
   the fitted rows' predictions with no warning (the caller meant
   `newdata =`); `predict(fit, type = "bogus")` raised `simpleError` from
   `match.arg()`. Cause: `...` unread, `match.arg()`. Change: R/predict.R:8-13
   uses `.latents_match_arg()` (still partial-matches, so `type = "post"`
   works) and `.latents_warn_unused_arguments()` (`latents_unused_argument`,
   as the growth-mixture `predict()` already does). Doc: R/predict-docs.R
   Conditions section. Test: "predict() of a profile model refuses a bad type
   and warns about a misnamed argument" (sabotage: reverting both lines made
   it fail and error).
2. **`get_results(pooled, "all")` was refused** (`latents_bad_argument`)
   although ?get_results promises `"all"` for every object, including pooled
   imputations ("Tables of the other model families"); and
   `get_results(pooled, "estimates", level = 0.9)` silently ignored `level`,
   where every other `get_results()` method refuses an unused argument.
   Change: R/pool-imputations.R:362-373 accepts `"all"` (the three tables,
   named) and refuses extra arguments through
   `.multilpa_check_table_arguments()`; `as.data.frame.latents_pooled()`
   (R/pool-imputations.R:377-381) keeps returning one table and refuses
   `what = "all"`. Test: "get_results() of pooled imputations offers \"all\"
   and refuses unused arguments" (sabotage: reverting made it error).
3. **`vcov(pooled, type = "robust")` silently returned the pooled covariance**
   while `coef()`/`confint()` of the same object warn on unused arguments.
   Change: R/pool-imputations.R:426-431 warns `latents_unused_argument`.
   Test: same block (sabotage: removing the warning failed the
   `expect_warning`).
4. **Stale zero-weight documentation.** ?get_results ("assignments") and
   ?predict said units of sampling weight zero are unclassified; since 1.1.3
   (ML-1, owner rule) they keep their row, posteriors and assignment
   (reproduced: rows 5 and 70 of a weight-0 fit have profile 1 and 2 in
   `get_results(fit, "assignments")` and `predict(fit)`). Changed the
   wording in R/get-data.R:106-108 and R/predict-docs.R:34-36. Test: "units of
   sampling weight zero keep their assignment, as documented".
5. **`vcov()`/`confint()` @return (follow-up 1).** "the same parameters" as
   `coef()` was wrong for held blocks: a covariate fit with
   `fixed = "measurement"` has 14 `coef()` entries, 2 `vcov()`/`confint()`
   rows; a two-level fit 18 vs 6. R/methods-documentation.R:59-63 now says
   `coef()` reports every parameter including held ones and `vcov()`/
   `confint()` the estimated ones. Test: "vcov() and confint() cover the
   estimated parameters only; coef() every one".
6. **`"replication"` columns undocumented** (`n_screened`, `n_converged`,
   `n_failed`): R/get-data.R:149-153. Test: "the replication table has the
   columns ?get_results lists".
7. **Diagnostics `"checks"` table missing from the ?get_results catalogue**:
   R/get-data.R:359-361.
8. **Condition-registry test was vacuous** (not a package defect, but the
   guard did not guard): test-tidy-core.R's "every condition class raised in
   R/ is documented" still matches `multilpa_` and compares two empty sets.
   A working `latents_` version is in the new test file (it failed on
   `latents_single_level` until `.multilpa_notice()` raising sites were
   recognised, showing it is sensitive). test-tidy-core.R itself is not
   mine: see "Needs another slice".

## Logged, not fixed

- Follow-up 3 decided as a defect but not fixable here: the fit does not
  store the original value of a zero-weight-only ordinal value (see Needs
  another slice). The numeric→ordered-factor type of ordinal columns in
  `get_results(fit, "data")` is the established representation (simulate()
  and tests assert ordered factors); left as is.
- `descriptives(fit, by = "profile")` puts zero-weight rows in the `NA`
  stratum (`.multilpa_assignment_values()` uses the counting mask) while
  `get_results(fit, "assignments")` gives them a profile. Consistent with
  "absent from counts", inconsistent with "keep their assignment";
  test-empty-classification-consumers.R:30 asserts the current behaviour.
  Owner decision needed.
- `get_results(fit, "profiles", scale = "standardized")` divides by the
  unweighted SD over all rows, zero-weight rows included
  (`.multilpa_standardization()`); the plot uses the same helper, so table and
  plot agree. Changing it moves a reported number; owner decision.
- `"counts"` of a single-level fit has a `groups` level row (one group class,
  count = n); frozen layout.
- `latents_partial_pair` is listed in ?latents-conditions as "Not raised";
  harmless but a registry entry for a class no one can catch.
- `multilpa_indicator_reused` keeps the pre-rename prefix in code and docs
  (consistent; renaming would break callers).
- `.latents_draw_rows()` returns a length-1 `NA` for a zero-row matrix; no
  caller passes zero rows.
- `print.latents_pooled()` prints "1 group classes"; `max(fmi)` in print and
  summary is `NA` if any SE is `NA`; `.multilpa_covariate_fit_frame()` would
  fail on a pre-`extreme_logits` saved fit. Cosmetic / old objects only.

## Needs another slice

- R/extra-indicators.R:61-164 (`.latents_prepare_extra()`) and :947
  (`.latents_draw_extra_columns()`): store the sparse original values of
  zero-weight-only ordinal categories (as `.multilpa_categorical_values()`
  does with the `latents_excluded` attribute, R/categorical.R:263-267) and
  restore them, appending the value as an extra factor level. Reproduction:
  `o = 9` on a weight-0 row shows `NA` in `get_results(fit, "data")` (and is
  counted as missing by `descriptives(fit)`), while categorical `"zzz"` on a
  weight-0 row is restored. Weighted rows' values and levels are unchanged by
  such a fix.
- tests/testthat/test-tidy-core.R:269-308: the registry test matches
  `multilpa_` classes and is vacuous; switch it to `latents_` (working version
  in test-audit-1101-accessors.R).
- `devtools::document()` also rewrote man/get_candidate.Rd from another
  slice's pending roxygen in R/enumeration.R (wording of its conditions).

## Test runs

```
test-audit-1101-accessors.R: PASS 49 FAIL 0 SKIP 0 ERROR 0
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 ERROR 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 ERROR 0
test-convenience.R: PASS 86 FAIL 0 SKIP 0 ERROR 0
test-get-data.R: PASS 230 FAIL 0 SKIP 0 ERROR 0
test-predict.R: PASS 19 FAIL 0 SKIP 0 ERROR 0
test-predict-empty-rows.R: PASS 31 FAIL 0 SKIP 0 ERROR 0
test-pool-imputations.R: PASS 39 FAIL 0 SKIP 0 ERROR 0
test-descriptives-safety.R: PASS 37 FAIL 0 SKIP 0 ERROR 0
test-data.R: PASS 74 FAIL 0 SKIP 0 ERROR 0
test-data-student-esm.R: PASS 21 FAIL 0 SKIP 0 ERROR 0
test-tidy-core.R: PASS 64 FAIL 0 SKIP 0 ERROR 0
test-review-1015-accessors.R: PASS 34 FAIL 0 SKIP 0 ERROR 0
test-empty-classification-consumers.R: PASS 39 FAIL 0 SKIP 0 ERROR 0
test-release-surface.R: PASS 88 FAIL 0 SKIP 0 ERROR 0
test-review-1027-error-contract.R: PASS 216 FAIL 0 SKIP 0 ERROR 0
test-audit-services-fixes.R: PASS 31 FAIL 0 SKIP 0 ERROR 0
test-release-inference.R: PASS 104 FAIL 0 SKIP 0 ERROR 0
test-audit-inference-repairs.R: PASS 86 FAIL 0 SKIP 0 ERROR 0
```

Engine golden (`Rscript equivalence/engine-golden/compare.R`):

```
identical: 225; within_tolerance: 0; DIFFERENT: 35; ERROR: 0
```

13 are the known bootstrap attribute rename (`replicate_rows` ->
`replicate_identity`: profiles_two_level_VEI, profiles_lca_single,
profiles_count_two_level, lta_homogeneous/general/general_weights_bootstrap,
7 svc_inference_*). The other 22 (all 21 `reg_*` cases and growth_lcga) differ
only at `tables/value/group_classes@display`, "length 4 vs 5": the empty
`"group_classes"` table of a single-level mixture regression now has a fifth
column, `n_groups_assigned`, from another slice's pending change at
R/mixture-regression-methods.R:174-177. None comes from this slice's files.
