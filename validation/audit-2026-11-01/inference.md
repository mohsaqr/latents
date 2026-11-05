# Audit 2026-11-01 — slice `inference`

Files: R/inference.R, R/inference-service.R, R/structure-inference.R,
R/classification-entropy.R, R/fixed-measurement.R, R/noise.R, R/weights.R,
R/bootstrap-inference.R (126 functions). Probe scripts (synthetic data, all
run under `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`) were kept
in the session scratchpad (`inf/probe1.R` .. `probe8.R`); their checks that
guard a contract are now in `tests/testthat/test-audit-1101-inference.R`.

## Independent numerical checks (all passed; no wrong number found)

| check | result |
|---|---|
| analytic score vs `numDeriv::grad`, all 14 structures (2 profiles x 2 group classes, 3 indicators) at a displaced point | max relative error 8.9e-9 .. 1.3e-7 |
| same, weighted (two zero-weight groups) VVI/VVV/VEV | 2.4e-8 .. 3.4e-8; group scores sum to the total score to 1e-13 |
| same, mixed Gaussian + ordinal + count + categorical with FIML, weighted and unweighted; continuous local dependence (varying/equal); categorical pair under FIML; single level | 1.4e-7 or better |
| Wald `covariance_unconstrained` vs `solve(numDeriv::hessian)` (data units, free coordinates) | 9e-10 .. 1.6e-7 (VEV/EEV rotation charts 5.8e-6; FIML mixed 4.7e-5, limited by numDeriv) |
| robust = `A^-1 B A^-1 G/(G-1)` by hand, G = groups of positive weight (38 of 40) | 6e-7 .. 2.3e-6 (weighted), 8.4e-5 (FIML mixed) |
| OPG = `solve(crossprod(group scores))`; `scaling_correction = tr(A^-1 B)/q` | 0.937882 both |
| delta-method natural covariance vs `numDeriv::jacobian` of decode -> natural | 1e-11 .. 1.5e-10, all structures, weighted, mixed |
| grand centring + scaling: covariance vs a stacked M-estimation sandwich computed independently (numDeriv Jacobian of model score + mean/SD equations) | 1.9e-8 (robust), 1.8e-8 (observed) |
| direct-effect SEs after shifting the covariate by 1000 | max difference 2.1e-13 |
| held blocks (`fixed = "measurement"` / `"means"`): free count == `n_parameters`; no held row in the table; vcov vs numDeriv on free coordinates | 4.8e-5, 1.5e-4 (near-boundary logit); `confint(parm = held)` -> `latents_held_parameter` |
| entropy: `-sum p log p` with `0 log 0 = 0`; relative entropy `1 - E/(n log K)`; weighted: E and n use row weights, zero-weight rows absent (135 of 150 classified) but kept in posteriors (150 rows) | exact to 1e-12; two-group-class complete entropy equals a brute-force `H(C|y) + sum_h P(h) sum_i H(X_i|h)` |
| BH adjustment counts only the non-NA tests | equal to `p.adjust` on the 6 tests of 22 rows |
| bootstrap table: SD, percentile limits, covariance from `attr(, "replicates")`; `vcov()`/`confint()` equal the table; reordered copy refused | all equal; `latents_bad_argument` |
| noise hypervolume vs `mclust::hypvol` (1 and 2 dimensions) | 77.41161 = 77.41161; 7.741161 = 7.741161 |

## Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.latents_parameter_component` | R/inference.R:5 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.latents_parameter_decode` | R/inference.R:23 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.latents_term_component` | R/inference.R:45 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.latents_term` | R/inference.R:58 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.multilpa_parameter_names` | R/inference.R:78 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.multilpa_coefficient_labels` | R/inference.R:94 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.multilpa_coefficients` | R/inference.R:104 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.multilpa_coefficients_raw` | R/inference.R:121 | ran synthetic: via the decode/encode round trip inside every score and vcov check (probes 1-7) | ok |
| `.multilpa_response_coordinates` | R/inference.R:186 | ran synthetic: via the decode/encode round trip inside every score and vcov check (probes 1-7) | ok |
| `.multilpa_response_decode` | R/inference.R:209 | ran synthetic: via the decode/encode round trip inside every score and vcov check (probes 1-7) | ok |
| `.multilpa_response_width` | R/inference.R:229 | ran synthetic: via the decode/encode round trip inside every score and vcov check (probes 1-7) | ok |
| `.multilpa_covariance_coordinates` | R/inference.R:240 | ran synthetic: via the decode/encode round trip inside every score and vcov check (probes 1-7) | ok |
| `.multilpa_coordinate_widths` | R/inference.R:281 | ran synthetic: via the decode/encode round trip inside every score and vcov check (probes 1-7) | ok |
| `.multilpa_free_index` | R/inference.R:322 | ran synthetic: fixed means/measurement, local dependence (equal/varying), boundary fix; free count == n_parameters, vcov vs numDeriv Hessian on free coordinates (probes 2,3,7) | ok |
| `.multilpa_restore_held` | R/inference.R:367 | ran synthetic: held blocks restored in fixed fits, vcov of free part equals inverse numerical Hessian (probe 2) | ok |
| `.multilpa_coordinate_scale` | R/inference.R:395 | ran synthetic: Wald vcov equals inverse numDeriv Hessian in data units across structures; direct-effect SEs invariant to a covariate shift of 1000 (diff 2e-13, probe 2) | ok |
| `.multilpa_decode` | R/inference.R:428 | ran synthetic: via the decode/encode round trip inside every score and vcov check (probes 1-7) | ok |
| `.multilpa_response_scores` | R/inference.R:506 | ran synthetic: analytic score vs numDeriv::grad on 14 structures, weighted, mixed categorical/ordinal/count, FIML, local dependence (rel. error <= 1.3e-7); test-audit-1101-inference.R | ok |
| `.multilpa_score` | R/inference.R:554 | ran synthetic: analytic score vs numDeriv::grad on 14 structures, weighted, mixed categorical/ordinal/count, FIML, local dependence (rel. error <= 1.3e-7); test-audit-1101-inference.R | ok |
| `.multilpa_full_measurement_score` | R/inference.R:613 | ran synthetic: analytic score vs numDeriv::grad on 14 structures, weighted, mixed categorical/ordinal/count, FIML, local dependence (rel. error <= 1.3e-7); test-audit-1101-inference.R | ok |
| `.multilpa_inference_jacobian` | R/inference.R:653 | ran synthetic: J V J' vs numDeriv::jacobian of decode->natural, all 14 structures, weighted, mixed (rel <= 1.5e-10) | ok |
| `parameter_inference` | R/inference.R:1068 | read (generic); roxygen @return corrected | fixed (doc) |
| `parameter_inference.default` | R/inference.R:1079 | read; covered by test-inference.R (refusal classes) | ok |
| `parameter_inference.multilpa` | R/inference.R:1101 | ran synthetic: observed/robust/opg vs hand formulas; weighted CR1 with G = positive-weight groups; boundary fix; noise refusal; bootstrap SD/percentiles/cov by hand (probes 1-7) | fixed (weights message, via R/weights.R) |
| `.multilpa_response_bound` | R/inference.R:1388 | ran synthetic: floored category and floored last category (probe 6) | ok |
| `.multilpa_probability_floor` | R/inference.R:1424 | read; exercised by probe 6 | ok |
| `coef.multilpa` | R/inference.R:1431 | ran synthetic; covered by test-inference.R | ok |
| `vcov.multilpa` | R/inference.R:1439 | ran synthetic: Wald and bootstrap (equals replicate covariance), unconstrained bootstrap refused with class (probe 6) | ok |
| `confint.multilpa` | R/inference.R:1464 | ran synthetic: bootstrap confint equals table; held parm refused latents_held_parameter (probes 2,6); test-review-1101-confint-fixed.R | ok |
| `.multilpa_reject_inference_arguments` | R/inference.R:1525 | read; covered by test-review-1026-inference.R | ok |
| `.multilpa_check_regularity` | R/inference.R:1544 | ran synthetic: boundary categorical, NB dispersion floor, robust group count (probes 5,6) | logged |
| `.multilpa_number_single_level` | R/inference.R:1608 | ran synthetic: single-level fit inference (probe 7) | ok |
| `.multilpa_inference_matrix` | R/inference.R:1625 | ran synthetic: every probe; bad data classes in test-inference.R | ok |
| `.multilpa_validate_bootstrap_data` | R/inference.R:1742 | ran synthetic: weighted bootstrap (probe 6) | ok |
| `.multilpa_tidy_labels` | R/inference.R:1801 | ran synthetic: names/labels of every table in probes 1-7; covered by test-tidy-inference.R, test-inference.R | ok |
| `.inference_score_information` | R/inference-service.R:60 | read; shared by other families; covered by test-transition-inference.R, test-growth-parameter-inference.R | ok |
| `.inference_positive_definite` | R/inference-service.R:120 | read; covered by other families' tests | ok |
| `.inference_invert` | R/inference-service.R:161 | read | ok |
| `.inference_require_rank` | R/inference-service.R:180 | ran synthetic: OPG rank path (probe 2) | ok |
| `.multilpa_inference_groups` | R/inference-service.R:205 | ran synthetic: zero-weight groups excluded (38 of 40), empty groups under drop_empty (28 of 30) (probes 1,6) | ok |
| `.inference_cluster_factor` | R/inference-service.R:227 | ran synthetic: CR1 G/(G-1) reproduced by hand (probes 1,2) | ok |
| `.inference_sandwich` | R/inference-service.R:244 | ran synthetic: robust vcov = A^-1 B A^-1 G/(G-1) by hand (probes 1,2; rel <= 8e-5 vs numDeriv) | ok |
| `.inference_delta_covariance` | R/inference-service.R:257 | ran synthetic: delta check (probe 1) | ok |
| `.inference_delta_se` | R/inference-service.R:270 | read; used by other families | ok |
| `.inference_numeric_jacobian` | R/inference-service.R:294 | read; used by other families | ok |
| `.inference_simplex_jacobian` | R/inference-service.R:326 | ran synthetic: inside the delta check (probe 1) | ok |
| `.inference_critical` | R/inference-service.R:340 | read | ok |
| `.inference_wald` | R/inference-service.R:357 | read; used by other families | ok |
| `.latents_memo` | R/inference-service.R:386 | read; covered by test-inference-cache.R | ok |
| `.multilpa_observed_hessian` | R/inference-service.R:415 | ran synthetic: failing objective gave an unclassed error | fixed |
| `.multilpa_opg_inverse` | R/inference-service.R:465 | ran synthetic: OPG = inverse cross-product by hand (probe 2); near-singular column refused with class | ok |
| `.multilpa_cross_product` | R/inference-service.R:481 | ran synthetic (probe 2) | ok |
| `.multilpa_require_clusters` | R/inference-service.R:525 | read; covered by test-cluster-robust-cr1.R | ok |
| `.multilpa_wald_bounds` | R/inference-service.R:575 | ran synthetic: logit/log intervals in probe 6 tables | ok |
| `.multilpa_interval_kind` | R/inference-service.R:603 | read; probes | ok |
| `.multilpa_adjust_p` | R/inference-service.R:624 | ran synthetic: BH over the non-NA tests only equals p.adjust on them (probe 5) | ok |
| `.multilpa_transformed_vcov` | R/inference-service.R:659 | read; covered by test-review-1026-inference-r2.R | ok |
| `.multilpa_grand_transform` | R/inference-service.R:694 | ran synthetic: stacked M-estimation sandwich computed independently with numDeriv (probe 4, rel 1.9e-8 robust, 1.8e-8 observed) | ok |
| `.multilpa_transform_covariance` | R/inference-service.R:787 | ran synthetic: as above (probe 4) | ok |
| `.multilpa_separable_structures` | R/structure-inference.R:36 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_uses_chart` | R/structure-inference.R:42 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_spread_shared` | R/structure-inference.R:60 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_chart` | R/structure-inference.R:76 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_chart_cholesky_positions` | R/structure-inference.R:122 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_chart_rotation_pairs` | R/structure-inference.R:131 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_fitted_blocks` | R/structure-inference.R:139 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_chart_anchor` | R/structure-inference.R:162 | ran synthetic (probe 1); shared-orientation anchor with tied eigenvalues not exercised | logged |
| `.multilpa_chart_encode` | R/structure-inference.R:180 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_chart_map` | R/structure-inference.R:274 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_chart_decode` | R/structure-inference.R:389 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.multilpa_profile_gradients` | R/structure-inference.R:418 | ran synthetic: analytic score vs numDeriv::grad on 14 structures, weighted, mixed categorical/ordinal/count, FIML, local dependence (rel. error <= 1.3e-7); test-audit-1101-inference.R | ok |
| `.multilpa_chart_scores` | R/structure-inference.R:467 | ran synthetic: analytic score vs numDeriv::grad on 14 structures, weighted, mixed categorical/ordinal/count, FIML, local dependence (rel. error <= 1.3e-7); test-audit-1101-inference.R | ok |
| `.multilpa_chart_natural_jacobian` | R/structure-inference.R:486 | ran synthetic: all ten constrained structures, score/vcov/delta vs numDeriv (probe 1); covered by test-structure-inference.R | ok |
| `.latents_row_conditional_entropy` | R/classification-entropy.R:19 | ran synthetic: equals brute-force within-class entropy row by row (probe 3) | ok |
| `.multilpa_complete_entropy` | R/classification-entropy.R:56 | ran synthetic: -sum p log p (1 class), weighted with zero-weight rows, 2 group classes vs brute force (probe 3); test-audit-1101-inference.R | ok |
| `.multilpa_level_relative_entropy` | R/classification-entropy.R:82 | ran synthetic: 1 - E/(n log K), weighted n = sum of positive row weights (probe 3); test-audit-1101-inference.R | ok |
| `.multilpa_level_avepp` | R/classification-entropy.R:96 | read; covered by test-empty-classification-consumers.R | ok |
| `.multilpa_level_posteriors` | R/classification-entropy.R:108 | read | ok |
| `.multilpa_fixed_blocks` | R/fixed-measurement.R:14 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `.multilpa_validate_fixed` | R/fixed-measurement.R:29 | ran synthetic: LD + fixed variances and VEI + fixed variances refused with class (probe 3) | ok |
| `.multilpa_held_parameters` | R/fixed-measurement.R:88 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `.multilpa_apply_held` | R/fixed-measurement.R:116 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `.multilpa_fixed_parameters` | R/fixed-measurement.R:141 | ran synthetic: free count == n_parameters for fixed means/measurement (probes 2,3) | ok |
| `.multilpa_complete_start` | R/fixed-measurement.R:173 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `fit_staged` | R/fixed-measurement.R:302 | ran synthetic: n_group_classes = 1, measurement = list(), variance_model = 'bogus' raised unclassed errors | fixed |
| `.multilpa_stage_start` | R/fixed-measurement.R:426 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `.multilpa_requested_levels` | R/fixed-measurement.R:446 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `.multilpa_category_counts` | R/fixed-measurement.R:458 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `.multilpa_check_measurement` | R/fixed-measurement.R:482 | read; exercised by fixed fits in probes 2-3; covered by test-fixed-measurement.R | ok |
| `.multilpa_log_hypervolume` | R/noise.R:24 | ran synthetic: equals mclust::hypvol in 1 and 2 dimensions (probe 5) | ok |
| `.multilpa_check_noise` | R/noise.R:46 | read; noise fit refused by parameter_inference (both methods, probe 6); covered by test-noise.R | ok |
| `.multilpa_noise_start` | R/noise.R:84 | read; noise fit refused by parameter_inference (both methods, probe 6); covered by test-noise.R | ok |
| `.multilpa_add_noise` | R/noise.R:114 | read; noise fit refused by parameter_inference (both methods, probe 6); covered by test-noise.R | ok |
| `.multilpa_split_noise` | R/noise.R:127 | read; noise fit refused by parameter_inference (both methods, probe 6); covered by test-noise.R | ok |
| `.multilpa_class_posteriors` | R/noise.R:151 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_informative_rows` | R/noise.R:168 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_carrying_rows` | R/noise.R:184 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_carrying_groups` | R/noise.R:198 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_group_sample_size` | R/noise.R:210 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_check_drop_empty` | R/noise.R:220 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_informative_groups` | R/noise.R:239 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_classified_group_posteriors` | R/noise.R:256 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_classified_posteriors` | R/noise.R:273 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_reported_rows` | R/noise.R:288 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_reported_groups` | R/noise.R:302 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_reported_posteriors` | R/noise.R:314 | ran synthetic: zero-weight rows absent from classified posteriors (135 of 150) but kept in output rows (probe 3); covered by test-review-1020-zero-weights.R | ok |
| `.multilpa_refuse_noise` | R/noise.R:334 | read; noise fit refused by parameter_inference (both methods, probe 6); covered by test-noise.R | ok |
| `.latents_sampling_weights` | R/weights.R:21 | ran synthetic: weights scale to the 38/27 positive units (probes 1,3); covered by test-weights.R | ok |
| `.latents_refuse_weights` | R/weights.R:69 | read | ok |
| `.latents_weighted_vcov` | R/weights.R:95 | ran synthetic: refusal for a weighted mixture regression named method = 'bootstrap', which it has not | fixed |
| `.latents_is_weighted` | R/weights.R:116 | read | ok |
| `.latents_print_weights` | R/weights.R:129 | read; covered by test-weights.R | ok |
| `.multilpa_refit_arguments` | R/bootstrap-inference.R:31 | ran synthetic: weighted bootstrap (probe 6); covered by test-bootstrap-inference.R, test-staged-bootstrap.R | ok |
| `.multilpa_refit_measurement` | R/bootstrap-inference.R:73 | ran synthetic: weighted bootstrap (probe 6); covered by test-bootstrap-inference.R, test-staged-bootstrap.R | ok |
| `.multilpa_staged_refit_arguments` | R/bootstrap-inference.R:101 | ran synthetic: weighted bootstrap (probe 6); covered by test-bootstrap-inference.R, test-staged-bootstrap.R | ok |
| `.multilpa_bootstrap_replicates` | R/bootstrap-inference.R:126 | ran synthetic: weighted bootstrap (probe 6); covered by test-bootstrap-inference.R, test-staged-bootstrap.R | ok |
| `.multilpa_bootstrap_inference` | R/bootstrap-inference.R:196 | ran synthetic: SD, percentile limits and covariance equal hand computation from the replicates (probe 6) | ok |
| `.latents_bootstrap_table` | R/bootstrap-inference.R:275 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |
| `.latents_bootstrap_identity` | R/bootstrap-inference.R:295 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |
| `.latents_bootstrap_columns` | R/bootstrap-inference.R:307 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |
| `.latents_bootstrap_matrix` | R/bootstrap-inference.R:325 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |
| `.latents_bootstrap_replicates` | R/bootstrap-inference.R:344 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |
| `as.data.frame.latents_bootstrap` | R/bootstrap-inference.R:358 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |
| `plot.latents_bootstrap` | R/bootstrap-inference.R:371 | read; covered by test-bootstrap-replicates.R (no visual change) | ok |
| `.latents_bootstrap_spread` | R/bootstrap-inference.R:457 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |
| `.latents_bootstrap_layout` | R/bootstrap-inference.R:484 | ran synthetic: replicates accessor, reordered copy refused latents_bad_argument (probe 6); test-review-1101-bootstrap-tables.R | ok |

## Fixed

**F1. Weighted refusal named a bootstrap the caller does not have**
(coordinator item; R/weights.R:95 `.latents_weighted_vcov()`, caller
R/inference.R:1120). Before, the message for every family was the same:

```
> .latents_weighted_vcov(TRUE, "observed", FALSE)   # as mixture_regression calls it
class: latents_unsupported_weights
message: A fit with sampling weights is a pseudo likelihood fit, whose `vcov_type = "observed"` errors are not valid; use `vcov_type = "robust"` (the default for weighted fits) or `method = "bootstrap"`.
```

Cause: the helper receives only a logical, so it could not know whether its
caller offers a bootstrap; the mixture-regression and covariate paths have
none (`mixture_regression()` inference has no `method`; covariate fits refuse
`method = "bootstrap"` with `latents_unsupported_inference`; the additive
families have none). Change: an internal argument `bootstrap = FALSE`; the
sentence `or method = "bootstrap"` is added only when the caller passes
`TRUE`. `parameter_inference.multilpa()` (whose weighted bootstrap was run in
probe 6) passes `TRUE`. Class unchanged. After:

```
> vcov(<weighted mixture regression>, type = "observed")
class: latents_unsupported_weights
message: ... use `vcov_type = "robust"` (the default for weighted fits).
```

The LTA and growth-mixture callers (other slices) now get the robust-only
message; they may pass `bootstrap = TRUE` where their weighted bootstrap is
verified (see "Needs another slice"). Test: "a weighted fit's refusal names
the bootstrap only where its inference offers one" (direct helper, a weighted
`mixture_regression()` through `vcov()`, and a weighted profile fit through
`parameter_inference()`).

**F2. `.multilpa_observed_hessian()` raised unclassed errors**
(R/inference-service.R:415-460). When the likelihood cannot be evaluated a
step from the estimate, `optimHess()` fails and was re-raised with `stop()`:

```
class: simpleError, error, condition
message: Observed Hessian failed: the likelihood cannot be evaluated here
```

Every other unusable information in the Wald path is
`latents_singular_information` (documented on `parameter_inference()`), so a
caller catching that class missed this one. Change: both re-raises (the
`optimHess()` failure and the `solve()` failure) are `errorCondition(...,
class = "latents_singular_information")`; message names the cause. No number
changes. Test: "an observed information that cannot be formed or inverted is
refused with its class".

**F3. `fit_staged()` raised unclassed errors for its own arguments**
(R/fixed-measurement.R:310-325). Reproduction:

```
> fit_staged(d, c("a","b"), "school", n_profiles = 2, n_group_classes = 1)
class: simpleError    message: `n_group_classes` must be at least two; staging one leaves nothing to estimate
> fit_staged(d, c("a","b"), "school", n_profiles = 2, variance_model = "bogus")
class: simpleError    message: 'arg' should be one of "varying", "equal"
```

`multilpa()` raises `latents_bad_argument` for the same mistakes. Change:
`match.arg()` -> `.latents_match_formal()` for `variance_model`,
`covariance_model`, `missing`; `stopifnot()` -> `.latents_require()` (same
messages, class `latents_bad_argument`). Roxygen `@param measurement` now
names the class. Test: "fit_staged() classes its own argument errors"
(n_group_classes = 1, measurement = list(), variance_model, missing).

**F4. Documentation: probability floor** (R/inference.R roxygen of
`parameter_inference()`, @return). It said "A probability at or below 1e-6 is
treated as on its bound (refused, or held by `boundary = "fix"`)". Probe 8: a
profile (mixing) probability of 5.585e-07 is not refused; only categorical
response probabilities are floored at 1e-6 (`.multilpa_check_regularity()`),
and `boundary = "fix"` holds response probabilities only. Now: "A categorical
response probability at or below 1e-6 is treated as on its bound ...; a
mixing probability is refused only at zero." (man/parameter_inference.Rd
regenerated.)

Each fix was broken on purpose once (message always naming the bootstrap;
unclassed `stop()` restored; `match.arg()`/`stopifnot()` restored): the four
corresponding tests failed (2 failures, 2 errors), then the files were
restored.

## Logged, not fixed

- Mixing probability near zero (5.6e-7 in probe 8) passes the regularity check; its Wald row is honest but uninformative (SE 0.0017, logit interval [0, 1]). Refusing it would remove numbers from existing valid fits (frozen); doc corrected instead (F4).
- OPG on that same fit is refused (`latents_singular_information`) because the near-zero logit's score column is ~1e-11: correct behaviour, recorded for context.
- `.multilpa_chart_anchor()` for a shared orientation takes eigenvectors of a weighted sum of the profile covariances; an accidental tie in that sum's eigenvalues (not in the profiles') would give an arbitrary basis. Not reachable on realistic data; not exercised.
- `.multilpa_inference_matrix()` centres with `colMeans(x, na.rm = TRUE)` over all rows including zero-weight ones; a translation only, harmless except for junk magnitudes (LEARNINGS 2026-10-22, 1e20 row).
- Wald covariance vs numDeriv agrees to 5e-5 / 1.5e-4 (relative, max entry) on FIML mixed fits and a fit with a logit near 14; attributed to numDeriv's Richardson error near flat directions, not to optimHess (scores agree to 1e-7).
- `parameter_inference()` @return does not list the bootstrap attributes `replicate_identity` and `adjust` or the Wald attributes `hessian_conditioned`, `group_scores_conditioned`; internal plumbing, deliberately not advertised.
- `.multilpa_log_hypervolume()` gives `-Inf` for a constant indicator (zero range), as `mclust::hypvol()` gives 0; a constant indicator is refused at fit time.

## Needs another slice

- R/simulation-service.R:141 `.latents_cluster_bootstrap()`: draws `n_groups` groups from all groups, including zero-weight (and, under `drop_empty`, empty) groups, so the number of sampled units that count varies across resamples (Binomial(G, G+/G) instead of G+). Under the zero-weight rule the resample should draw G+ units from the positive-weight ones. Affects bootstrap SEs of weighted fits with zero-weight units only.
- R/lta-methods.R:408, 553, 592 and R/growth-mixture-methods.R:209, R/mixture-regression.R:417 (growth path, `random`): callers of `.latents_weighted_vcov()`; where a weighted bootstrap exists and is verified (LTA; growth mixture except multilevel), pass `bootstrap = TRUE` to name it in the refusal. Mixture regression without `random`, covariate (R/covariate-inference.R:361, 649), additive (R/additive-methods.R) and pool-imputations callers must keep the default `FALSE`.

## Test runs (verbatim)

test-audit-1101-inference.R: PASS 25 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-inference.R: PASS 54 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-inference-cache.R: PASS 21 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-inference-fixed-measurement.R: PASS 33 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-structure-inference.R: PASS 42 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-fixed-measurement.R: PASS 71 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-noise.R: PASS 17 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-weights.R: PASS 76 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-cluster-robust-cr1.R: PASS 6 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-bootstrap-inference.R: PASS 66 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-bootstrap-replicates.R: PASS 50 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-review-1101-bootstrap-tables.R: PASS 67 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-review-1101-confint-fixed.R: PASS 77 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-audit-inference-repairs.R: PASS 86 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-release-inference.R: PASS 104 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-review-families.R: PASS 72 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-staged-bootstrap.R: PASS 57 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-staged-identity.R: PASS 54 FAIL 0 SKIP 0 ERROR 0 WARN 0
test-audit-1101-inference.R before the VVV tolerance was set to 1e-5 (numDeriv error 6.3e-6 at that displaced point): 1 failure; after: PASS 25.

Baseline before any change (same files): every file PASS with 0 FAIL / 0 SKIP / 0 ERROR (test-inference 54, -cache 21, -fixed-measurement 33, structure 42, fixed-measurement 71, noise 17, weights 76, cr1 6, bootstrap-inference 66, bootstrap-replicates 50, review-1101-bootstrap-tables 67, review-1101-confint-fixed 77, audit-inference-repairs 86, release-inference 104, api-compat-0921 44, compat-0921 10).

## devtools::document()

Rewrote man/parameter_inference.Rd and man/fit_staged.Rd (this slice's
roxygen) and man/diagnostics.Rd (roxygen changed by another slice in
R/diagnostics.R). NAMESPACE unchanged.

## Engine golden (`Rscript equivalence/engine-golden/compare.R`)

20 cases differ (this run rewrote equivalence/engine-golden/last-compare.csv):

- 13 bootstrap cases, `replicate_rows` -> `replicate_identity` attribute (expected, F2/F9 of the bootstrap review): profiles_two_level_VEI, profiles_lca_single, profiles_count_two_level, lta_homogeneous_bootstrap, lta_general_bootstrap, lta_general_weights_bootstrap, svc_inference_lpa, svc_inference_lpa_full, svc_inference_two_level, svc_inference_bootstrap_lca, svc_inference_bootstrap_mixed_extra, svc_inference_bootstrap_vei, svc_inference_bootstrap_staged.
- 6 from F1, refusal message text only (class unchanged): cov_weights, family_additive_weights, lta_general_weights, lta_general_mixed_extensions, reg_gaussian_weights_observed_refused, growth_weights (`/probes/inference_observed/value/message` or `/refusal/message`).
- 1 from F3, error class `simpleError` -> `latents_bad_argument`: svc_fit_staged_one_class (`/refusal/error_class`).
- No reg_*/growth_lcga `group_classes@display` difference appeared in this run; no other difference.
