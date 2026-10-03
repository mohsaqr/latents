# Engine golden baseline (phase E0)

Phase 0 of the engine consolidation (`validation/ENGINE_DESIGN.md`). Before any
estimation engine is moved onto the shared kernel, this directory records what
every engine and documented option produces **today**, so each later phase can
show that nothing moved.

A *case* fits one model (or runs one service verb) on a small slice of a
bundled dataset with a fixed seed. Its *fingerprint* is everything a user can
observe about the result, reduced to plain data:

- `print(x)` and `print(summary(x))` as text;
- every table: `get_results(x, "all")` for every fit class (and again with
  `vcov_type = "robust"` for the classes whose `get_results()` takes it:
  additive, `lta()` general engine, mixture regression, growth mixture);
  `"candidates"`/`"best"` for family and transition enumerations, `"fit"` for
  regression enumerations, `"estimates"`/`"imputations"`/`"fits"` for pooled
  imputations;
- `as.data.frame(x)`, `logLik` (value, df, nobs), `nobs`, `coef`, `vcov`,
  `confint` where the class has a method;
- `parameter_inference(x, vcov_type = "observed")` and `"robust"` where the
  class has a method (the table with all its attributes: covariance, Hessian,
  gradient, ...);
- `simulate(x, seed = 1)` where the class has a method;
- extra *service* probes named by the case (diagnostics, three-step, R3STEP,
  bootstrap inference, bootstrap LRT, OPG, `boundary = "fix"`, predict,
  starting values, report, TNA conversion, ...), each fingerprinted the same
  way;
- every warning and message raised while fitting and while probing (class and
  text);
- a refusal: when the package refuses a fit or a probe, the condition class
  and message are recorded instead. A refusal is behaviour to preserve.

Errors that are *not* one of the package's classed refusals
(`latents_*`/`multilpa_*`) are recorded too, marked `expected = FALSE`, and
listed in `unexpected-errors.csv` and in the record run's output.

Functions, environments and plot objects are replaced by placeholders and
formulas/calls are stored deparsed, so an RDS compares equal across sessions.
Each probe starts from `set.seed(<case seed>)`, so no probe depends on the
random draws of an earlier one. The case seed is derived from the case name
(`golden_case_seed()`); the fits also pass explicit `seed =` arguments.

## Files

| File | Role |
|---|---|
| `catalogue.R` | `golden_case()`, shared data slices, `golden_catalogue()` joining the sections |
| `catalogue-profiles.R` | `lpa()`, `multilpa()`, `lca()`, `multilca()`: structures, indicator types, weights, noise, prior, fixed, centering, sequences, acceleration |
| `catalogue-covariates-families.R` | membership covariates; Houle families (additive, dispersion, cross-level); family enumeration and BLRT |
| `catalogue-transitions.R` | `lta()`: homogeneous and general engines, orders, covariates, mover-stayer, indicators; LTA bootstrap, BLRT, enumeration |
| `catalogue-regression.R` | `mixture_regression()` and growth mixtures; `enumerate_regressions()`, `compare_models()` |
| `catalogue-services.R` | enumeration, BLRT, inference, three-step, R3STEP, sensitivity, pooling, staged fits, starting values, diagnostics, predict, descriptives |
| `fingerprint.R` | `golden_fingerprint()`, `golden_probe()`, `golden_services()`, `golden_compare()` |
| `setup.R` | loads the package (`devtools::load_all()`), pins session settings, runs one case |
| `record.R` | writes the baseline |
| `compare.R` | re-fits and compares against the baseline |
| `golden/<case>.rds` | one fingerprint per case |
| `timing.csv` | fit time per case (`case`, `seconds`) |
| `catalogue.csv` | case, engine, fitted/refused, services, what it covers |
| `unexpected-errors.csv` | non-refusal errors seen while recording |
| `session-info.txt` | package version and `sessionInfo()` of the recording |
| `last-compare.csv` | the most recent comparison |

## Usage

From the package root:

```sh
Rscript equivalence/engine-golden/record.R              # record every case
Rscript equivalence/engine-golden/record.R '^lta_'      # re-record matching cases only
Rscript equivalence/engine-golden/compare.R             # compare every case
Rscript equivalence/engine-golden/compare.R '^growth_'  # compare matching cases
```

`compare.R` prints one row per case: `case`, `status`, `max_rel_diff`,
`first_diff_path`, `seconds_now`, `seconds_baseline`, `slowdown`, `slow`, and
writes `last-compare.csv` (which also has a `detail` column with the first
differing values). It exits with status 1 when any case is `DIFFERENT` or
`ERROR`.

Setting `GOLDEN_OUT=/some/dir` makes either script read and write there
instead of here (a scratch baseline that leaves the committed one alone).

Re-record only for an *announced, intended* change of results, and only the
cases that change; the diff of `golden/` then shows exactly what moved.

## Statuses and tolerances

| Status | Meaning |
|---|---|
| `identical` | the fingerprint is bit-for-bit the baseline |
| `within_tolerance` | only doubles differ, all within tolerance |
| `DIFFERENT` | anything else differs (fails) |
| `ERROR` | no baseline, or the harness itself failed for the case (fails) |

- Character, logical and integer values, names, classes, dimensions and every
  attribute must be **identical**. This includes the printed output and the
  text of warnings and refusals.
- Doubles: relative difference `|a - b| / max(|a|, |b|)` at most **1e-10** for
  any path naming a log likelihood (`log_lik`, `loglik`, `log_likelihood`),
  **1e-8** for everything else. An absolute difference of at most **1e-10**
  always passes (values at or near zero). `NA`/`NaN`/`Inf` patterns must match
  exactly.
- Timing is a warning, not a failure: `slow` is `TRUE` when the fit is more
  than **1.2x** its baseline time *and* more than 0.05 s slower. Only the fit
  is timed (not the probes); a fit under one second is timed twice and the
  shorter time kept.

Caveats for later phases. Printed output rounds to about seven significant
digits, so an estimate moving within tolerance can still flip a printed digit
and turn the case `DIFFERENT` through its `print` path; check `detail` and
judge. The `gradient`, `scaled_score` and `score_displacement` attributes of
inference tables are numerically near zero at a maximum, so their relative
differences are large under any change to how scores are summed; a difference
confined to those paths is not a change in results.

## Session pinning

`setup.R` sets `options(width = 80, digits = 7, scipen = 0, OutDec = ".")`,
`LC_COLLATE = "C"`, the default `RNGkind()`, and sends plots to a null device.
The baseline was recorded on R 4.5.2 / macOS (see `session-info.txt`); BLAS,
platform or R version changes can move the last bits of a fit, so compare on
the platform that recorded.

## Cases

251 cases recorded on 0.9.12 (213 fitted, 38 recorded refusals), plus 9 added in phase E7 for the new models (multilevel growth, negative-binomial and ordinal regression mixtures): 260 in all. Services in brackets.

| Case | Engine | Result | Covers |
|---|---|---|---|
| `profiles_two_level_basic` | multilpa | fitted | multilpa two-level, 2 profiles x 2 group classes, diagonal varying [diagnostics three_step_bch bootstrap_lrt] |
| `profiles_refusal_noise_two_level` | multilpa | refused | noise = TRUE with two group classes is refused |
| `profiles_lpa_single_varying` | multilpa | fitted | lpa() single-level srl (all 5 indicators), 3 profiles, varying diagonal [descriptives predict starting_values inference_opg report] |
| `profiles_lpa_single_equal` | multilpa | fitted | lpa() single-level, variance_model = "equal" |
| `profiles_lpa_single_full` | multilpa | fitted | lpa() single-level, covariance_model = "full", varying |
| `profiles_lpa_single_full_equal` | multilpa | fitted | lpa() single-level, full covariance shared across profiles |
| `profiles_lpa_engagement_id_null` | multilpa | fitted | multilpa(id = NULL) on engagement rows: single-level message, 2 profiles |
| `profiles_lpa_one_profile` | multilpa | fitted | lpa() with a single profile |
| `profiles_two_level_equal` | multilpa | fitted | multilpa two-level, variance_model = "equal" [predict starting_values] |
| `profiles_two_level_full` | multilpa | fitted | multilpa two-level, covariance_model = "full" (40 students, so robust/OPG apply) [inference_opg] |
| `profiles_two_level_three_profiles` | multilpa | fitted | multilpa two-level, 3 profiles x 2 group classes |
| `profiles_two_level_one_class` | multilpa | fitted | multilpa two-level data with one group class [starting_values_measurement diagnostics] |
| `profiles_structure_EII` | multilpa | fitted | lpa() single-level, 2 profiles, model = "EII" (srl, 150 rows) |
| `profiles_structure_VII` | multilpa | fitted | lpa() single-level, 2 profiles, model = "VII" (srl, 150 rows) |
| `profiles_structure_EEI` | multilpa | fitted | lpa() single-level, 2 profiles, model = "EEI" (srl, 150 rows) |
| `profiles_structure_VEI` | multilpa | fitted | lpa() single-level, 2 profiles, model = "VEI" (srl, 150 rows) |
| `profiles_structure_EVI` | multilpa | fitted | lpa() single-level, 2 profiles, model = "EVI" (srl, 150 rows) |
| `profiles_structure_VVI` | multilpa | fitted | lpa() single-level, 2 profiles, model = "VVI" (srl, 150 rows) |
| `profiles_structure_EEE` | multilpa | fitted | lpa() single-level, 2 profiles, model = "EEE" (srl, 150 rows) |
| `profiles_structure_VEE` | multilpa | fitted | lpa() single-level, 2 profiles, model = "VEE" (srl, 150 rows) |
| `profiles_structure_EVE` | multilpa | fitted | lpa() single-level, 2 profiles, model = "EVE" (srl, 150 rows) |
| `profiles_structure_VVE` | multilpa | fitted | lpa() single-level, 2 profiles, model = "VVE" (srl, 150 rows) |
| `profiles_structure_EEV` | multilpa | fitted | lpa() single-level, 2 profiles, model = "EEV" (srl, 150 rows) |
| `profiles_structure_VEV` | multilpa | fitted | lpa() single-level, 2 profiles, model = "VEV" (srl, 150 rows) |
| `profiles_structure_EVV` | multilpa | fitted | lpa() single-level, 2 profiles, model = "EVV" (srl, 150 rows) |
| `profiles_structure_VVV` | multilpa | fitted | lpa() single-level, 2 profiles, model = "VVV" (srl, 150 rows) |
| `profiles_structure_pieces_VEV` | multilpa | fitted | structure by volume/shape/orientation (varying/equal/varying = VEV) |
| `profiles_structure_pieces_EII` | multilpa | fitted | structure by pieces: equal volume, spherical shape (EII) |
| `profiles_structure_three_profiles_EVE` | multilpa | fitted | EVE (minorize-maximize orientation step) with 3 profiles |
| `profiles_two_level_VEI` | multilpa | fitted | two-level multilpa with model = "VEI" [inference_bootstrap] |
| `profiles_two_level_EEV` | multilpa | fitted | two-level multilpa with model = "EEV" |
| `profiles_two_level_VVE` | multilpa | fitted | two-level multilpa with model = "VVE" |
| `profiles_refusal_model_and_pieces` | multilpa | refused | model = together with volume/shape/orientation is refused |
| `profiles_fiml_two_level` | multilpa | fitted | missing = "fiml", two-level diagonal [descriptives] |
| `profiles_fiml_single_full` | multilpa | fitted | missing = "fiml", single-level full covariance |
| `profiles_refusal_missing_error` | multilpa | refused | missing values under missing = "error" are refused |
| `profiles_lca_single` | multilca | fitted | lca() single-level on 4 binary student_esm items, 2 classes [predict starting_values inference_bootstrap] |
| `profiles_multilca_two_level` | multilca | fitted | multilca() two-level on binary items, 2 classes x 2 group classes [diagnostics] |
| `profiles_mixed_categorical` | multilpa | fitted | multilpa mixed: 2 continuous (happy, relaxed) + 2 categorical items |
| `profiles_ordinal_single` | multilpa | fitted | lpa() with ordinal items (worried, exhausted) + continuous happy |
| `profiles_ordinal_two_level` | multilpa | fitted | multilpa two-level with an ordinal item and a categorical item |
| `profiles_count_poisson` | multilpa | fitted | lpa() with a Poisson count (questions) + continuous score [predict] |
| `profiles_count_nb_varying` | multilpa | fitted | negative binomial count, count_dispersion = "varying" |
| `profiles_count_nb_equal` | multilpa | fitted | negative binomial count, count_dispersion = "equal" |
| `profiles_count_two_level` | multilpa | fitted | multilpa two-level: Poisson count + ordinal + continuous (mixed) [inference_bootstrap] |
| `profiles_mixed_all_types` | multilpa | fitted | single-level: continuous + categorical + ordinal + count |
| `profiles_extra_fiml` | multilpa | fitted | ordinal + count indicators with missing = "fiml" |
| `profiles_refusal_ordinal_fixed` | multilpa | refused | ordinal indicators with fixed = measurement are refused |
| `profiles_weights_single` | multilpa | fitted | lpa() with row sampling weights |
| `profiles_weights_two_level` | multilpa | fitted | multilpa two-level with group sampling weights |
| `profiles_weights_categorical` | multilpa | fitted | weighted two-level fit with categorical indicators |
| `profiles_refusal_weights_prior` | multilpa | refused | weights together with prior are refused |
| `profiles_noise_single` | multilpa | fitted | noise = TRUE single-level, VVI |
| `profiles_noise_VVV` | multilpa | fitted | noise = TRUE single-level, model = "VVV" |
| `profiles_prior_VVI` | multilpa | fitted | prior = prior_control(), VVI |
| `profiles_prior_EEE` | multilpa | fitted | prior = prior_control(), EEE |
| `profiles_prior_VEV_shrinkage` | multilpa | fitted | prior_control(shrinkage = 0.1), VEV |
| `profiles_prior_noise_EII` | multilpa | fitted | prior and noise together, EII |
| `profiles_refusal_prior_VEE` | multilpa | refused | prior with VEE (no mclust prior) is refused |
| `profiles_fixed_measurement` | multilpa | fitted | fixed = "measurement" from a one-class stage, 2 group classes [inference_fix bootstrap_lrt] |
| `profiles_fixed_means` | multilpa | fitted | fixed = "means" only |
| `profiles_fixed_categorical` | multilpa | fitted | fixed = "response_probabilities" on a latent class model |
| `profiles_evaluate_only` | multilpa | fitted | max_iter = 0 with start: evaluate a supplied parameter set |
| `profiles_user_start` | multilpa | fitted | user start list (means, variances, probabilities), two-level |
| `profiles_refusal_fixed_without_start` | multilpa | refused | fixed = "means" without start is refused |
| `profiles_centering_person` | multilpa | fitted | centering = "person", two-level |
| `profiles_centering_grand` | multilpa | fitted | centering = "grand", two-level |
| `profiles_time_sequences` | multilpa | fitted | time = "sequence": sequence tables, get_tna / get_group_tna [get_tna get_group_tna] |
| `profiles_acceleration_none` | multilpa | fitted | acceleration = "none" (plain EM), two-level |
| `profiles_acceleration_squarem_full` | multilpa | fitted | acceleration = "squarem" explicit, single-level full covariance |
| `profiles_acceleration_none_full` | multilpa | fitted | acceleration = "none", single-level full covariance |
| `profiles_select_converged` | multilpa | fitted | select_start = "converged": max_iter = 60 leaves the best start unconverged, a worse one converged |
| `profiles_select_likelihood` | multilpa | fitted | select_start = "likelihood" on the same short run (picks the unconverged best) |
| `profiles_min_variance_tol` | multilpa | fitted | non-default tol and min_variance |
| `profiles_refusal_two_classes_no_id` | multilpa | refused | n_group_classes > 1 with id = NULL is refused |
| `profiles_refusal_missing_id` | multilpa | refused | omitting id is refused |
| `profiles_refusal_bad_model` | multilpa | refused | an unknown model code is refused |
| `cov_profile_shared` | covariates | fitted | two-level, profile_covariates only (row covariate), shared slopes [inference_opg diagnostics inference_bootstrap bootstrap_lrt] |
| `cov_group_only` | covariates | fitted | two-level, group_covariates only (student-constant covariate) [inference_opg diagnostics] |
| `cov_both_shared` | covariates | fitted | two-level, profile and group covariates, shared slopes [inference_opg diagnostics] |
| `cov_slopes_group_class` | covariates | fitted | two-level, profile_slopes = group_class (cross-level interaction) [inference_opg] |
| `cov_one_group_class` | covariates | fitted | two-level with one group class, profile covariate (pooled LPA logit) [inference_opg diagnostics] |
| `cov_single_level` | covariates | fitted | single level (lpa, id = NULL) with a profile covariate [inference_opg diagnostics] |
| `cov_three_profiles_equal` | covariates | fitted | three profiles, variance_model equal, profile covariate |
| `cov_fiml` | covariates | fitted | profile covariate with missing indicators, missing = fiml [inference_opg] |
| `cov_full_covariance` | covariates | fitted | profile covariate, covariance_model = full [inference_opg] |
| `cov_categorical_esm` | covariates | fitted | two continuous + binary categorical (ESM) with a row covariate [inference_opg inference_bootstrap] |
| `cov_categorical_boundary` | covariates | fitted | categorical indicator at its probability bound with a covariate (Wald refused) [inference_opg] |
| `cov_ordinal` | covariates | fitted | ordinal indicator (adjacent-category logit) with a profile covariate |
| `cov_count_poisson` | covariates | fitted | Poisson count indicator with a profile covariate |
| `cov_esm_mixed` | covariates | fitted | ESM: gaussian + ordinal + binary categorical, row covariate day, two levels |
| `cov_weights` | covariates | fitted | sampling weights (constant within student) with profile covariate |
| `cov_refusal_start` | covariates | refused | refusal: start = with membership covariates |
| `cov_refusal_fixed` | covariates | refused | refusal: fixed measurement with membership covariates |
| `cov_refusal_missing_covariate` | covariates | refused | refusal: missing values in a membership covariate |
| `cov_refusal_slopes_without_covariates` | covariates | refused | refusal: profile_slopes = group_class without profile_covariates |
| `cov_refusal_group_covariate_varies` | covariates | refused | refusal: a group covariate that varies within a group |
| `family_additive_varying` | additive | fitted | family additive, between_variance varying, 2 group classes (simulated, interior) [inference_opg inference_bootstrap inference_fix] |
| `family_additive_equal` | additive | fitted | family additive, between_variance equal [inference_opg] |
| `family_additive_three_classes` | additive | fitted | family additive, 3 group classes (a between variance at zero: boundary) |
| `family_additive_weights` | additive | fitted | family additive with sampling weights (robust default) |
| `family_additive_engagement_boundary` | additive | fitted | family additive on course_engagement: between variance at zero, Wald refused [inference_fix] |
| `family_dispersion_equal` | additive | fitted | family dispersion (between_variance equal, the only choice) [inference_opg] |
| `family_dispersion_three_classes` | additive | fitted | family dispersion, 3 group classes (singular information) |
| `family_dispersion_varying_refused` | additive | refused | refusal: family dispersion with between_variance varying |
| `family_additive_dispersion_varying` | additive | fitted | family additive_dispersion, between_variance varying [inference_opg] |
| `family_additive_dispersion_equal` | additive | fitted | family additive_dispersion, between_variance equal |
| `family_additive_dispersion_three_classes` | additive | fitted | family additive_dispersion, varying, 3 group classes (interior) |
| `family_restricted_cross_level_varying` | cross_level | fitted | family restricted_cross_level, between varying, profiles varying |
| `family_restricted_cross_level_equal` | cross_level | fitted | family restricted_cross_level, between equal, profiles equal |
| `family_full_cross_level_varying` | cross_level | fitted | family full_cross_level, between varying, profiles varying |
| `family_full_cross_level_equal` | cross_level | fitted | family full_cross_level, between equal, profiles equal |
| `family_refusal_categorical` | additive | refused | refusal: additive family with a categorical indicator |
| `family_refusal_fiml` | additive | refused | refusal: additive family with missing = fiml |
| `family_refusal_covariates` | additive | refused | refusal: additive family with membership covariates |
| `family_refusal_cross_level_categorical` | cross_level | refused | refusal: full_cross_level with a categorical indicator |
| `family_enumeration` | additive | fitted | enumerate_classes over additive/dispersion/additive_dispersion x 1:2 classes |
| `family_bootstrap_lrt` | additive | fitted | bootstrap_lrt additive 1 vs 2 group classes (iter 5, seeded) |
| `family_bootstrap_lrt_refused` | additive | refused | refusal: bootstrap_lrt between additive and dispersion (not nested) |
| `lta_homogeneous_basic` | lta (homogeneous) | fitted | lta() 2 profiles, 1 group class, diagonal varying, 20 students x ~14 courses [diagnostics descriptives inference_opg get_tna get_group_tna inference_bootstrap] |
| `lta_homogeneous_group_classes` | lta (homogeneous) | fitted | lta() 2 profiles x 2 group classes (mixture of transition matrices) [get_tna get_group_tna diagnostics] |
| `lta_homogeneous_equal_variance` | lta (homogeneous) | fitted | lta() variance_model = equal |
| `lta_homogeneous_full_covariance` | lta (homogeneous) | fitted | lta() covariance_model = full (varying) [inference_opg] |
| `lta_homogeneous_model_EEI` | lta (homogeneous) | fitted | lta() model = EEI |
| `lta_homogeneous_model_VVI` | lta (homogeneous) | fitted | lta() model = VVI |
| `lta_homogeneous_model_EEE` | lta (homogeneous) | fitted | lta() model = EEE |
| `lta_homogeneous_model_VVV` | lta (homogeneous) | fitted | lta() model = VVV |
| `lta_homogeneous_model_VEI` | lta (homogeneous) | fitted | lta() model = VEI (inference refused for this structure) |
| `lta_homogeneous_categorical` | lta (homogeneous) | fitted | lta() categorical binary indicators (ESM activities) [inference_fix] |
| `lta_homogeneous_mixed` | lta (homogeneous) | fitted | lta() mixed indicators: two Gaussian, one categorical factor |
| `lta_homogeneous_fiml` | lta (homogeneous) | fitted | lta() missing = fiml with missing indicator values |
| `lta_homogeneous_occasions_observed` | lta (homogeneous) | fitted | lta() occasions = observed on a panel with skipped positions |
| `lta_homogeneous_occasions_grid` | lta (homogeneous) | fitted | lta() occasions = grid on a panel with skipped positions |
| `lta_homogeneous_select_converged` | lta (homogeneous) | fitted | lta() select_start = converged, 3 profiles, 2 group classes [inference_fix] |
| `lta_homogeneous_bootstrap` | lta (homogeneous) | fitted | parameter_inference(method = bootstrap, iter = 5) on a homogeneous fit |
| `lta_general_occasion_transitions` | lta (general) | fitted | lta() transitions = occasion (one matrix per move) [inference_opg] |
| `lta_general_order2` | lta (general) | fitted | lta() order = 2 (second-order Markov chain) [inference_fix] |
| `lta_general_mover_stayer` | lta (general) | fitted | lta() mover_stayer = TRUE |
| `lta_general_transition_covariates` | lta (general) | fitted | lta() transition_covariates = previous_grade [inference_opg get_tna] |
| `lta_general_initial_covariates` | lta (general) | fitted | lta() initial_covariates = previous_grade |
| `lta_general_occasion_measurement` | lta (general) | fitted | lta() measurement = occasion (non-invariant profiles) |
| `lta_general_group_classes` | lta (general) | fitted | lta() n_group_classes = 2 with transition covariates |
| `lta_general_full_covariance` | lta (general) | fitted | lta() covariance_model = full with occasion transitions |
| `lta_general_model_VEI` | lta (general) | fitted | lta() model = VEI with transition covariates (EM, inference refused) |
| `lta_general_categorical` | lta (general) | fitted | lta() categorical indicators with occasion transitions |
| `lta_general_ordinal` | lta (general) | fitted | lta() ordinal indicator (adjacent-category logit) with a Gaussian one |
| `lta_general_count_poisson` | lta (general) | fitted | lta() Poisson count indicator with a Gaussian one |
| `lta_general_count_nb_varying` | lta (general) | fitted | lta() negative binomial count, count_dispersion = varying |
| `lta_general_count_nb_equal` | lta (general) | fitted | lta() negative binomial count, count_dispersion = equal |
| `lta_general_weights` | lta (general) | fitted | lta() weights (pseudo ML; robust default, observed refused) |
| `lta_general_fiml_covariates` | lta (general) | fitted | lta() missing = fiml with transition covariates |
| `lta_general_mixed_extensions` | lta (general) | fitted | lta() combination: mixed categorical/ordinal/count indicators, weights, fiml, occasion transitions |
| `lta_general_order2_mover_stayer` | lta (general) | fitted | lta() combination: order 2, mover-stayer, occasion transitions |
| `lta_general_bootstrap` | lta (general) | fitted | parameter_inference(method = bootstrap, iter = 5) on a covariate transition fit |
| `lta_general_weights_bootstrap` | lta (general) | fitted | parameter_inference(method = bootstrap) on a weighted transition fit |
| `lta_refusal_one_profile` | lta (refusal) | refused | lta() n_profiles = 1 is refused (latents_bad_transition) |
| `lta_refusal_order3` | lta (refusal) | refused | lta() order = 3 is refused |
| `lta_refusal_no_time` | lta (refusal) | refused | lta() without time is refused |
| `lta_refusal_missing_error` | lta (refusal) | refused | lta() missing indicators under missing = error are refused |
| `lta_refusal_bad_model` | lta (refusal) | refused | lta() an unknown covariance structure code is refused |
| `lta_bootstrap_lrt_covariates` | lta (compare) | fitted | bootstrap_lrt() homogeneous null vs transition-covariate alternative |
| `lta_bootstrap_lrt_group_classes` | lta (compare) | fitted | bootstrap_lrt() 1 vs 2 group classes of homogeneous transition fits |
| `lta_refusal_bootstrap_lrt_reversed` | lta (compare) | refused | bootstrap_lrt() with the models in the wrong order is refused |
| `lta_enumeration` | lta (compare) | fitted | enumerate_classes(time =) over 2:3 profiles x 1:2 group classes |
| `lta_enumeration_structures` | lta (compare) | fitted | enumerate_classes(time =) crossing structures VVI and VEI |
| `lta_enumeration_candidate` | lta (compare) | fitted | candidate_fit() picks the 3-profile, 1-class fit out of a transition grid |
| `lta_refusal_enumeration_one_profile` | lta (compare) | refused | enumerate_classes(time =) with a one-profile candidate is refused |
| `reg_gaussian_single` | mixture_regression | fitted | gaussian, single level, 2 classes, varying variance, observed vcov; predict [predict] |
| `reg_gaussian_equal_variance` | mixture_regression | fitted | gaussian, variance = equal |
| `reg_gaussian_common` | mixture_regression | fitted | gaussian, common = ~ sleep (coefficient shared across classes) [predict] |
| `reg_gaussian_membership` | mixture_regression | fitted | gaussian, row-level membership = ~ sleep (concomitant logit) [predict] |
| `reg_gaussian_select_converged` | mixture_regression | fitted | gaussian, select_start = converged |
| `reg_gaussian_vcov_robust_clustered` | mixture_regression | fitted | gaussian, observation classes with id and one group class: robust vcov clustered on id |
| `reg_gaussian_vcov_opg` | mixture_regression | fitted | gaussian, vcov_type = opg |
| `reg_gaussian_vcov_none` | mixture_regression | fitted | gaussian, vcov_type = none (no inference stored) |
| `reg_gaussian_weights_rows` | mixture_regression | fitted | gaussian, single-level sampling weights per row (robust default) |
| `reg_gaussian_weights_observed_refused` | mixture_regression | refused | weights with vcov_type = observed is refused |
| `reg_gaussian_missing_omit` | mixture_regression | fitted | missing = omit drops rows with NA outcome (latents_rows_dropped warning) |
| `reg_refusal_missing_error` | mixture_regression | refused | missing = error with NA outcome is refused (latents_missing_data) |
| `reg_group_level` | mixture_regression | fitted | gaussian, class_level = group (one class per student) [predict] |
| `reg_group_membership` | mixture_regression | fitted | class_level = group with group-constant membership = ~ motivation |
| `reg_group_weights` | mixture_regression | fitted | class_level = group with per-student sampling weights |
| `reg_two_level` | mixture_regression | fitted | two-level: observation classes, id, n_group_classes = 2 [predict] |
| `reg_two_level_group_membership` | mixture_regression | fitted | two-level with group_membership = ~ motivation and variance = equal |
| `reg_two_level_robust` | mixture_regression | fitted | two-level, vcov_type = robust |
| `reg_binomial_group` | mixture_regression | fitted | binomial (0/1 outcome), class_level = group [predict] |
| `reg_binomial_cbind` | mixture_regression | fitted | binomial with cbind(successes, failures), single level [predict] |
| `reg_poisson` | mixture_regression | fitted | poisson, single level [predict] |
| `reg_poisson_offset` | mixture_regression | fitted | poisson with offset(log(exposure)) and class_level = group |
| `reg_refusal_binary_rows` | mixture_regression | refused | binary outcome with one trial per class assignment is refused |
| `reg_refusal_group_without_id` | mixture_regression | refused | class_level = group without id is refused |
| `reg_refusal_poisson_noninteger` | mixture_regression | refused | poisson on a non-count outcome is refused |
| `reg_enumerate` | enumerate_regressions | fitted | enumerate_regressions over 1:2 classes, no bootstrap |
| `reg_enumerate_bootstrap` | enumerate_regressions | fitted | enumerate_regressions with a 3-replicate BLRT (bootstrap_starts = 2) |
| `reg_compare_models` | compare_models | fitted | compare_models over 1-, 2- and 3-class single-level regressions |
| `growth_lcga` | mixture_regression | fitted | latent class growth analysis: score ~ wave, class_level = group, no random [predict] |
| `growth_random_intercept` | growth_mixture | fitted | growth mixture, random = intercept, covariance varying [predict] |
| `growth_random_slope` | growth_mixture | fitted | growth mixture, random intercept + slope (random = wave), varying; compare_models vs LCGA [compare_models] |
| `growth_random_slope_only` | growth_mixture | fitted | growth mixture, random = ~ 0 + wave (slope without random intercept) |
| `growth_covariance_equal` | growth_mixture | fitted | growth mixture, random = wave, random_covariance = equal |
| `growth_covariance_proportional` | growth_mixture | fitted | growth mixture, random = wave, random_covariance = proportional |
| `growth_random_diagonal` | growth_mixture | fitted | growth mixture, random = wave, random_diagonal = TRUE |
| `growth_variance_equal` | growth_mixture | fitted | growth mixture, random intercept, residual variance = equal |
| `growth_membership` | growth_mixture | fitted | growth mixture, random intercept, membership = ~ motivation |
| `growth_weights` | growth_mixture | fitted | growth mixture, random intercept, per-student sampling weights |
| `growth_vcov_robust` | growth_mixture | fitted | growth mixture, random = wave, vcov_type = robust |
| `growth_vcov_opg` | growth_mixture | fitted | growth mixture, random = wave, vcov_type = opg |
| `growth_refusal_poisson_random` | growth_mixture | refused | random effects with a non-gaussian family are refused |
| `growth_refusal_observation_level` | growth_mixture | refused | random effects with class_level = observation are refused |
| `svc_enumerate_lpa_basic` | enumeration | fitted | enumerate_lpa(): 1-3 profiles x the four basic structures (EEI/VVI/EEE/VVV) |
| `svc_enumerate_lpa_all_structures` | enumeration | fitted | enumerate_lpa(model = "all"): all 14 structures at 2 profiles |
| `svc_enumerate_lca` | enumeration | fitted | enumerate_lca(): 1-3 classes on four binary ESM items |
| `svc_enumerate_classes_two_level` | enumeration | fitted | enumerate_classes(): two-level grid, 1-2 profiles x 1-2 group classes x EEI/VVI |
| `svc_candidate_fit` | enumeration | fitted | candidate_fit() picks the 2x2 VVI fit out of a two-level grid |
| `svc_candidate_fit_ambiguous` | enumeration | refused | candidate_fit() by counts alone over a multi-structure grid is refused |
| `svc_bootstrap_lrt_lpa` | bootstrap_lrt | fitted | bootstrap_lrt(): single-level LPA 2 vs 3 profiles, 5 replicates [bootstrap_lrt] |
| `svc_bootstrap_lrt_lca` | bootstrap_lrt | fitted | bootstrap_lrt(): LCA 1 vs 2 classes on binary items, 5 replicates [bootstrap_lrt] |
| `svc_bootstrap_lrt_fiml` | bootstrap_lrt | fitted | bootstrap_lrt() with missing = "fiml": LPA 1 vs 2 profiles [bootstrap_lrt] |
| `svc_bootstrap_lrt_profiles` | bootstrap_lrt | fitted | bootstrap_lrt(): two-level 1 vs 2 profiles (one group class) [bootstrap_lrt] |
| `svc_bootstrap_lrt_covariates` | bootstrap_lrt | fitted | bootstrap_lrt(): covariate model (profile covariate z), 2 vs 3 profiles [bootstrap_lrt] |
| `svc_bootstrap_lrt_mixed_families` | bootstrap_lrt | fitted | bootstrap_lrt() of a plain null against a covariate alternative: refused (currently by an unclassed error) [bootstrap_lrt] |
| `svc_bootstrap_lrt_fixed` | bootstrap_lrt | fitted | bootstrap_lrt(): fixed measurement, 1 vs 2 group classes (fit_staged alternative, fixed = "measurement" null) [bootstrap_lrt] |
| `svc_fit_staged_one_class` | fit_staged | refused | fit_staged(n_group_classes = 1) is refused (currently by an unclassed error) |
| `svc_inference_lpa` | inference | fitted | parameter_inference(): Wald observed/robust/OPG and bootstrap (5) on a single-level LPA [inference_opg inference_bootstrap] |
| `svc_inference_lpa_full` | inference | fitted | parameter_inference(): Wald and bootstrap on a VVV single-level LPA [inference_opg inference_bootstrap] |
| `svc_inference_two_level` | inference | fitted | parameter_inference(): Wald (all three), bootstrap and boundary = "fix" on a two-level fit [inference_opg inference_bootstrap inference_fix] |
| `svc_inference_boundary_lca` | inference | fitted | a covariate LCA with a response probability at its floor: Wald refused, boundary = "fix" holds it [inference_fix inference_opg] |
| `svc_inference_bootstrap_lca` | inference | fitted | parameter_inference(method = "bootstrap") on a categorical LCA (5 replicates) [inference_bootstrap] |
| `svc_inference_bootstrap_mixed_extra` | inference | fitted | bootstrap inference on a mixed ordinal + count (Poisson) + continuous LPA [inference_bootstrap] |
| `svc_inference_bootstrap_vei` | inference | fitted | bootstrap inference on a constrained VEI structure [inference_bootstrap inference_opg] |
| `svc_inference_bootstrap_staged` | inference | fitted | bootstrap inference on a fit_staged() fit (measurement held) [inference_bootstrap] |
| `svc_inference_adjust_level` | inference | fitted | parameter_inference(level = 0.9, adjust = "holm", step = 1e-5) |
| `svc_three_step_two_level` | three_step | fitted | three_step(): BCH, proportional, modal, pairwise contrasts and a group-level outcome [three_step_bch three_step_proportional three_step_modal three_step_pairs three_step_groups] |
| `svc_three_step_independent` | three_step | fitted | three_step(vcov_type = "independent", ci_level = 0.9, adjust = "holm") with contrasts |
| `svc_three_step_single_level` | three_step | fitted | three_step() on a single-level 3-profile LPA: BCH, modal, pairs [three_step_bch three_step_modal three_step_pairs] |
| `svc_three_step_indicator_refused` | three_step | refused | three_step() with an outcome that is one of the indicators is refused |
| `svc_r3step_two_level` | r3step | fitted | r3step(): individual-level and group-level covariates [r3step r3step_groups] |
| `svc_r3step_by_group_class` | r3step | fitted | r3step(by_group_class = TRUE) |
| `svc_r3step_observed` | r3step | fitted | r3step(vcov_type = "observed", adjust = "holm", ci_level = 0.9) |
| `svc_sensitivity_two_level` | sensitivity | fitted | sensitivity() over seeds 1:2 on a two-level fit [sensitivity] |
| `svc_sensitivity_options` | sensitivity | fitted | sensitivity(seeds = 3:4, n_starts = 1, max_iter = 50, tol = 1e-6) on an LPA |
| `svc_pool_covariates` | pool_imputations | fitted | pool_imputations(): three hand-made imputations of a covariate model (Rubin's rules) |
| `svc_pool_robust` | pool_imputations | fitted | pool_imputations(vcov_type = "robust", level = 0.9) of a plain two-level LPA |
| `svc_pool_mice` | pool_imputations | fitted | pool_imputations() of a seeded mice::mice() `mids` (pmm, m = 3); covariate with missing values |
| `svc_fit_staged` | fit_staged | fitted | fit_staged(): measurement at one group class, then two group classes [starting_values starting_values_measurement] |
| `svc_fit_staged_categorical` | fit_staged | fitted | fit_staged() with categorical indicators and full covariance continuous parts |
| `svc_starting_values_full` | starting_values | fitted | starting_values(): full covariance, covariance auto/drop/keep and measurement only |
| `svc_diagnostics_report` | diagnostics | fitted | diagnostics(), descriptives(), report() and predict() on a two-level fit [diagnostics descriptives report predict] |
| `svc_diagnostics_overall` | diagnostics | fitted | diagnostics(by = "overall") on a single-level LPA |
| `svc_predict_types` | predict | fitted | predict(type = "posterior") and predict(type = "density") on new rows |
| `svc_descriptives_data` | descriptives | fitted | descriptives() of a data frame: overall, with id (ICC), and by a factor |

## Known findings recorded in the baseline

These unclassed (non-`latents_*`) errors are current behaviour and are recorded as such (`unexpected-errors.csv`). Fixing one is an intended change: re-record the affected case.

- `get_results(x, "bch_weights")` on a noise fit (`noise = TRUE`) fails with "arguments imply differing number of rows"; because `"all"` and `summary()` build that table, both fail for every noise fit (`profiles_noise_single`, `profiles_noise_VVV`, `profiles_prior_noise_EII`). The fingerprint then records each table separately (`tables_each`), so the other tables are still covered.
- `get_tna()` / `get_group_tna()` have no method for a general-engine `lta()` fit (`multilpa_lta`): bare S3 dispatch error (`lta_general_transition_covariates`).
- `lta(time = NULL)` refuses through an unclassed `stopifnot()` (`lta_refusal_no_time`).
- A `group_covariate` that varies within a group is refused unclassed (`cov_refusal_group_covariate_varies`).
- `bootstrap_lrt()` of a plain null against a covariate alternative is refused unclassed (`svc_bootstrap_lrt_mixed_families`).
- `fit_staged(n_group_classes = 1)` is refused unclassed (`svc_fit_staged_one_class`).
- Bootstrap `parameter_inference()` is refused (classed, `latents_unsupported_inference`) for covariate and additive-family fits; the `inference_bootstrap` service passes no `data`, which transition fits require, so the `lta_*_bootstrap` cases call it with `data =`.
