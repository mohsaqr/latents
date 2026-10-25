# Conditions raised by the package

Lists the classes of the errors, warnings and messages the package
raises, and what raises each. A condition can be caught by its class
with [`tryCatch()`](https://rdrr.io/r/base/conditions.html) or
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html). The
classes are stable; the message text may change between versions.

## Value

Nothing; this page documents conditions rather than a function.

## Input and data contracts

- `latents_bad_data`:

  The data frame, the indicator columns or the group column are invalid:
  too few rows, a duplicated or absent column name, an indicator that is
  not a numeric vector, an indicator that is constant or entirely
  unobserved, missing or non-finite values under `missing = "error"`, or
  a missing group identifier.

- `latents_bad_argument`:

  An argument has an invalid value: a scalar such as `n_profiles`,
  `n_group_classes`, `n_starts`, `max_iter`, `tol` or `min_variance`; a
  `seed` or `seeds` that is not a whole number accepted by
  [`set.seed()`](https://rdrr.io/r/base/Random.html); a `what` naming a
  table or view the object does not have; or an argument the verb or
  table does not use.

- `latents_bad_time`:

  `time` is missing, or repeats a position within a group.

- `latents_bad_categorical`:

  A categorical indicator cannot be encoded, or `min_probability` is too
  large for its number of categories.

- `latents_bad_start`:

  Supplied starting values lack a block the model needs, have the wrong
  shape, or carry labels naming an indicator or category the fit does
  not have.

- `latents_bad_fixed`:

  A `fixed` request names a block the model does not have or that
  `start` does not supply, or `start` is absent.

- `latents_bad_stage`:

  A `measurement` fit given to
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  already has more than one group class; was fitted to different
  indicators, profiles or measurement options than the stage requests;
  or was fitted with `weights`, `centering`, or ordinal or count
  indicators.

- `latents_bad_transition`:

  A transition model was requested with fewer than two profiles, or with
  no group observed at two occasions.

- `latents_bad_inference_data`:

  The data passed to a verb with a fit do not reproduce the fit: a
  different number of rows, different values or a different group
  ordering.

- `latents_no_indicator_data`:

  A fit carries no stored indicators, so the request needs the original
  data.

- `latents_incomplete_fit`:

  A fit lacks fields the request needs, usually because it was made by
  an older version.

## Estimation and identification

- `latents_unidentified`:

  The requested model is not identified from the data: more profiles
  than distinct observed indicator rows, more group classes than groups,
  or more than one group class with only one profile or only singleton
  groups.

- `latents_all_starts_failed`:

  Every EM start failed, so no fit was produced. The message gives the
  distinct start errors.

- `latents_empty_profile`:

  A profile has no observed responses for an indicator, so its response
  probabilities are not identified; or a group class has no effective
  membership, so its profile prevalences are not defined.

- `latents_too_few_groups`:

  There are too few independent groups to form the covariance of the
  reported quantities. Raised by robust and OPG
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  and by
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  and
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
  whose cluster-robust contributions sum to zero at the estimate and so
  span at most one dimension fewer than the number of groups.

- `latents_inseparable_classes`:

  The classes are not separated well enough for the requested three-step
  correction.

- `latents_no_converge`:

  A fit required by the request did not converge.

- `latents_boundary_fit`:

  A variance or response probability is at its bound, a mixing
  probability is zero, or a membership logit coefficient diverges
  (separation). The parameter is on the edge of its space, where the
  normal approximation of Wald inference does not hold, so no standard
  error is reported.

- `latents_singular_information`:

  The observed information is not positive definite or is numerically
  singular, so no standard error can be formed.

- `latents_unknown_candidate`:

  The requested profile and group-class combination was not enumerated.

- `latents_failed_candidate`:

  An enumerated candidate could not be fitted. The message gives that
  candidate's error.

- `latents_no_group_classes`:

  A table, level or view that needs discrete group classes was requested
  from a model without them, for example `level = "groups"` for a fit
  without group classes, or a group-level plot of a fit with fewer than
  two group classes.

## Model comparison

- `latents_bad_nesting`:

  Two models are not nested as the comparison requires, including a
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  pair whose held measurement blocks differ, are held at different
  values, or are held in only one of the two models.

- `latents_incomparable_models`:

  Two models were fitted to different data, or differ in covariance
  structure, centering, scaling mode or `drop_empty` where the
  comparison requires them to match.

- `latents_parallel_failure`:

  The worker cluster could not be started, or a broken cluster's
  surviving processes could not be stopped safely. A process lost during
  fitting is otherwise retried and reported with `latents_worker_lost`
  instead.

- `latents_worker_lost`:

  One or more worker tasks returned no result, possibly because their
  processes were killed for memory. Missing results are retried once
  with fewer workers. The warning reports the numbers lost, recovered
  and still failed, and suggests fewer `cores`. Persistent losses become
  failed starts or candidate rows, retaining other results. A socket
  connection failure makes its current bounded chunk uncertain; earlier
  chunks are retained, and the uncertain chunk is retried.

- `latents_suspect_categorical`:

  A numeric indicator fitted as categorical has non-integer values or
  more than ten categories, so it may be continuous. Every distinct
  value is fitted as a category. The warning names the indicators (field
  `indicators`) and is given once per call, including an enumeration.

- `latents_progress`:

  A message from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  with `verbose = TRUE` (the default when `cores > 1`): one as the
  candidates are prepared, one as their starts are fitted, and one when
  the grid is done. Muffle it with
  `suppressMessages(classes = "latents_progress")` or `verbose = FALSE`.

- `latents_candidate_failed`:

  A candidate of an enumeration that was fitted as one worker task lost
  that worker twice. It is recorded as the candidate's error, in its row
  of the grid, and the other candidates are returned.

- `latents_reversed_likelihood`:

  The model with more parameters has the lower likelihood, so at least
  one fit is at a local optimum.

- `latents_bad_outcome`:

  A distal outcome given to
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  is of a type the requested method cannot handle, varies within a group
  when a group-level outcome is requested, or is a measurement
  indicator.

- `latents_bad_covariate`:

  A predictor given to
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  varies within a group when group-level membership is predicted, makes
  the design rank deficient, or is a measurement indicator; or a
  predictor given to
  `multilpa(profile_covariates =, group_covariates =)` is an indicator
  or the group identifier, or makes the membership design rank
  deficient.

- `multilpa_indicator_reused`:

  A measurement indicator was given as an outcome to
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  or as a predictor to
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md).
  Raised together with `latents_bad_outcome` or `latents_bad_covariate`.

- `latents_bad_scores`:

  Score contributions for the robust sandwich could not be formed.

## Unsupported requests

These are raised when the package does not implement the requested
analysis.

- `latents_unsupported_inference`:

  Standard errors are not available for this fit. Raised by
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  for a fit made with `prior` (Wald path), for the cross-level families,
  for a covariate fit made by an older version, and for a fit recording
  a covariance structure code this version does not know; and by
  `parameter_inference(method = "bootstrap")` for an additive or
  covariate fit, a multilevel growth mixture, or a fit holding a
  measurement block through `fixed`. Also raised by
  `vcov(scale = "unconstrained")` for a bootstrap covariance.

- `latents_unsupported_prior`:

  `multilpa(prior = )` was combined with an option the conjugate prior
  is not defined for: the VEE, EVE, VVE or EVV covariance structure
  (mclust defines none), categorical indicators, a `fixed` block,
  membership covariates, or missing indicator values.

- `latents_unsupported_noise`:

  `multilpa(noise = TRUE)` was combined with more than one group class,
  categorical indicators, a `fixed` block, membership covariates or
  missing indicator values; or a fit with a noise component was given to
  a verb that does not support one:
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md),
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md),
  bivariate residuals,
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  or a posterior plot.

- `latents_unsupported_structure`:

  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  was asked for a covariance structure other than EEI, VVI, EEE or VVV
  together with `profile_covariates` or `group_covariates`.

- `latents_bad_weights`:

  `weights` does not name a column, or its values are missing, negative,
  non-finite, vary within a unit, or are positive for fewer than two
  units.

- `latents_unsupported_weights`:

  A weighted fit was given to a verb or option that does not support
  sampling weights:
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  `enumerate_regressions(bootstrap = )`,
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md),
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
  `prior` or `noise`; or a `vcov_type` other than `"robust"` was
  requested for it.

- `latents_unsupported_indicator`:

  `ordinal` or `count` indicators were combined with an option that does
  not support them: membership covariates, `start`, `fixed`, `prior` or
  `noise`.

- `latents_pooling_failed`:

  A fit or its inference failed for one imputation in
  [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md),
  or the imputations' fits report different parameters. No imputation is
  dropped; the message names the imputation and gives the original
  reason.

- `latents_bootstrap_failed`:

  Fewer than two resamples produced a usable fit. Raised by
  `parameter_inference(method = "bootstrap")`; the message gives the
  first reason a resample failed.

- `latents_no_free_parameters`:

  The fit holds every parameter fixed, so there is no standard error to
  report. Raised by
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  and [`vcov()`](https://rdrr.io/r/stats/vcov.html).

- `latents_held_parameter`:

  A parameter named in `parm` was held fixed by the fit and has no
  confidence interval. Raised by
  [`confint()`](https://rdrr.io/r/stats/confint.html).

- `latents_unsupported_sensitivity`:

  [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  was called on a model family whose refit needs arguments the shared
  refit does not carry, whose profile labels cannot yet be aligned
  between two fits, or on a directly fixed fit whose original free
  starts cannot be replayed.

- `latents_unsupported_three_step`:

  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  or [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  was called on a membership-covariate fit, whose classification errors
  depend on the predictors; the correction uses one unconditional error
  matrix.

- `latents_unsupported_bootstrap`:

  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  was called on a person-centred fit, for which raw data cannot be
  simulated under the null model because the group baselines were
  removed before estimation.

- `latents_unsupported_tna`:

  [`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md)
  or
  [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  was called on a general transition fit (an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit with
  occasion-varying or covariate-dependent transitions, a second order,
  movers and stayers, or ordinal, count or weighted data), which has no
  single transition matrix; its transitions are given by
  `get_results(x, "transitions")`.

- `latents_no_plot`:

  `plot(x, what = "all")` was called on an object whose plot method has
  no views to draw.

- `latents_no_time`:

  A sequence table or verb was requested for a fit made without `time`.

- `latents_unsupported_prediction`:

  [`predict()`](https://rdrr.io/r/stats/predict.html) was called on a
  growth mixture, a membership-covariate fit or a latent transition fit,
  for which prediction on new rows is not implemented; results for the
  fitted rows are available from
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).

- `latents_missing_package`:

  A required package from `Suggests` is not installed. Raised by
  [`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md)
  and
  [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  when `tna` is absent; the message names the package and how to install
  it.

## Plotting

- `latents_no_continuous`, `latents_no_categorical`:

  The requested panel needs indicators of a kind the fit does not have.

- `latents_nothing_to_describe`:

  [`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
  has no variables to describe: `vars` named none, or the frame has no
  numeric column outside `id` and `by`.

- `latents_nothing_to_plot`:

  The requested panel has no rows to draw.

- `latents_inaccurate_probability`:

  The numerical integral behind a Lo-Mendell-Rubin p-value did not
  settle, so the p-value is not reported.

- `latents_no_local_dependence`:

  `get_results(fit, "local_dependence")` on a fit made without
  `local_dependence`.

- `latents_no_direct_effects`:

  `get_results(fit, "direct_effects")` on a fit made without
  `direct_effects`.

- `latents_partial_pair`:

  A warning: rows observe one item of a `local_dependence` pair and not
  the other, so the pair is treated as missing in those rows.

- `latents_no_lmr`:

  `get_results(x, "lmr")` on a grid with no two converged candidates one
  profile or one group class apart. Also of class
  `latents_nothing_to_plot`.

- `latents_invalid_lmr`:

  The alternative of a nested LMR comparison has a lower likelihood than
  the null beyond numerical roundoff. The comparison table retains its
  row with unavailable statistics and the reason in `note`; improve the
  alternative's fit before comparing.

- `latents_no_tree`:

  A tree would have fewer than two rows, or the group level of a grid
  with one group class was asked for. Also of class
  `latents_nothing_to_plot`; `get_results(x, "all")` leaves the `"tree"`
  table out.

- `latents_unsupported_view`:

  A grid verb
  ([`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_selection()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  or [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a grid)
  was asked for a view that this kind of grid cannot draw, such as the
  tree of a group-class family or regression grid. The message gives the
  reason and the views the grid offers; `plot_views(x)` lists them. Also
  of class `latents_bad_argument`.

- `latents_unknown_category`, `latents_unknown_criterion`:

  A named category or information criterion does not exist in the fit.
  An unknown criterion of a grid is also of class
  `latents_bad_argument`.

- `latents_bad_scale`:

  An unrecognised plotting scale was requested, or
  `scale = "standardized"` met an indicator with zero or undefined
  standard deviation.

## Warnings

These are raised when a result is returned but qualified.

- `latents_failed_starts`:

  Some, but not all, EM starts failed. The remaining starts are in
  `get_results(fit, "starts")`.

- `latents_not_identified`:

  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  with a binary outcome and one trial per class assignment, where a
  mixture of Bernoulli regressions is not identified.

- `latents_missing_data`:

  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  with `missing = "error"` met a missing value in a variable the model
  uses.

- `latents_rows_dropped`:

  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  with `missing = "omit"` dropped incomplete rows; the message states
  how many.

- `latents_no_valid_start`:

  Every
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  start degenerated: a class emptied or its variance collapsed.

- `latents_degenerate_start`:

  Some
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  starts degenerated and were set aside.

- `latents_separation`:

  A
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  logit or log-link coefficient diverged (perfect separation); it and
  its standard error are unreliable.

- `latents_em_decrease`:

  The
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  EM likelihood decreased beyond rounding error. This indicates a defect
  and should be reported.

- `latents_unconverged`:

  The best start had not converged when `max_iter` was reached, so the
  returned estimate is not a maximum.

- `latents_boundary`:

  A variance, or a covariance eigenvalue, is at `min_variance`, so the
  fit is bound-active and its standard errors are not valid at that
  parameter.

- `latents_small_classes`:

  A profile or group class has an effective membership below one
  observation; or, in
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  with `method = "ml"`, a class is carried by a single unit, so its
  standard error is reported as `NA`.

- `latents_thin_profiles`:

  A Gaussian component with free variances of its own has posterior
  membership mass below the number of free parameters in its Gaussian
  density. This is a heuristic prompt to inspect the solution, not a
  proof of singularity or nonidentification. Covariances held at zero
  and held means are excluded, and held variances do not trigger it.
  Memberships use normalized sampling weights and only rows with an
  observed Gaussian indicator; this mass is not an independent sample
  count or a weight effective sample size. Transition models check
  occasion-specific measurement separately; cross-level models check
  both Gaussian levels. The estimate is returned unchanged. See
  `get_results(fit, "gaussian_diagnostics")` for continuous covariance
  and support diagnostics; none supplies an automatic spike
  classification or an identification guarantee.

- `latents_all_starts_boundary`:

  `boundary = "exclude"` found every start with a variance at
  `min_variance`, so none was set aside and the best start across all of
  them is reported.

- `latents_weak_class`:

  `multilpa(family = "additive")` only: a group class is supported by
  fewer than 50 effective groups (see `effective_groups` in
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)),
  because it is small, poorly separated, or both. Standard errors are
  reported, but intervals for class means and weights may be
  miscalibrated.

- `latents_failed_replicates`:

  Some bootstrap replicates failed validation, so the bootstrap p-value
  is `NA`.

- `latents_single_level`:

  A message, not a warning: `id = NULL` was passed, so the fit has one
  observation per unit and no second level. Raised by
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  on every such fit and once by an
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  grid; [`suppressMessages()`](https://rdrr.io/r/base/message.html) or a
  `latents_single_level` handler silences it.

- `latents_bootstrap_dropped`:

  Some resamples did not produce a usable fit and were excluded from the
  interval. Raised by `parameter_inference(method = "bootstrap")`; the
  message states how many.

- `latents_sensitivity_dropped`:

  One or more seeds did not produce a fit in
  [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md).
  Their rows are `NA`, and the warning states how many failed and the
  first reason.

- `latents_extreme_coefficients`:

  A membership logit coefficient, measured per root-mean-square unit of
  its covariate, is large enough that the class is close to separated,
  so the estimate is driven by scaling, a sparse class or separation
  rather than by the data.

- `latents_empty_transition_row`:

  A profile is never occupied before a final occasion, so its transition
  row is uniform by construction rather than estimated;
  `get_results(fit, "transitions", estimated = FALSE)` lists such rows.
  The transition-network handoff also warns when it includes such a row.

- `latents_unverified_alignment`:

  A supplied `data` frame has too little fitted information to verify
  row order (for example only a repeated group identifier). The result
  assumes the rows are in fitting order.

- `latents_random_boundary`:

  A growth mixture's random-effect covariance is degenerate in a class
  (a variance at its floor, or perfectly correlated random effects), or,
  with `residual_correlation`, the residual variance, process variance
  or process rate is on its bound; the standard errors of those
  parameters are not reported.

- `latents_singular_covariance`:

  A growth mixture with correlated residuals met a marginal covariance
  that is numerically singular (a process variance far larger than the
  residual variance). During fitting this is treated as a rejected step;
  it reaches the caller only for parameter values supplied by hand.

- `latents_unused_argument`:

  A method that accepts `...` (several
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods,
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
  of a data frame) was given an argument it does not use. The call
  proceeds and the warning names the argument.

- `latents_unused_category`:

  An ordinal indicator is a factor with levels never observed in the
  data. Their intercepts diverge, so Wald inference is refused; call
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html) before
  fitting. It is a warning rather than an error because bootstrap refits
  keep the fitted categories.

- `latents_no_standard_errors`:

  A message, not a warning: a table that includes standard errors by
  default could not compute them for this fit (a bound-active or
  unconverged fit, singular information, or a family or noise component
  without inference), so they are `NA` and the message gives the reason.
  Passing `data` to `get_results(x, "profiles")` or `"responses"` raises
  the underlying condition instead.

## Examples

``` r
set.seed(1)
example_data <- data.frame(
  school = rep(seq_len(8), each = 6),
  score_a = stats::rnorm(48), score_b = stats::rnorm(48)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#> Warning: Profile 1: 2.0 Gaussian membership mass for 4 free means and variances. This is a conservative diagnostic, not an identification test; inspect the Gaussian covariance diagnostics.
# A sequence verb on a fit made without `time` is catchable by class.
tryCatch(get_results(fit, "sequences"),
         latents_no_time = function(condition) {
  "this fit carries no ordering"
})
#> [1] "this fit carries no ordering"
```
