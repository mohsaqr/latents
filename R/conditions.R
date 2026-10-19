#' Conditions raised by the package
#'
#' Lists the classes of the errors, warnings and messages the package raises,
#' and what raises each. A condition can be caught by its class with
#' [tryCatch()] or [withCallingHandlers()]. The classes are stable; the
#' message text may change between versions.
#'
#' @section Input and data contracts:
#' \describe{
#'   \item{`latents_bad_data`}{The data frame, the indicator columns or the
#'     group column are invalid: too few rows, a duplicated or absent column
#'     name, an indicator that is not a numeric vector, an indicator that is
#'     constant or entirely unobserved, missing or non-finite values under
#'     `missing = "error"`, or a missing group identifier.}
#'   \item{`latents_bad_argument`}{An argument has an invalid value: a scalar
#'     such as `n_profiles`, `n_group_classes`, `n_starts`, `max_iter`, `tol`
#'     or `min_variance`; a `seed` or `seeds` that is not a whole number
#'     accepted by `set.seed()`; a `what` naming a table or view the object
#'     does not have; or an argument the verb or table does not use.}
#'   \item{`latents_bad_time`}{`time` is missing, or repeats a position within
#'     a group.}
#'   \item{`latents_bad_categorical`}{A categorical indicator cannot be
#'     encoded, or `min_probability` is too large for its number of
#'     categories.}
#'   \item{`latents_bad_start`}{Supplied starting values lack a block the
#'     model needs, have the wrong shape, or carry labels naming an indicator
#'     or category the fit does not have.}
#'   \item{`latents_bad_fixed`}{A `fixed` request names a block the model does
#'     not have or that `start` does not supply, or `start` is absent.}
#'   \item{`latents_bad_stage`}{A `measurement` fit given to [fit_staged()]
#'     already has more than one group class; was fitted to different
#'     indicators, profiles or measurement options than the stage requests; or
#'     was fitted with `weights`, `centering`, or ordinal or count indicators.}
#'   \item{`latents_bad_transition`}{A transition model was requested with
#'     fewer than two profiles, or with no group observed at two occasions.}
#'   \item{`latents_bad_inference_data`}{The data passed to a verb with a fit
#'     do not reproduce the fit: a different number of rows, different values
#'     or a different group ordering.}
#'   \item{`latents_no_indicator_data`}{A fit carries no stored indicators, so
#'     the request needs the original data.}
#'   \item{`latents_incomplete_fit`}{A fit lacks fields the request needs,
#'     usually because it was made by an older version.}
#' }
#'
#' @section Estimation and identification:
#' \describe{
#'   \item{`latents_unidentified`}{The requested model is not identified from
#'     the data: more profiles than distinct observed indicator rows, more
#'     group classes than groups, or more than one group class with only one
#'     profile or only singleton groups.}
#'   \item{`latents_all_starts_failed`}{Every EM start failed, so no fit was
#'     produced. The message gives the distinct start errors.}
#'   \item{`latents_empty_profile`}{A profile has no observed responses for an
#'     indicator, so its response probabilities are not identified; or a group
#'     class has no effective membership, so its profile prevalences are not
#'     defined.}
#'   \item{`latents_too_few_groups`}{There are too few independent groups to
#'     form the covariance of the reported quantities. Raised by robust and
#'     OPG [parameter_inference()], and by [three_step()] and [r3step()], whose
#'     cluster-robust contributions sum to zero at the estimate and so span at
#'     most one dimension fewer than the number of groups.}
#'   \item{`latents_inseparable_classes`}{The classes are not separated well
#'     enough for the requested three-step correction.}
#'   \item{`latents_no_converge`}{A fit required by the request did not
#'     converge.}
#'   \item{`latents_boundary_fit`}{A variance or response probability is at
#'     its bound, a mixing probability is zero, or a membership logit
#'     coefficient diverges (separation). The parameter is on the edge of its
#'     space, where the normal approximation of Wald inference does not hold,
#'     so no standard error is reported.}
#'   \item{`latents_singular_information`}{The observed information is not
#'     positive definite or is numerically singular, so no standard error can
#'     be formed.}
#'   \item{`latents_unknown_candidate`}{The requested profile and group-class
#'     combination was not enumerated.}
#'   \item{`latents_failed_candidate`}{An enumerated candidate could not be
#'     fitted. The message gives that candidate's error.}
#'   \item{`latents_no_group_classes`}{A table, level or view that needs
#'     discrete group classes was requested from a model without them, for
#'     example `level = "groups"` for a fit without group classes, or a
#'     group-level plot of a fit with fewer than two group classes.}
#' }
#'
#' @section Model comparison:
#' \describe{
#'   \item{`latents_bad_nesting`}{Two models are not nested as the comparison
#'     requires, including a [bootstrap_lrt()] pair whose held measurement
#'     blocks differ, are held at different values, or are held in only one of
#'     the two models.}
#'   \item{`latents_incomparable_models`}{Two models were fitted to different
#'     data, or differ in covariance structure, centering, scaling mode or
#'     `drop_empty` where the comparison requires them to match.}
#'   \item{`latents_parallel_failure`}{With `cores` above one, a worker
#'     process failed as a whole (for example by running out of memory), so
#'     some candidates of an enumeration were lost. An estimation error of a
#'     single candidate is recorded in the grid instead.}
#'   \item{`latents_reversed_likelihood`}{The model with more parameters has
#'     the lower likelihood, so at least one fit is at a local optimum.}
#'   \item{`latents_bad_outcome`}{A distal outcome given to [three_step()] is
#'     of a type the requested method cannot handle, varies within a group
#'     when a group-level outcome is requested, or is a measurement
#'     indicator.}
#'   \item{`latents_bad_covariate`}{A predictor given to [r3step()] varies
#'     within a group when group-level membership is predicted, makes the
#'     design rank deficient, or is a measurement indicator; or a predictor
#'     given to `multilpa(profile_covariates =, group_covariates =)` is an
#'     indicator or the group identifier, or makes the membership design rank
#'     deficient.}
#'   \item{`multilpa_indicator_reused`}{A measurement indicator was given as
#'     an outcome to [three_step()] or as a predictor to [r3step()]. Raised
#'     together with `latents_bad_outcome` or `latents_bad_covariate`.}
#'   \item{`latents_bad_scores`}{Score contributions for the robust sandwich
#'     could not be formed.}
#' }
#'
#' @section Unsupported requests:
#' These are raised when the package does not implement the requested
#' analysis.
#' \describe{
#'   \item{`latents_unsupported_inference`}{Standard errors are not available
#'     for this fit. Raised by [parameter_inference()] for a fit made with
#'     `prior` (Wald path), for the cross-level families, for a covariate fit
#'     made by an older version, and for a fit recording a covariance
#'     structure code this version does not know; and by
#'     `parameter_inference(method = "bootstrap")` for an additive or
#'     covariate fit, a multilevel growth mixture, or a fit holding a
#'     measurement block through `fixed`. Also raised by
#'     `vcov(scale = "unconstrained")` for a bootstrap covariance.}
#'   \item{`latents_unsupported_prior`}{`multilpa(prior = )` was combined
#'     with an option the conjugate prior is not defined for: the VEE, EVE,
#'     VVE or EVV covariance structure (mclust defines none), categorical
#'     indicators, a `fixed` block, membership covariates, or missing
#'     indicator values.}
#'   \item{`latents_unsupported_noise`}{`multilpa(noise = TRUE)` was combined
#'     with more than one group class, categorical indicators, a `fixed`
#'     block, membership covariates or missing indicator values; or a fit with
#'     a noise component was given to a verb that does not support one:
#'     [parameter_inference()], [three_step()], [r3step()], [bootstrap_lrt()],
#'     [starting_values()], bivariate residuals, [fit_staged()] or a
#'     posterior plot.}
#'   \item{`latents_unsupported_structure`}{`multilpa()` was asked for a
#'     covariance structure other than EEI, VVI, EEE or VVV together with
#'     `profile_covariates` or `group_covariates`.}
#'   \item{`latents_bad_weights`}{`weights` does not name a column, or its
#'     values are missing, negative, non-finite, vary within a unit, or are
#'     positive for fewer than two units.}
#'   \item{`latents_unsupported_weights`}{A weighted fit was given to a verb
#'     or option that does not support sampling weights: [bootstrap_lrt()],
#'     `enumerate_regressions(bootstrap = )`, [three_step()], [r3step()],
#'     `prior` or `noise`; or a `vcov_type` other than `"robust"` was
#'     requested for it.}
#'   \item{`latents_unsupported_indicator`}{`ordinal` or `count` indicators
#'     were combined with an option that does not support them: membership
#'     covariates, `start`, `fixed`, `prior` or `noise`.}
#'   \item{`latents_pooling_failed`}{A fit or its inference failed for one
#'     imputation in [pool_imputations()], or the imputations' fits report
#'     different parameters. No imputation is dropped; the message names the
#'     imputation and gives the original reason.}
#'   \item{`latents_bootstrap_failed`}{Fewer than two resamples produced a
#'     usable fit. Raised by `parameter_inference(method = "bootstrap")`; the
#'     message gives the first reason a resample failed.}
#'   \item{`latents_no_free_parameters`}{The fit holds every parameter fixed,
#'     so there is no standard error to report. Raised by
#'     [parameter_inference()] and [vcov()].}
#'   \item{`latents_held_parameter`}{A parameter named in `parm` was held
#'     fixed by the fit and has no confidence interval. Raised by
#'     [confint()].}
#'   \item{`latents_unsupported_sensitivity`}{[sensitivity()] was called on a
#'     model family whose refit needs arguments the shared refit does not
#'     carry, whose profile labels cannot yet be aligned between two fits, or
#'     on a directly fixed fit whose original free starts cannot be
#'     replayed.}
#'   \item{`latents_unsupported_three_step`}{[three_step()] or [r3step()] was
#'     called on a membership-covariate fit, whose classification errors
#'     depend on the predictors; the correction uses one unconditional error
#'     matrix.}
#'   \item{`latents_unsupported_bootstrap`}{[bootstrap_lrt()] was called on a
#'     person-centred fit, for which raw data cannot be simulated under the
#'     null model because the group baselines were removed before
#'     estimation.}
#'   \item{`latents_unsupported_tna`}{[get_tna()] or [get_group_tna()] was
#'     called on a general transition fit (an [lta()] fit with
#'     occasion-varying or covariate-dependent transitions, a second order,
#'     movers and stayers, or ordinal, count or weighted data), which has no
#'     single transition matrix; its transitions are given by
#'     `get_results(x, "transitions")`.}
#'   \item{`latents_no_plot`}{`plot(x, what = "all")` was called on an object
#'     whose plot method has no views to draw.}
#'   \item{`latents_no_time`}{A sequence table or verb was requested for a
#'     fit made without `time`.}
#'   \item{`latents_unsupported_prediction`}{`predict()` was called on a
#'     growth mixture, a membership-covariate fit or a latent transition fit,
#'     for which prediction on new rows is not implemented; results for the
#'     fitted rows are available from [get_results()].}
#'   \item{`latents_missing_package`}{A required package from `Suggests` is
#'     not installed. Raised by [get_tna()] and [get_group_tna()] when `tna`
#'     is absent; the message names the package and how to install it.}
#' }
#'
#' @section Plotting:
#' \describe{
#'   \item{`latents_no_continuous`, `latents_no_categorical`}{The requested
#'     panel needs indicators of a kind the fit does not have.}
#'   \item{`latents_nothing_to_describe`}{[descriptives()] has no variables to
#'     describe: `vars` named none, or the frame has no numeric column outside
#'     `id` and `by`.}
#'   \item{`latents_nothing_to_plot`}{The requested panel has no rows to draw.}
#'   \item{`latents_inaccurate_probability`}{The numerical integral behind a
#'     Lo-Mendell-Rubin p-value did not settle, so the p-value is not reported.}
#'   \item{`latents_no_local_dependence`}{`get_results(fit, "local_dependence")`
#'     on a fit made without `local_dependence`.}
#'   \item{`latents_no_direct_effects`}{`get_results(fit, "direct_effects")`
#'     on a fit made without `direct_effects`.}
#'   \item{`latents_partial_pair`}{A warning: rows observe one item of a
#'     `local_dependence` pair and not the other, so the pair is treated as
#'     missing in those rows.}
#'   \item{`latents_no_lmr`}{`get_results(x, "lmr")` on a grid with no two
#'     converged candidates one profile or one group class apart. Also of class
#'     `latents_nothing_to_plot`.}
#'   \item{`latents_no_tree`}{A tree would have fewer than two rows, or the
#'     group level of a grid with one group class was asked for. Also of
#'     class `latents_nothing_to_plot`; `get_results(x, "all")` leaves the
#'     `"tree"` table out.}
#'   \item{`latents_unsupported_view`}{A grid verb ([plot_enumeration()],
#'     [plot_selection()], [plot_tree()] or `plot()` of a grid) was asked for
#'     a view that this kind of grid cannot draw, such as the tree of a
#'     group-class family or regression grid. The message gives the reason
#'     and the views the grid offers; `plot_views(x)` lists them. Also of
#'     class `latents_bad_argument`.}
#'   \item{`latents_unknown_category`, `latents_unknown_criterion`}{A named
#'     category or information criterion does not exist in the fit. An
#'     unknown criterion of a grid is also of class `latents_bad_argument`.}
#'   \item{`latents_bad_scale`}{An unrecognised plotting scale was requested,
#'     or `scale = "standardized"` met an indicator with zero or undefined
#'     standard deviation.}
#' }
#'
#' @section Warnings:
#' These are raised when a result is returned but qualified.
#' \describe{
#'   \item{`latents_failed_starts`}{Some, but not all, EM starts failed. The
#'     remaining starts are in `get_results(fit, "starts")`.}
#'   \item{`latents_not_identified`}{[mixture_regression()] with a binary
#'     outcome and one trial per class assignment, where a mixture of
#'     Bernoulli regressions is not identified.}
#'   \item{`latents_missing_data`}{[mixture_regression()] with
#'     `missing = "error"` met a missing value in a variable the model uses.}
#'   \item{`latents_rows_dropped`}{[mixture_regression()] with
#'     `missing = "omit"` dropped incomplete rows; the message states how
#'     many.}
#'   \item{`latents_no_valid_start`}{Every [mixture_regression()] start
#'     degenerated: a class emptied or its variance collapsed.}
#'   \item{`latents_degenerate_start`}{Some [mixture_regression()] starts
#'     degenerated and were set aside.}
#'   \item{`latents_separation`}{A [mixture_regression()] logit or log-link
#'     coefficient diverged (perfect separation); it and its standard error
#'     are unreliable.}
#'   \item{`latents_em_decrease`}{The [mixture_regression()] EM likelihood
#'     decreased beyond rounding error. This indicates a defect and should be
#'     reported.}
#'   \item{`latents_unconverged`}{The best start had not converged when
#'     `max_iter` was reached, so the returned estimate is not a maximum.}
#'   \item{`latents_boundary`}{A variance, or a covariance eigenvalue, is at
#'     `min_variance`, so the fit is bound-active and its standard errors are
#'     not valid at that parameter.}
#'   \item{`latents_small_classes`}{A profile or group class has an effective
#'     membership below one observation.}
#'   \item{`latents_all_starts_boundary`}{`boundary = "exclude"` found every
#'     start with a variance at `min_variance`, so none was set aside and the
#'     best start across all of them is reported.}
#'   \item{`latents_weak_class`}{`multilpa(family = "additive")` only: a
#'     group class is supported by fewer than 50 effective groups (see
#'     `effective_groups` in [get_results()]), because it is small, poorly
#'     separated, or both. Standard errors are reported, but intervals for
#'     class means and weights may be miscalibrated.}
#'   \item{`latents_failed_replicates`}{Some bootstrap replicates failed
#'     validation, so the bootstrap p-value is `NA`.}
#'   \item{`latents_single_level`}{A message, not a warning: `id = NULL` was
#'     passed, so the fit has one observation per unit and no second level.
#'     Raised by [multilpa()] on every such fit and once by an
#'     [enumerate_classes()] grid; `suppressMessages()` or a
#'     `latents_single_level` handler silences it.}
#'   \item{`latents_bootstrap_dropped`}{Some resamples did not produce a
#'     usable fit and were excluded from the interval. Raised by
#'     `parameter_inference(method = "bootstrap")`; the message states how
#'     many.}
#'   \item{`latents_sensitivity_dropped`}{One or more seeds did not produce a
#'     fit in [sensitivity()]. Their rows are `NA`, and the warning states how
#'     many failed and the first reason.}
#'   \item{`latents_extreme_coefficients`}{A membership logit coefficient,
#'     measured per root-mean-square unit of its covariate, is large enough
#'     that the class is close to separated, so the estimate is driven by
#'     scaling, a sparse class or separation rather than by the data.}
#'   \item{`latents_empty_transition_row`}{A profile is never occupied before
#'     a final occasion, so its transition row is uniform by construction
#'     rather than estimated; `get_results(fit, "transitions", estimated =
#'     FALSE)` lists such rows. The transition-network handoff also warns when
#'     it includes such a row.}
#'   \item{`latents_unverified_alignment`}{A supplied `data` frame has too
#'     little fitted information to verify row order (for example only a
#'     repeated group identifier). The result assumes the rows are in fitting
#'     order.}
#'   \item{`latents_random_boundary`}{A growth mixture's random-effect
#'     covariance is degenerate in a class (a variance at its floor, or
#'     perfectly correlated random effects), or, with `residual_correlation`,
#'     the residual variance, process variance or process rate is on its
#'     bound; the standard errors of those parameters are not reported.}
#'   \item{`latents_singular_covariance`}{A growth mixture with correlated
#'     residuals met a marginal covariance that is numerically singular (a
#'     process variance far larger than the residual variance). During
#'     fitting this is treated as a rejected step; it reaches the caller only
#'     for parameter values supplied by hand.}
#'   \item{`latents_unused_argument`}{A method that accepts `...` (several
#'     `plot()` methods, `vcov()`, `predict()` and [descriptives()] of a data
#'     frame) was given an argument it does not use. The call proceeds and
#'     the warning names the argument.}
#'   \item{`latents_unused_category`}{An ordinal indicator is a factor with
#'     levels never observed in the data. Their intercepts diverge, so Wald
#'     inference is refused; call `droplevels()` before fitting. It is a
#'     warning rather than an error because bootstrap refits keep the fitted
#'     categories.}
#'   \item{`latents_no_standard_errors`}{A message, not a warning: a table
#'     that includes standard errors by default could not compute them for
#'     this fit (a bound-active or unconverged fit, singular information, or a
#'     family or noise component without inference), so they are `NA` and the
#'     message gives the reason. Passing `data` to
#'     `get_results(x, "profiles")` or `"responses"` raises the underlying
#'     condition instead.}
#' }
#'
#' @return Nothing; this page documents conditions rather than a function.
#' @examples
#' set.seed(1)
#' example_data <- data.frame(
#'   school = rep(seq_len(8), each = 6),
#'   score_a = stats::rnorm(48), score_b = stats::rnorm(48)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' # A sequence verb on a fit made without `time` is catchable by class.
#' tryCatch(get_results(fit, "sequences"),
#'          latents_no_time = function(condition) {
#'   "this fit carries no ordering"
#' })
#' @name latents-conditions
NULL
