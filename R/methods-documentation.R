#' Print latents results
#'
#' Print methods for the objects the package returns: fitted models,
#' enumerations, bootstrap tests, covariate and transition models, pooled
#' imputations, their summaries, result tables and plot lists. Each prints a
#' compact overview; [get_results()] returns the tables themselves as data
#' frames.
#'
#' @param x An object returned by a latents function, or the `summary()` of
#'   one.
#' @param digits Number of significant digits printed.
#' @param rows How many rows of each table to print. A longer table is shown to
#'   that depth, with its remaining row count and the `get_results()` call that
#'   returns it whole.
#' @param n The most rows of a result table printed.
#' @param ... Passed to the underlying data frame printing, or ignored.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @seealso [get_results()] for the tables, [latents-summary] and
#'   [latents-as-data-frame].
#' @examples
#' fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
#'                 "student", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 2, seed = 1)
#' print(fit, rows = 3)
#' @name latents-print
NULL

#' Summarise latents results
#'
#' Summary methods for fitted models, enumerations, bootstrap tests, covariate
#' and transition models, mixture regressions and pooled imputations. A summary
#' collects the model's main tables and diagnostics; printing it shows each
#' table to a set depth.
#'
#' @param object An object returned by a latents function.
#' @param level Confidence level for the reported intervals, strictly between
#'   zero and one.
#' @param vcov_type `"observed"`, `"robust"` or `"opg"`, as for
#'   [parameter_inference()]. For mixture regressions and growth mixtures,
#'   `NULL` reuses the inference stored with the fit, or, when there is none,
#'   uses robust errors for a weighted fit and observed otherwise.
#' @param ... Ignored; present for compatibility with [summary()].
#' @return An object of class `summary_<class>` holding the summarised tables.
#'   Print it, or read any table with [get_results()] or [as.data.frame()].
#' @seealso [latents-print], [latents-as-data-frame], [get_results()].
#' @examples
#' fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
#'                 "student", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 2, seed = 1)
#' print(summary(fit), rows = 3)
#' @name latents-summary
NULL

#' Convert latents results to a data frame
#'
#' `as.data.frame()` returns the primary table of a latents result: the
#' measurement parameters of a fitted model, the candidates of an enumeration,
#' the test of a bootstrap, and so on. It is the same table as `get_results(x)`;
#' [get_results()] with `what =` returns the others.
#'
#' @param x An object returned by a latents function, or the `summary()` of
#'   one.
#' @param row.names Passed to [data.frame()]; `NULL` gives default row names.
#' @param optional Ignored; present for compatibility with [as.data.frame()].
#' @param ... Must be empty. An argument here raises `latents_bad_argument`
#'   naming it, rather than being dropped.
#' @return A base `data.frame`, one row per parameter, candidate or test
#'   depending on the object.
#' @seealso [get_results()], [latents-print], [latents-summary].
#' @examples
#' fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
#'                 "student", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 2, seed = 1)
#' head(as.data.frame(fit))
#' @name latents-as-data-frame
NULL

#' Coefficients, covariance, intervals and likelihood of latents fits
#'
#' `coef()`, `vcov()`, `confint()`, `logLik()` and `nobs()` for multilevel
#' profile and class models ([multilpa()]), covariate models, latent transition
#' models and, for `vcov()`, pooled imputations. Mixture regressions, growth
#' mixtures and the additive, cross-level and general transition families
#' document the same methods with their [get_results()] page.
#'
#' @param object A fitted model: from [multilpa()] (class `multilpa`), a
#'   covariate model (`multilpa_covariates`), a latent transition model
#'   (`multilpa_transitions`), or pooled imputations (`latents_pooled`,
#'   `vcov()` only).
#' @param scale `"natural"`, the default, or `"unconstrained"`. Natural
#'   coefficients are in the units `coef()` reports: means and variances in the
#'   indicators' units, probabilities as probabilities, with the covariance
#'   carried from the estimation scale by the delta method. Unconstrained is the
#'   scale the model is estimated on: log variances (diagonal), log-Cholesky
#'   coordinates (full covariance), log volumes, shapes and orientations (a
#'   structure that constrains them across profiles), and baseline-category
#'   logits for mixing and response probabilities. In a covariate model, means
#'   and membership coefficients are the same on both scales.
#' @param data Optional. The data frame the model was fitted to. When omitted
#'   it is rebuilt from what the fit stores (indicators, identifiers, occasions
#'   and designs), which round-trips exactly; supplying it is the stronger check
#'   that the caller still holds that frame.
#' @param parm Optional coefficient names or indices; defaults to every
#'   estimated coefficient. Naming a coefficient that `fixed` held raises
#'   `latents_held_parameter`, because a held value has no interval.
#' @param level Confidence level strictly between zero and one.
#' @param step Finite-difference step for the observed information (covariate
#'   models).
#' @param vcov_type `"observed"`, `"robust"` or `"opg"`, as for
#'   [parameter_inference()] (covariate models).
#' @param boundary `"error"` or `"fix"`, as for [parameter_inference()]:
#'   `"fix"` holds response probabilities on their bound, whose rows and
#'   columns are then zero (covariate models).
#' @param ... For `vcov()` and `confint()` of a `multilpa` fit, and `confint()`
#'   of a covariate model, passed to [parameter_inference()] (for example
#'   `method = "bootstrap"`, `vcov_type = "robust"` or `step`). Ignored
#'   otherwise.
#' @return `coef()`: a named numeric vector. `vcov()`: a named symmetric matrix
#'   over the same parameters. `confint()`: a matrix with one row per parameter
#'   and the lower and upper bounds as columns. `logLik()`: a `"logLik"` object
#'   with `df` and `nobs` attributes. `nobs()`: the number of groups, the
#'   independent units of the likelihood.
#' @seealso [parameter_inference()] for the tidy table of estimates, standard
#'   errors and intervals.
#' @examples
#' fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
#'                 "student", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 2, seed = 1)
#' coef(fit)
#' logLik(fit)
#' nobs(fit)
#' @name latents-model-methods
NULL
