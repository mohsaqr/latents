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
#'   and membership coefficients are the same on both scales. With
#'   `method = "bootstrap"` only the natural scale exists, and
#'   `scale = "unconstrained"` raises `latents_unsupported_inference`.
#' @param data Optional. The data frame the model was fitted to. When omitted
#'   it is rebuilt from what the fit stores (indicators, identifiers, occasions
#'   and designs), which round-trips to floating-point precision; supplying it is the stronger check
#'   that the caller still holds that frame.
#' @param parm Optional coefficient names or indices; defaults to every
#'   estimated coefficient. Naming a coefficient that `fixed` held raises
#'   `latents_held_parameter`, because a held value has no interval; an
#'   unknown name or an out-of-range index raises `latents_bad_argument`.
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
