#' Coefficients, covariance, intervals and likelihood of latents fits
#'
#' `coef()`, `vcov()`, `confint()`, `logLik()` and `nobs()` for multilevel
#' profile and class models ([multilpa()]), covariate models, latent transition
#' models and, for `coef()`, `vcov()` and `confint()`, pooled imputations.
#' Mixture regressions, growth
#' mixtures and the additive, cross-level and general transition families
#' document the same methods with their [get_results()] page.
#'
#' @param object A fitted model: from [multilpa()] (class `multilpa`), a
#'   covariate model (`multilpa_covariates`), a latent transition model
#'   (`multilpa_transitions`), or pooled imputations (`latents_pooled`:
#'   `coef()`, `vcov()` and `confint()`).
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
#'   and designs), to floating-point precision.
#' @param parm Optional coefficient names or indices; defaults to every
#'   estimated coefficient. A coefficient held by `fixed` has no interval and
#'   raises `latents_held_parameter`; an unknown name or an out-of-range index
#'   raises `latents_bad_argument`.
#' @param level Confidence level, strictly between 0 and 1 (default 0.95).
#' @param step Finite-difference step for the observed information (covariate
#'   models).
#' @param vcov_type `"observed"`, `"robust"` or `"opg"`, as for
#'   [parameter_inference()] (covariate models).
#' @param boundary `"error"` or `"fix"`, as for [parameter_inference()]:
#'   `"fix"` holds response probabilities on their bound, whose rows and
#'   columns are then zero (covariate models).
#' @param ... For `vcov()` and `confint()` of a `multilpa` fit, and `confint()`
#'   of a covariate model, passed to [parameter_inference()] (for example
#'   `method = "bootstrap"`, `vcov_type = "robust"` or `step`); an argument
#'   [parameter_inference()] does not take is refused with
#'   `latents_bad_argument`. `vcov()` of a covariate model warns
#'   `latents_unused_argument` for an argument it does not take, and refuses
#'   `method = "bootstrap"` with `latents_unsupported_inference`. Ignored
#'   otherwise.
#' @return `coef()`: a named numeric vector. `vcov()`: a named symmetric matrix
#'   over the same parameters. `confint()`: a matrix with one row per parameter
#'   and the lower and upper bounds as columns. `logLik()`: a `"logLik"` object
#'   with `df` and `nobs` attributes. `nobs()`: the number of independent groups
#'   counted by the fitted BIC. With `drop_empty = TRUE`, groups without any
#'   observed indicator are excluded; with `drop_empty = FALSE`, they are counted.
#'   For a single-level fit these groups are its individual observations.
#' @details Categorical parameter names escape reserved characters within
#'   indicator and category components (`%`, `:`, `,`, `[` and `]`) using percent
#'   codes. This keeps distinct parameters distinct when their original labels
#'   contain delimiters. Ordinary names are unchanged; use `names(coef(object))`
#'   or numeric positions when selecting parameters for `confint()`.
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
