# Single-level latent profile and latent class analysis.
#
# The package's models are two-level: observations nested in groups, with
# latent classes at both levels. A single-level model is the special case with
# one observation per unit and no group classes. `multilpa(id = NULL)` fits it;
# these two verbs name it, so that the ordinary analysis reads as what it is.

#' Fit a latent profile model
#'
#' @description
#' Latent profile analysis describes heterogeneity in continuous indicators
#' through a specified number of latent profiles. Each profile represents a
#' pattern of indicator means and within-profile variation. `lpa()` estimates
#' these patterns using a finite mixture of multivariate normal distributions.
#'
#' By default, indicators are independent within profiles and their variances
#' differ across profiles. Alternative covariance structures allow associations
#' between indicators within profiles and equality constraints across profiles.
#' The fitted model provides profile means, variances, proportions, and
#' posterior membership probabilities for each observation, together with
#' model-fit and classification summaries.
#'
#' The model treats observations as independent. For observations nested within
#' groups, use [multilpa()]. Additional measurement options support mixed
#' continuous and discrete indicators.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the indicator columns. Numeric columns are continuous
#'   indicators; categorical ones are named in `categorical`.
#' @param n_profiles Number of profiles; a positive integer.
#' @param ... Further arguments for [multilpa()]: the covariance structure
#'   (`model`, or `variance_model` and `covariance_model`, or `volume`,
#'   `shape` and `orientation`), `missing`, `categorical`, `n_starts`,
#'   `seed`, `boundary`, `prior`, `noise` and `profile_covariates`.
#'   Arguments `id` and `n_group_classes` are not accepted.
#' @return A fitted object as returned by [multilpa()], usually of class
#'   `multilpa`, or `multilpa_covariates` with membership covariates. Its
#'   `call` records the original `lpa()` call.
#' @section Conditions:
#'   `latents_bad_argument` when `id` or `n_group_classes` is passed; the
#'   conditions [multilpa()] raises otherwise.
#' @references
#' Oberski, D. (2016). Mixture models: Latent profile and latent class
#' analysis. In J. Robertson & M. Kaptein (Eds.), *Modern Statistical Methods
#' for HCI* (pp. 275--287). Springer.
#' @seealso [lca()] for categorical indicators; [multilpa()] for the two-level
#'   model; [enumerate_lpa()] to compare numbers of profiles.
#' @examples
#' fit <- lpa(srl, c("cognitive_strategies", "intrinsic_value",
#'                   "self_efficacy", "self_regulation", "test_anxiety"),
#'            n_profiles = 2, variance_model = "equal",
#'            covariance_model = "full", n_starts = 3, seed = 1)
#' fit
#' plot(fit, what = "raincloud")
#' @export
lpa <- function(data, vars, n_profiles, ...) {
  .latents_check_data_frame(data, substitute(data))
  .multilpa_single_level_call("lpa", "multilpa", list(...))
  fit <- .multilpa_quiet_single_level(
    multilpa(data, vars, id = NULL, n_profiles = n_profiles, ...))
  # The call the caller wrote, so that `update()` refits from the caller's
  # objects rather than from this function's argument names.
  fit$call <- match.call()
  fit
}

#' Fit a latent class model
#'
#' @description
#' Latent class analysis describes heterogeneity in categorical response
#' patterns through a specified number of latent classes. Each class is
#' characterized by the probabilities of the response categories for each
#' indicator. These probabilities describe how response patterns differ
#' between classes.
#'
#' `lca()` estimates the class-specific response probabilities and class
#' proportions, together with posterior membership probabilities for each
#' observation. Indicators are assumed independent conditional on class
#' membership unless associations between selected indicators are explicitly
#' specified. The posterior probabilities describe uncertainty in
#' classification; an observation's assigned class is a summary of these
#' probabilities.
#'
#' The model treats observations as independent. For observations nested within
#' groups, use [multilca()]. Ordinal and count indicators can also be fitted
#' using the structured measurement models described in [multilpa()].
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the discrete indicator columns. Categories are taken
#'   from the data. Columns named in `ordinal` or `count` (through `...`) are
#'   modelled as such; the rest are categorical.
#' @param n_classes Number of classes; a positive integer.
#' @param ... Further arguments for [multilca()], such as `missing`,
#'   `n_starts`, `seed`, `min_probability` and `profile_covariates`, except
#'   `id` and `n_group_classes`.
#' @return A fitted object as returned by [multilca()], usually of class
#'   `multilpa`, or `multilpa_covariates` with membership covariates. Its
#'   `call` records the original `lca()` call. Classes are labelled
#'   `profile_1`, `profile_2`, and so on in the returned tables.
#' @section Conditions:
#'   `latents_bad_argument` when `id`, `n_group_classes` or `n_profiles` is
#'   passed; the conditions [multilca()] raises otherwise.
#' @references
#' Lazarsfeld, P. F., & Henry, N. W. (1968). *Latent Structure Analysis*.
#' Houghton Mifflin.
#' @seealso [lpa()] for continuous indicators; [multilca()] for the
#'   two-level model; [enumerate_lca()] to compare numbers of classes.
#' @examples
#' fit <- lca(student_esm, c("time_with_friends", "on_social_media",
#'                           "tv_video_games", "sports"),
#'            n_classes = 2, n_starts = 3, seed = 1)
#' fit
#' plot(fit, what = "heatmap")
#' @export
lca <- function(data, vars, n_classes, ...) {
  .latents_check_data_frame(data, substitute(data))
  .multilpa_single_level_call("lca", "multilca", list(...))
  if ("n_profiles" %in% names(list(...))) {
    stop(errorCondition(
      "`lca()` takes the number of classes as `n_classes`, not `n_profiles`.",
      class = "latents_bad_argument", call = NULL))
  }
  fit <- .multilpa_quiet_single_level(
    multilca(data, vars, id = NULL, n_profiles = n_classes, ...))
  fit$call <- match.call()
  fit
}

#' Compare latent profile models
#'
#' @description
#' `enumerate_lpa()` fits latent profile models with different numbers of
#' profiles and specified covariance structures. It brings together evidence
#' about fit, classification, and estimation stability to support selection
#' of the number of profiles and their covariance specification. Each fitted
#' candidate remains available for detailed examination.
#'
#' @details
#' This is the single-level interface to [enumerate_classes()].
#' For observations nested in groups, use [enumerate_classes()] with the
#' grouping column as `id`, which also varies the number of group classes.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the indicator columns.
#' @param n_profiles Numbers of profiles to compare, specified as positive
#'   integers.
#' @param model The covariance structures: `"basic"`, the default, for the
#'   four that combine equal or varying variances with covariances absent or
#'   present (`EEI`, `VVI`, `EEE`, `VVV`); `"all"` for all 14; or their
#'   three-letter codes. See [enumerate_classes()].
#' @param ... Further arguments for [multilpa()], such as `missing`,
#'   `n_starts`, `seed`, `boundary` and `categorical`, except `id` and
#'   `n_group_classes`.
#' @return An object of class `multilpa_enumeration`, as [enumerate_classes()]
#'   returns: [as.data.frame()] gives one row per candidate, [summary()] the
#'   candidate minimizing each criterion, [plot()] the criteria against the
#'   number of profiles, and [get_candidate()] one fitted candidate.
#' @section Conditions:
#'   `latents_bad_argument` when `id` or `n_group_classes` is passed, or when
#'   `model` names an unknown structure.
#' @references
#' Celeux, G., & Govaert, G. (1995). Gaussian parsimonious clustering models.
#' *Pattern Recognition*, 28(5), 781--793.
#' @seealso [lpa()] to fit one model; [enumerate_lca()] for categorical
#'   indicators; [bootstrap_lrt()] to test one number of profiles against
#'   another.
#' @examples
#' models <- enumerate_lpa(srl, c("cognitive_strategies", "intrinsic_value",
#'                                "self_efficacy", "self_regulation",
#'                                "test_anxiety"),
#'                         n_profiles = 1:3, n_starts = 2, seed = 1)
#' summary(models)
#' plot(models)
#' @export
enumerate_lpa <- function(data, vars, n_profiles = 1:4, model = "basic", ...) {
  .multilpa_single_level_call("enumerate_lpa", "enumerate_classes", list(...))
  # Preserve omission: the default gives way to legacy covariance switches
  # and to all-categorical data, whereas an explicit model does not.
  data_expression <- substitute(data)
  .latents_check_data_frame(data, data_expression)
  if (missing(model)) {
    return(.latents_candidate_data(.multilpa_quiet_single_level(
      enumerate_classes(data, vars, id = NULL, n_profiles = n_profiles, ...)),
      data_expression))
  }
  .latents_candidate_data(.multilpa_quiet_single_level(
    enumerate_classes(data, vars, id = NULL, n_profiles = n_profiles,
                      model = model, ...)), data_expression)
}

#' Compare latent class models
#'
#' @description
#' `enumerate_lca()` fits latent class models with different numbers of classes
#' to examine heterogeneity in categorical response patterns. It brings
#' together evidence about fit, classification, and estimation stability to
#' support selection of the number of classes. Each fitted candidate remains
#' available for examination of its class-specific response probabilities and
#' class proportions.
#'
#' @details
#' This is the single-level interface to [enumerate_classes()] for categorical
#' indicators. Every indicator is categorical unless named in `ordinal` or
#' `count`, as in [lca()].
#'
#' For observations nested in groups, use [enumerate_classes()] with the
#' grouping column as `id` and the indicators as `categorical`.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the categorical indicator columns.
#' @param n_classes Numbers of classes to compare, specified as positive
#'   integers.
#' @param ... Further arguments for [multilca()], such as `missing`,
#'   `n_starts`, `seed` and `min_probability`, except `id` and
#'   `n_group_classes`. Every indicator is categorical unless named in
#'   `ordinal` or `count`. `categorical` is accepted but has no further
#'   effect; its names must be in `vars`.
#' @return An object of class `multilpa_enumeration`, as [enumerate_classes()]
#'   returns: [as.data.frame()] gives one row per candidate, [summary()] the
#'   candidate minimizing each criterion, [plot()] the criteria against the
#'   number of classes, and [get_candidate()] one fitted candidate. The number
#'   of classes is reported in the `n_profiles` column of the enumeration
#'   table and selected through the `n_profiles` argument of [get_candidate()].
#' @section Conditions:
#'   `latents_bad_argument` when `id` or `n_group_classes` is passed, or when
#'   `model` is passed (no indicator is continuous); `latents_bad_categorical`
#'   when `categorical` names a column not in `vars`.
#' @seealso [lca()] to fit one model; [enumerate_lpa()] for continuous
#'   indicators.
#' @examples
#' models <- enumerate_lca(student_esm, c("time_with_friends",
#'                                        "on_social_media", "sports"),
#'                         n_classes = 1:3, n_starts = 2, seed = 1)
#' summary(models)
#' @export
enumerate_lca <- function(data, vars, n_classes = 1:4, ...) {
  .latents_check_data_frame(data, substitute(data))
  extra <- list(...)
  .multilpa_single_level_call("enumerate_lca", "enumerate_classes", extra)
  # As in multilca(), a `categorical` adds to the set; every indicator stays
  # categorical.
  categorical <- .latents_lca_categorical(vars, extra)
  extra$categorical <- NULL
  .latents_candidate_data(.multilpa_quiet_single_level(
    do.call(enumerate_classes, c(list(data = data, vars = vars, id = NULL,
                                      n_profiles = n_classes,
                                      categorical = categorical), extra))),
    substitute(data))
}

#' Refuse the two-level arguments in a single-level verb
#' @param verb,two_level The verb's name and the two-level verb to point to.
#' @param extra `list(...)` as the verb received it.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_single_level_call <- function(verb, two_level, extra) {
  refused <- intersect(c("id", "n_group_classes"), names(extra))
  if (length(refused) > 0L) {
    stop(errorCondition(sprintf(paste(
      "`%s()` fits a single-level model, which has no %s. For observations",
      "nested in groups use `%s()` with the grouping column as `id`."),
      verb, paste(sprintf("`%s`", refused), collapse = " or "), two_level),
      class = "latents_bad_argument", call = NULL))
  }
  invisible(NULL)
}

#' Evaluate a single-level fit without its single-level notice
#'
#' The notice tells a caller of `multilpa(id = NULL)` that the model has one
#' level. A caller of `lpa()` or `lca()` asked for exactly that by name.
#' @param expr The fitting call.
#' @return Its value.
#' @noRd
.multilpa_quiet_single_level <- function(expr) {
  withCallingHandlers(expr, latents_single_level = function(notice) {
    invokeRestart("muffleMessage")
  })
}
