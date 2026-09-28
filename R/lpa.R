# Single-level latent profile and latent class analysis.
#
# The package's models are two-level: observations nested in groups, with
# latent classes at both levels. A single-level model is the special case with
# one observation per unit and no group classes. `multilpa(id = NULL)` fits it;
# these two verbs name it, so that the ordinary analysis reads as what it is.

#' Latent profile analysis
#'
#' Fits a single-level latent profile model: a finite mixture of multivariate
#' normal distributions, in which each observation belongs to one of
#' `n_profiles` unobserved profiles and each profile has its own mean on every
#' indicator. It is the model [multilpa()] fits with `id = NULL`, and it
#' returns the same object, so every accessor, plot and inference verb of the
#' package applies to it.
#'
#' For observations nested in groups (students in schools, occasions in
#' persons), use [multilpa()] with the grouping column as `id`, which adds
#' latent classes of groups above the profiles.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the indicator columns. Numeric columns are continuous
#'   indicators; name categorical ones in `categorical`.
#' @param n_profiles Number of profiles, a positive whole number.
#' @param ... Further arguments for [multilpa()]: the covariance structure
#'   (`variance_model`, `covariance_model`, or `volume`, `shape`,
#'   `orientation`), `missing`, `categorical`, `n_starts`, `seed`, `prior`,
#'   `noise`, `profile_covariates` and the others it documents, except `id`
#'   and `n_group_classes`, which a single-level model does not have.
#' @return A fitted `multilpa` object, as [multilpa()] returns.
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
  .multilpa_single_level_call("lpa", "multilpa", list(...))
  .multilpa_quiet_single_level(
    multilpa(data, vars, id = NULL, n_profiles = n_profiles, ...))
}

#' Latent class analysis
#'
#' Fits a single-level latent class model: each observation belongs to one of
#' `n_classes` unobserved classes, and within a class the categorical
#' indicators are independent, each with its own response probabilities. It
#' is the model [multilca()] fits with `id = NULL`, and it returns the same
#' object, so every accessor, plot and inference verb of the package applies
#' to it. In the results the classes are labelled `profile_1`, `profile_2`,
#' and so on, as throughout the package.
#'
#' For observations nested in groups, use [multilca()] with the grouping
#' column as `id`, which adds latent classes of groups above the classes of
#' observations.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the categorical indicator columns. Their categories are
#'   taken from the data; no recoding is needed.
#' @param n_classes Number of classes, a positive whole number.
#' @param ... Further arguments for [multilca()], such as `missing`,
#'   `n_starts`, `seed`, `min_probability` and `profile_covariates`, except
#'   `id` and `n_group_classes`, which a single-level model does not have.
#' @return A fitted `multilpa` object, as [multilca()] returns.
#' @section Conditions:
#'   `latents_bad_argument` when `id` or `n_group_classes` is passed; the
#'   conditions [multilca()] raises otherwise.
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
  .multilpa_single_level_call("lca", "multilca", list(...))
  .multilpa_quiet_single_level(
    multilca(data, vars, id = NULL, n_profiles = n_classes, ...))
}

#' Compare latent profile models
#'
#' Fits single-level latent profile models with each number of profiles in
#' `n_profiles` under each covariance structure in `model`, and returns them
#' in one table for comparison by information criteria, entropy and the
#' diagnostics of each fit. It is [enumerate_classes()] with `id = NULL`, and
#' it returns the same object.
#'
#' For observations nested in groups, use [enumerate_classes()] with the
#' grouping column as `id`, which also crosses the number of group classes.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the indicator columns.
#' @param n_profiles Positive whole numbers of profiles to compare.
#' @param model The covariance structures: `"basic"`, the default, for the
#'   four that combine equal or varying variances with covariances absent or
#'   present (`EEI`, `VVI`, `EEE`, `VVV`); `"all"` for all 14; or their
#'   three-letter codes. See [enumerate_classes()].
#' @param ... Further arguments for [multilpa()], such as `missing`,
#'   `n_starts`, `seed` and `categorical`, except `id` and `n_group_classes`.
#' @return An object of class `multilpa_enumeration`, as [enumerate_classes()]
#'   returns: [as.data.frame()] gives one row per candidate, [summary()] the
#'   candidate each criterion prefers, [plot()] the criteria against the
#'   number of profiles, and [candidate_fit()] one fitted candidate.
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
  .multilpa_quiet_single_level(
    enumerate_classes(data, vars, id = NULL, n_profiles = n_profiles,
                      model = model, ...))
}

#' Compare latent class models
#'
#' Fits single-level latent class models with each number of classes in
#' `n_classes` and returns them in one table for comparison by information
#' criteria, entropy and the diagnostics of each fit. Every indicator is
#' categorical, as in [lca()]. It is [enumerate_classes()] with `id = NULL`,
#' and it returns the same object; the number of classes is its `n_profiles`
#' column.
#'
#' For observations nested in groups, use [enumerate_classes()] with the
#' grouping column as `id` and the indicators as `categorical`.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the categorical indicator columns.
#' @param n_classes Positive whole numbers of classes to compare.
#' @param ... Further arguments for [multilca()], such as `missing`,
#'   `n_starts`, `seed` and `min_probability`, except `id`,
#'   `n_group_classes` and `categorical`.
#' @return An object of class `multilpa_enumeration`, as [enumerate_classes()]
#'   returns: [as.data.frame()] gives one row per candidate, [summary()] the
#'   candidate each criterion prefers, [plot()] the criteria against the
#'   number of classes, and [candidate_fit()] one fitted candidate, named by
#'   its `n_profiles`.
#' @section Conditions:
#'   `latents_bad_argument` when `id`, `n_group_classes`, `categorical` or
#'   `model` is passed.
#' @seealso [lca()] to fit one model; [enumerate_lpa()] for continuous
#'   indicators.
#' @examples
#' models <- enumerate_lca(student_esm, c("time_with_friends",
#'                                        "on_social_media", "sports"),
#'                         n_classes = 1:3, n_starts = 2, seed = 1)
#' summary(models)
#' @export
enumerate_lca <- function(data, vars, n_classes = 1:4, ...) {
  extra <- list(...)
  .multilpa_single_level_call("enumerate_lca", "enumerate_classes", extra)
  if ("categorical" %in% names(extra)) {
    stop(errorCondition(paste(
      "`categorical` cannot be supplied to enumerate_lca(): every indicator in",
      "`vars` is categorical. Use enumerate_lpa() for mixed indicators."),
      class = "latents_bad_argument", call = NULL))
  }
  .multilpa_quiet_single_level(
    enumerate_classes(data, vars, id = NULL, n_profiles = n_classes,
                      categorical = vars, ...))
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
