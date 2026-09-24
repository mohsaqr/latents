#' Fit a two-level latent class model
#'
#' Two-level latent class analysis for categorical indicators: [multilpa()]
#' with every indicator in `vars` treated as categorical. Each profile is
#' described by an unrestricted probability for every category of every item,
#' and the group classes differ in how probable each profile is for their
#' observations. For models that combine categorical and continuous
#' indicators, call [multilpa()] and name the categorical ones in its
#' `categorical` argument.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the categorical indicator columns. Numeric, integer,
#'   logical, character and factor columns are accepted; factor levels keep
#'   their declared order and other types are sorted.
#' @param id Name of the group identifier column, or `NULL` for a single-level
#'   model.
#' @param n_profiles Number of observation-level latent classes (profiles).
#' @param n_group_classes Number of group-level latent classes.
#' @param ... Further arguments to [multilpa()], such as `n_starts`, `seed`,
#'   `missing`, `tol` or `profile_covariates`. `categorical` is set to `vars`
#'   and cannot be supplied.
#' @return A fitted model of class `multilpa` (or `multilpa_covariates` when
#'   membership covariates are given), exactly as [multilpa()] returns it with
#'   `categorical = vars`. The response probabilities are in
#'   `get_results(fit, "responses")`, one row per profile, item and category.
#' @section Conditions:
#'   `latents_bad_argument` when `categorical` is supplied, since every
#'   indicator is categorical by definition here. Every condition of
#'   [multilpa()] can also be raised.
#' @seealso [multilpa()] for continuous and mixed indicators, and
#'   `vignette("lca", package = "latents")`.
#' @examples
#' activities <- c("time_with_friends", "on_social_media", "tv_video_games")
#' fit <- multilca(subset(student_esm, day <= 1), vars = activities,
#'                 id = "student", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 1, seed = 1)
#' get_results(fit, "responses")
#' get_results(fit, "profile_probabilities")
#' @export
multilca <- function(data, vars, id, n_profiles, n_group_classes = 2L, ...) {
  stopifnot("`vars` must be a character vector of column names" =
              is.character(vars) && length(vars) >= 1L && !anyNA(vars))
  if ("categorical" %in% names(list(...))) {
    stop(errorCondition(paste(
      "`categorical` cannot be supplied to multilca(): every indicator in",
      "`vars` is categorical. Use multilpa() for mixed indicators."),
      class = "latents_bad_argument", call = NULL))
  }
  fit <- multilpa(data, vars = vars, id = id, n_profiles = n_profiles,
                  n_group_classes = n_group_classes, categorical = vars, ...)
  # Report the call the caller wrote, as multilpa() does for its own calls.
  fit$call <- match.call()
  fit
}
