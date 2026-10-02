#' Fit a two-level latent class model
#'
#' Two-level latent class analysis for categorical indicators: [multilpa()]
#' with discrete indicators in `vars` treated as categorical unless named
#' in `ordinal` or `count` through `...`. Each profile is described by
#' unrestricted response probabilities for categorical items and the
#' corresponding structured parameters for ordinal and count items,
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
#' @param n_group_classes Number of group-level latent classes. Defaults to 2,
#'   or to 1 when `id = NULL` (a single-level model has no second level).
#' @param ... Further arguments to [multilpa()], such as `n_starts`, `seed`,
#'   `missing`, `tol` or `profile_covariates`. `ordinal` and `count` can name
#'   structured discrete indicators; the
#'   remaining indicators are categorical. `categorical` cannot be supplied.
#' @return A fitted model of class `multilpa` (or `multilpa_covariates` when
#'   membership covariates are given), exactly as [multilpa()] returns it with
#'   the corresponding discrete indicator types. Categorical response
#'   probabilities are in
#'   `get_results(fit, "responses")`, one row per profile, item and category.
#' @section Conditions:
#'   `latents_bad_argument` when `categorical` is supplied, since its
#'   columns are derived from `vars`, `ordinal` and `count` here. Every condition of
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
  if (missing(id)) .multilpa_missing_id("multilca", "`lca()`")
  stopifnot("`vars` must be a character vector of column names" =
              is.character(vars) && length(vars) >= 1L && !anyNA(vars))
  if ("categorical" %in% names(list(...))) {
    stop(errorCondition(paste(
      "`categorical` cannot be supplied to multilca(): every indicator in",
      "`vars` is categorical. Use multilpa() for mixed indicators."),
      class = "latents_bad_argument", call = NULL))
  }
  # Forward `n_group_classes` only when it was given, so that `id = NULL`
  # resolves to one group class exactly as it does in multilpa().
  extra <- list(...)
  categorical <- setdiff(vars, c(extra$ordinal, extra$count))
  arguments <- c(list(data = data, vars = vars, id = id, n_profiles = n_profiles,
                       categorical = categorical), extra)
  if (!missing(n_group_classes)) arguments$n_group_classes <- n_group_classes
  fit <- do.call(multilpa, arguments)
  # Report the call the caller wrote, as multilpa() does for its own calls.
  fit$call <- match.call()
  fit
}
