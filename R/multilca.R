#' Fit a two-level latent class model
#'
#' Fits a two-level latent class model for discrete indicators and returns a
#' fitted `multilpa` object. It is [multilpa()] with every indicator in `vars`
#' treated as categorical, except those named in `ordinal` or `count` through
#' `...`. Each profile has unrestricted response probabilities for the
#' categorical indicators and structured parameters for the ordinal and count
#' ones; the group classes differ in the probability of each profile. For
#' models that combine categorical and continuous indicators, use [multilpa()]
#' and name the categorical ones in its `categorical` argument.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the indicator columns. Numeric, integer, logical,
#'   character and factor columns are accepted; factor levels keep their
#'   declared order and other types are sorted.
#' @param id Name of the group identifier column, or `NULL` for a single-level
#'   model.
#' @param n_profiles Number of observation-level latent classes (profiles).
#' @param n_group_classes Number of group-level latent classes. Defaults to 2,
#'   or to 1 when `id = NULL`.
#' @param ... Further arguments to [multilpa()], such as `n_starts`, `seed`,
#'   `missing`, `tol` or `profile_covariates`. `ordinal` and `count` name
#'   ordinal and count indicators; the remaining indicators are categorical.
#'   `categorical` is accepted but has no further effect; its names must be in
#'   `vars`.
#' @return A fitted `multilpa` object (`multilpa_covariates` when membership
#'   covariates are given), as [multilpa()] returns it. Categorical response
#'   probabilities are in `get_results(fit, "responses")`, one row per
#'   profile, item and category.
#' @section Conditions:
#'   `latents_bad_categorical` when `categorical` names a column not in
#'   `vars`; otherwise the conditions of [multilpa()].
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
  # Forward `n_group_classes` only when it was given, so that `id = NULL`
  # resolves to one group class exactly as it does in multilpa().
  extra <- list(...)
  # A `categorical` written as for multilpa() or enumerate_classes() adds to
  # the set rather than replacing it, so every indicator stays categorical.
  categorical <- .latents_lca_categorical(vars, extra)
  extra$categorical <- NULL
  arguments <- c(list(data = data, vars = vars, id = id, n_profiles = n_profiles,
                       categorical = categorical), extra)
  if (!missing(n_group_classes)) arguments$n_group_classes <- n_group_classes
  fit <- do.call(multilpa, arguments)
  # Report the call the caller wrote, as multilpa() does for its own calls.
  fit$call <- match.call()
  fit
}

#' The categorical indicators of an all-categorical verb
#'
#' Every indicator in `vars` not named in `ordinal` or `count`, together with
#' any `categorical` the caller wrote as for [multilpa()]. A name in
#' `categorical` that is not in `vars` is refused here, before an enumeration
#' could turn it into a grid of failed candidates.
#' @param vars The indicator names.
#' @param extra `list(...)` as the verb received it.
#' @return A character vector of indicator names.
#' @noRd
.latents_lca_categorical <- function(vars, extra) {
  stopifnot("`categorical` must be a character vector of column names" =
              is.null(extra$categorical) ||
              (is.character(extra$categorical) && !anyNA(extra$categorical)))
  unknown <- setdiff(extra$categorical, vars)
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "`categorical` names %s, which %s not in `vars`.",
      paste(sprintf("`%s`", unknown), collapse = ", "),
      if (length(unknown) == 1L) "is" else "are"),
      class = "latents_bad_categorical", call = NULL))
  }
  union(setdiff(vars, c(extra$ordinal, extra$count)), extra$categorical)
}
