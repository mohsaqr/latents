#' Fit a two-level latent class model: classes of categorical responses and
#' classes of the groups they are nested in
#'
#' @description
#' `multilca()` fits a latent class model for categorical data with observations
#' nested within groups. It describes response patterns through latent classes
#' and differences between groups through the relative prevalence of those
#' classes. Groups with similar class compositions are represented by the same
#' latent group class.
#'
#' By default, the response probabilities defining the observation-level
#' classes are shared across groups. The model estimates these probabilities
#' together with class proportions and posterior membership probabilities,
#' allowing both observations and groups to be classified while retaining
#' information about classification uncertainty.
#'
#' @details
#' The function calls [multilpa()], treating
#' all indicators in `vars` as categorical except those specified through
#' `ordinal` or `count` in `...`. Categorical indicators have unrestricted
#' class-specific response probabilities; ordinal and count indicators use the
#' corresponding structured measurement models. For a combination of continuous
#' and categorical indicators, use [multilpa()] with `categorical`.
#'
#' @param data A data frame with one row per observation.
#' @param vars Names of the indicator columns. Numeric, integer, logical,
#'   character and factor columns are accepted, and every distinct value of a
#'   column is one category; factor levels retain their declared order and
#'   other types are sorted.
#' @param id Name of the group identifier column, or `NULL` for a single-level
#'   model ([lca()] fits that directly).
#' @param n_profiles Number of observation-level latent classes (profiles).
#' @param n_group_classes Number of group-level latent classes. Defaults to 2,
#'   or to 1 when `id = NULL`. With 1, the groups do not differ in their
#'   class proportions.
#' @param ... Further arguments, passed to [multilpa()], which documents all
#'   of them. The ones used most often:
#'   * `missing`: `"error"` (default) stops when an indicator has a missing
#'     value; `"fiml"` keeps every row and uses the responses it has.
#'   * `ordinal`, `count`: indicators in `vars` modelled as ordered
#'     categories or as counts; the remaining indicators are categorical.
#'   * `n_starts` (default 10), `seed`, `max_iter` and `tol`: the number of
#'     starts, the seed that makes them reproducible, and the limits of each
#'     start.
#'   * `min_probability`: the lower bound on a response probability (default
#'     `1e-10`).
#'   * `profile_covariates`, `group_covariates`: covariates that predict
#'     class and group-class membership.
#'   * `time`: the occasion column, used to order each group's rows in the
#'     sequence tables and [plot_sequences()].
#'
#'   `categorical` is accepted but has no further effect; its names must be in
#'   `vars`.
#' @return A fitted model of class `multilpa`, or `multilpa_covariates` when
#'   membership covariates are specified, as returned by [multilpa()]. The
#'   classes are numbered 1, 2, and so on in the `profile` column of its
#'   tables, and the group classes in `group_class`. [get_results()] returns
#'   the tables by name, among them:
#'   * `"responses"` (also `as.data.frame(fit)`): one row per class, item
#'     and category, with the response probability;
#'   * `"profile_probabilities"`: the proportion of each class within each
#'     group class, and the proportion of each group class;
#'   * `"assignments"`: the data with each row's most likely class, its
#'     group's most likely group class and the posterior probabilities;
#'   * `"group_posteriors"`: the posterior probability of each group class
#'     for each group;
#'   * `"classification"` and `"entropy"`: how clearly rows and groups are
#'     classified;
#'   * `"information_criteria"` and `"starts"`.
#' @section Conditions:
#'   Errors: `latents_bad_argument` when `id` is left out or an argument has
#'   a value that is not accepted; `latents_bad_categorical` when
#'   `categorical` names a column not in `vars`; `latents_bad_data` when a
#'   column is missing from `data` or an indicator has missing values under
#'   `missing = "error"`; `latents_unidentified` when more than one group
#'   class is requested with one class or with one row per group. Warnings:
#'   `latents_unconverged` when the reported start did not converge and
#'   `latents_small_classes` when a class or group class is nearly empty.
#'   Otherwise the conditions of [multilpa()]; [latents-conditions] lists
#'   them all.
#' @seealso [multilpa()] for continuous and mixed indicators, and
#'   `vignette("lca", package = "latents")`.
#' @examples
#' activities <- c("time_with_friends", "on_social_media", "tv_video_games")
#' fit <- multilca(subset(student_esm, day <= 1), vars = activities,
#'                 id = "student", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 2, seed = 1)
#' get_results(fit, "responses")
#' get_results(fit, "profile_probabilities")
#' @export
multilca <- function(data, vars, id, n_profiles, n_group_classes = 2L, ...) {
  if (missing(id)) .multilpa_missing_id("multilca", "`lca()`")
  .latents_check_data_frame(data, substitute(data))
  .latents_require("`vars` must be a character vector of column names" =
                     is.character(vars) && length(vars) >= 1L && !anyNA(vars))
  # Forward `n_group_classes` only when it was given, so that `id = NULL`
  # resolves to one group class exactly as it does in multilpa().
  extra <- list(...)
  # A `categorical` written as for multilpa() or enumerate_classes() adds to
  # the set rather than replacing it, so every indicator stays categorical.
  categorical <- .latents_lca_categorical(vars, extra)
  extra$categorical <- NULL
  # Use the promises already forced above and the evaluated dots. Replaying
  # the caller's expressions would run data generators and controls twice.
  # Bind these values to symbols so no forwarded condition carries the data.
  arguments <- c(list(data = data, vars = vars, id = id, n_profiles = n_profiles,
                       categorical = categorical), extra)
  if (!missing(n_group_classes)) arguments$n_group_classes <- n_group_classes
  fitting_env <- new.env(parent = parent.frame())
  fitting_env$multilpa <- multilpa
  fitting_call <- .latents_bind_values(as.call(c(list(quote(multilpa)), arguments)),
                                       fitting_env, parent.frame(), values = arguments)
  fit <- eval(fitting_call, fitting_env)
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
  # Exact names: `$` on a list matches partially, so `extra$count` read a
  # `count_model` or `count_dispersion` value as indicator names whenever
  # `count` itself was absent (audit 2026-11-01).
  categorical <- extra[["categorical"]]
  .latents_require("`categorical` must be a character vector of column names" =
                     is.null(categorical) ||
                     (is.character(categorical) && !anyNA(categorical)))
  unknown <- setdiff(categorical, vars)
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "`categorical` names %s, which %s not in `vars`.",
      paste(sprintf("`%s`", unknown), collapse = ", "),
      if (length(unknown) == 1L) "is" else "are"),
      class = "latents_bad_categorical", call = NULL))
  }
  union(setdiff(vars, c(extra[["ordinal"]], extra[["count"]])), categorical)
}
