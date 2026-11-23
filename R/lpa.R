# Single-level latent profile and latent class analysis.
#
# The package's models are two-level: observations nested in groups, with
# latent classes at both levels. A single-level model is the special case with
# one observation per unit and no group classes. `multilpa(id = NULL)` fits it;
# these two verbs name it, so that the ordinary analysis reads as what it is.

#' Fit a latent profile model that groups observations by their pattern on
#' continuous indicators
#'
#' @description
#' `lpa()` fits a latent profile model, which describes heterogeneity in
#' continuous indicators through a specified number of latent profiles. Each
#' profile represents a pattern of indicator means and within-profile variation,
#' estimated as a finite mixture of multivariate normal distributions.
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
#' @param ... Further arguments, passed to [multilpa()], which documents all
#'   of them. The ones used most often:
#'   * `variance_model`: `"varying"` (default) gives each profile its own
#'     variances; `"equal"` gives all profiles the same variances.
#'   * `covariance_model`: `"diagonal"` (default) treats the indicators as
#'     uncorrelated within a profile; `"full"` estimates their within-profile
#'     covariances.
#'   * `model`: the covariance structure as one of mclust's fourteen codes
#'     (`"EII"`, `"VII"`, `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`,
#'     `"VEE"`, `"EVE"`, `"VVE"`, `"EEV"`, `"VEV"`, `"EVV"`, `"VVV"`),
#'     instead of `variance_model` and `covariance_model`. The default
#'     structure is `"VVI"`.
#'   * `missing`: `"error"` (default) stops when an indicator has a missing
#'     value; `"fiml"` keeps every row and uses the values it has.
#'   * `categorical`, `ordinal`, `count`: indicators in `vars` that are not
#'     continuous.
#'   * `n_starts` (default 10) and `seed`: the number of starts and the seed
#'     that makes them reproducible.
#'   * `boundary`: `"allow"` (default) or `"exclude"`, whether a solution
#'     with a variance at its lower bound may be reported.
#'   * `prior`, `noise` and `profile_covariates`: a prior on the Gaussian
#'     parameters, a noise component for outliers, and covariates that
#'     predict profile membership.
#'
#'   `id` and `n_group_classes` are not accepted: they belong to the
#'   two-level model.
#' @return A fitted model of class `multilpa`, or `multilpa_covariates` when
#'   `profile_covariates` is given. Printing it shows the log likelihood,
#'   convergence, profile sizes and profile means; `summary()` prints every
#'   table. [get_results()] returns the tables by name, among them:
#'   * `"profiles"` (also `as.data.frame(fit)`): one row per profile and
#'     indicator, with the mean, variance and their standard errors;
#'   * `"profile_probabilities"`: the proportion of each profile;
#'   * `"assignments"`: the data with each row's most likely profile and its
#'     posterior probabilities;
#'   * `"classification"` and `"entropy"`: how clearly the rows are
#'     classified;
#'   * `"information_criteria"`: the log likelihood, AIC, BIC and related
#'     criteria;
#'   * `"starts"`: the log likelihood and convergence of every start.
#'
#'   `get_results(fit, "all")` returns every table. The fit's `call` records
#'   the original `lpa()` call.
#' @section Conditions:
#'   Errors: `latents_bad_argument` when `id` or `n_group_classes` is passed
#'   or an argument has a value that is not accepted; `latents_bad_data` when
#'   a column is missing from `data` or an indicator has missing values under
#'   `missing = "error"`; `latents_unidentified` when there are more profiles
#'   than distinct rows. Warnings: `latents_unconverged` when the reported
#'   start did not converge, `latents_boundary` when a variance is at its
#'   lower bound, and `latents_small_classes` when a profile is nearly empty.
#'   [multilpa()] raises the others, and [latents-conditions] lists them all.
#' @seealso [lca()] for categorical indicators; [multilpa()] for the two-level
#'   model; [enumerate_lpa()] to compare numbers of profiles.
#' @examples
#' fit <- lpa(srl, c("cognitive_strategies", "intrinsic_value",
#'                   "self_efficacy", "self_regulation", "test_anxiety"),
#'            n_profiles = 2, variance_model = "equal",
#'            covariance_model = "full", n_starts = 3, seed = 1)
#' fit
#' get_results(fit, "profile_probabilities")
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

#' Fit a latent class model that groups observations by their pattern of
#' categorical responses
#'
#' @description
#' `lca()` fits a latent class model, which describes heterogeneity in
#' categorical response patterns through a specified number of latent classes.
#' Each class is characterized by the probabilities of the response categories
#' for each indicator, which describe how response patterns differ between
#' classes.
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
#' @param vars Names of the discrete indicator columns. Numeric, integer,
#'   logical, character and factor columns are accepted, and every distinct
#'   value of a column is one category. Columns named in `ordinal` or `count`
#'   (through `...`) are modelled as such; the rest are categorical.
#' @param n_classes Number of classes; a positive integer.
#' @param ... Further arguments, passed through [multilca()] to [multilpa()],
#'   which documents all of them. The ones used most often:
#'   * `missing`: `"error"` (default) stops when an indicator has a missing
#'     value; `"fiml"` keeps every row and uses the responses it has.
#'   * `ordinal`, `count`: indicators in `vars` modelled as ordered categories
#'     or as counts instead of as unordered categories.
#'   * `local_dependence`: pairs of indicators, written `"item1:item2"`,
#'     that stay associated within a class.
#'   * `n_starts` (default 10) and `seed`: the number of starts and the seed
#'     that makes them reproducible.
#'   * `min_probability`: the lower bound on a response probability (default
#'     `1e-10`).
#'   * `profile_covariates`: covariates that predict class membership.
#'
#'   `id` and `n_group_classes` are not accepted (they belong to the
#'   two-level model), and neither is `n_profiles`: the number of classes is
#'   `n_classes`.
#' @return A fitted model of class `multilpa`, or `multilpa_covariates` when
#'   `profile_covariates` is given. The classes are numbered 1, 2, and so on
#'   in the `profile` column of its tables. [get_results()] returns the
#'   tables by name, among them:
#'   * `"responses"` (also `as.data.frame(fit)`): one row per class,
#'     indicator and category, with the response probability and its
#'     standard error;
#'   * `"profile_probabilities"`: the proportion of each class;
#'   * `"assignments"`: the data with each row's most likely class and its
#'     posterior probabilities;
#'   * `"classification"` and `"entropy"`: how clearly the rows are
#'     classified;
#'   * `"residuals"`: bivariate residuals, which show pairs of indicators
#'     the classes do not explain;
#'   * `"information_criteria"` and `"starts"`.
#'
#'   The fit's `call` records the original `lca()` call.
#' @section Conditions:
#'   Errors: `latents_bad_argument` when `id`, `n_group_classes` or
#'   `n_profiles` is passed; `latents_bad_data` when a column is missing from
#'   `data` or an indicator has missing values under `missing = "error"`;
#'   `latents_bad_categorical` when `categorical` names a column not in
#'   `vars`. Warnings: `latents_unidentified` when the model has more free
#'   parameters than the table of response patterns has free cells,
#'   `latents_unconverged` when the reported start did not converge, and
#'   `latents_small_classes` when a class is nearly empty. [multilpa()]
#'   raises the others, and [latents-conditions] lists them all.
#' @seealso [lpa()] for continuous indicators; [multilca()] for the
#'   two-level model; [enumerate_lca()] to compare numbers of classes.
#' @examples
#' fit <- lca(student_esm, c("time_with_friends", "on_social_media",
#'                           "tv_video_games", "sports"),
#'            n_classes = 2, n_starts = 3, seed = 1)
#' fit
#' get_results(fit, "profile_probabilities")
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

#' Fit latent profile models with different numbers of profiles and covariance
#' structures, and compare them
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
#'   integers. Default `1:4`.
#' @param model The covariance structures to fit at every number of
#'   profiles; each structure with each number of profiles is one candidate.
#'   * `"basic"` (default): four structures. `"EEI"`: the profiles share
#'     their variances and the indicators are uncorrelated within a profile.
#'     `"VVI"`: each profile has its own variances, still uncorrelated.
#'     `"EEE"`: the profiles share one covariance matrix. `"VVV"`: each
#'     profile has its own covariance matrix.
#'   * `"all"`: all fourteen structures of mclust.
#'   * One or more of the codes `"EII"`, `"VII"`, `"EEI"`, `"VEI"`, `"EVI"`,
#'     `"VVI"`, `"EEE"`, `"VEE"`, `"EVE"`, `"VVE"`, `"EEV"`, `"VEV"`,
#'     `"EVV"`, `"VVV"`: those structures only. [multilpa()] explains the
#'     three letters.
#'   * `NULL`: no structures are crossed, and every candidate has the
#'     structure [multilpa()] fits by default (`"VVI"`).
#'
#'   Giving `variance_model`, `covariance_model`, `volume`, `shape` or
#'   `orientation` in `...` sets one structure for every candidate; `model`
#'   must then be left out.
#' @param ... Further arguments for [multilpa()], used for every candidate,
#'   such as `missing` (`"error"` or `"fiml"`), `n_starts`, `seed`,
#'   `boundary` (`"allow"` or `"exclude"`) and `categorical`. `id` and
#'   `n_group_classes` are not accepted.
#' @return An object of class `multilpa_enumeration`, as [enumerate_classes()]
#'   returns. [get_candidates()] (also [as.data.frame()]) gives one row per
#'   candidate with its log likelihood, information criteria, entropy,
#'   convergence and replication; [get_best()] and [summary()] the candidate
#'   with the lowest value of each criterion; [plot()] the criteria against
#'   the number of profiles; and [get_candidate()] one fitted candidate, a
#'   `multilpa` fit. A candidate that could not be fitted keeps its row,
#'   with `NA` values.
#' @section Conditions:
#'   `latents_bad_argument` when `id` or `n_group_classes` is passed, when
#'   `model` names an unknown structure, and when `model` is given together
#'   with `variance_model`, `covariance_model`, `volume`, `shape` or
#'   `orientation`. A `latents_single_start` message names the candidates
#'   whose best log likelihood was reached from one start only.
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

#' Fit latent class models with different numbers of classes, and compare them
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
#'   integers. Default `1:4`.
#' @param ... Further arguments for [multilca()], used for every candidate,
#'   such as `missing` (`"error"` or `"fiml"`), `n_starts`, `seed` and
#'   `min_probability`. `id` and `n_group_classes` are not accepted, and
#'   neither is `model`: a latent class model has no covariance structure.
#'   Every indicator is categorical unless named in `ordinal` or `count`.
#'   `categorical` is accepted but has no further effect; its names must be
#'   in `vars`.
#' @return An object of class `multilpa_enumeration`, as [enumerate_classes()]
#'   returns. [get_candidates()] (also [as.data.frame()]) gives one row per
#'   candidate with its log likelihood, information criteria, entropy,
#'   convergence and replication; [get_best()] and [summary()] the candidate
#'   with the lowest value of each criterion; [plot()] the criteria against
#'   the number of classes; and [get_candidate()] one fitted candidate. The
#'   number of classes is in the `n_profiles` column of the candidates table,
#'   and [get_candidate()] selects it through its `n_profiles` argument. A
#'   candidate that could not be fitted keeps its row, with `NA` values.
#' @section Conditions:
#'   `latents_bad_argument` when `id` or `n_group_classes` is passed, or when
#'   `model` is passed (no indicator is continuous); `latents_bad_categorical`
#'   when `categorical` names a column not in `vars`. A candidate with more
#'   free parameters than the table of response patterns has free cells
#'   warns with class `latents_unidentified`. A `latents_single_start`
#'   message names the candidates whose best log likelihood was reached from
#'   one start only.
#' @seealso [lca()] to fit one model; [enumerate_lpa()] for continuous
#'   indicators.
#' @examples
#' models <- enumerate_lca(student_esm, c("time_with_friends",
#'                                        "on_social_media", "tv_video_games",
#'                                        "sports"),
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
