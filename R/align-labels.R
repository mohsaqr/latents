#' Number one fit's profiles like another's
#'
#' Profiles and latent classes have no natural order, so two fits of the
#' same model to different samples, waves or subgroups usually number them
#' differently. `align_labels()` renumbers the profiles of `fit`, and its
#' group classes in a two-level model, so that each takes the label of the
#' most similar profile (group class) of `reference`. Nothing is
#' re-estimated: the log likelihood, the estimates and every assignment are
#' the same, under new labels.
#'
#' Profiles are matched on their measurement parameters (means, variances
#' and, when both fits estimate them, correlations; response probabilities;
#' ordinal and count parameters), measured in the reference's units, with
#' each profile's share as a tie-breaker. Every permutation is tried for up
#' to seven profiles, and each reference profile takes its nearest unclaimed
#' match above that. Group classes are then matched on the profile mixture
#' each carries. A membership-covariate fit may be aligned to the same model
#' without covariates: its share is the mean posterior probability, and its
#' membership logits are re-expressed against the new last profile. This is
#' the alignment the cluster bootstrap, [sensitivity()] and
#' [pool_imputations()] use.
#'
#' @param fit A fit from [lpa()], [lca()], [multilpa()] or [multilca()],
#'   with or without membership covariates.
#' @param reference A fit of the same kind whose labels `fit` takes: the
#'   same indicators, coded the same way, and the same numbers of profiles
#'   and group classes.
#' @return `fit`, of the same class, with its profiles (and group classes)
#'   renumbered: every table of [get_results()], the assignments and the
#'   posterior probabilities follow the new labels. Raises
#'   `latents_bad_argument` when either object is not a profile or class
#'   model, or when the two fits differ in their indicators, how they are
#'   coded, or their numbers of profiles or group classes.
#' @seealso [compare_models()] to compare fits on their criteria.
#' @examples
#' first_days <- subset(student_esm, day <= 3)
#' last_days <- subset(student_esm, day > 3)
#' feelings <- c("happy", "exhausted")
#' early <- lpa(first_days, feelings, n_profiles = 2, n_starts = 2, seed = 1)
#' late <- lpa(last_days, feelings, n_profiles = 2, n_starts = 2, seed = 1)
#' # The later days' fit numbers the same two profiles the other way round.
#' get_results(early, "profiles")
#' get_results(late, "profiles")
#' get_results(align_labels(late, early), "profiles")
#' @export
align_labels <- function(fit, reference) {
  supported <- c("multilpa", "multilpa_covariates")
  invisible(Map(function(x, name) {
    if (!inherits(x, supported)) {
      stop(errorCondition(sprintf(paste(
        "`%s` must be a profile or class model from lpa(), lca(), multilpa()",
        "or multilca(), but it is %s."), name, .latents_describe(x)),
        class = "latents_bad_argument", call = NULL))
    }
  }, list(fit, reference), c("fit", "reference")))
  layout <- function(x) {
    list(profiles = x$n_profiles, group_classes = x$n_group_classes,
         continuous = colnames(x$means),
         categorical = lapply(x$response_probabilities %||% list(), colnames),
         ordinal = x$ordinal %||% character(), count = x[["count"]] %||% character(),
         count_dispersion = length(x$count_dispersion) > 0L)
  }
  differs <- !vapply(names(layout(fit)), function(part) {
    identical(layout(fit)[[part]], layout(reference)[[part]])
  }, logical(1))
  if (any(differs)) {
    described <- c(profiles = "numbers of profiles", group_classes = "numbers of group classes",
                   continuous = "continuous indicators",
                   categorical = "categorical indicators or their categories",
                   ordinal = "ordinal indicators", count = "count indicators",
                   count_dispersion = "count models")
    stop(errorCondition(sprintf(paste(
      "`fit` and `reference` must be the same model of the same indicators to",
      "be aligned; they differ in their %s."),
      paste(described[names(differs)[differs]], collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  .multilpa_align_labels(fit, reference)
}
