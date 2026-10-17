#' Does the solution survive a different seed?
#'
#' A latent profile model is fitted by EM from random starts, so what comes back
#' is a local maximum reached from one particular set of them. `n_starts`
#' reports how many starts *within a single seed's stream* reached the best
#' likelihood; it cannot say whether a different stream would have found a
#' different mode. This verb refits the model under each of several seeds and
#' reports what changed.
#'
#' Three things are worth reading off the result. Whether every seed reached the
#' same maximised log likelihood, which is what `optimum` counts. Whether the
#' seeds that reached it converged. And how many observations were assigned to a
#' different profile, which is what `agreement` measures, after the arbitrary
#' profile labels of each refit have been matched to the reference fit's -- two
#' fits of the same mixture can be identical and still number their profiles
#' differently, so comparing labels directly would report disagreement that is
#' not there.
#'
#' A seed whose refit fails contributes a row with `NA` estimates rather than
#' being dropped silently, and a `latents_sensitivity_dropped` warning names
#' how many failed, so the table is never quietly shorter than `seeds`.
#'
#' @param x The reference fit: a [multilpa()] model, an [lta()] model, or a
#'   growth mixture ([mixture_regression()] with `random`).
#'   A [fit_staged()] result varies the second-stage starts while holding its
#'   first-stage measurement solution fixed. A directly fitted model with
#'   measurement blocks held is refused: its original free starting
#'   values are not retained, so replaying the same seed would not replay the
#'   same search. A membership-covariate fit is refused with
#'   `latents_unsupported_sensitivity`. For an [lta()] fit, `agreement` is
#'   the share of occasions given the same profile, after matching the
#'   refit's profiles on the measurement; for a growth mixture it is the
#'   share of persons given the same class, after matching the refit's
#'   classes by their posterior memberships of the fitted persons.
#' @param data Optional for a [multilpa()] fit or a growth mixture, which carry
#'   the columns they were built from, so it is only needed to override them;
#'   required for an [lta()] fit (`latents_bad_argument` otherwise), whose
#'   refits need the occasions and covariates. It must be the data the model
#'   was fitted to, in its original row order; shared identifier and indicator
#'   columns are checked against the fit before assignment agreement is
#'   calculated (`latents_bad_inference_data`).
#' @param seeds Seeds to refit under: at least two distinct whole numbers of
#'   either sign, as `set.seed()` accepts. Anything else raises
#'   `latents_bad_argument`. The reference fit's own seed may
#'   be among them, in which case that row reproduces it and is the arithmetic
#'   check that the refit is the same model.
#' @param n_starts Starts per refit. Defaults to the number the reference fit
#'   used, so each seed gets the same search effort it did.
#' @param max_iter,tol Passed to each refit; default to the reference fit's
#'   (for an [lta()] fit, to the defaults of [lta()]).
#' @param tolerance Two log likelihoods within this of each other count as the
#'   same optimum. Defaults to the reference fit's own replication tolerance,
#'   which is what it used to decide whether its starts agreed (for a growth
#'   mixture `1e-6 * (1 + abs(logLik))`).
#' @return A base `data.frame`, one row per seed, with columns `seed`,
#'   `log_likelihood`, `converged`, `iterations`, `optimum` (1 for the best
#'   maximum found across the seeds, 2 for the next distinct one, and so on),
#'   `best` (whether this seed reached optimum 1) and `agreement` (the
#'   proportion of observations given the same profile as the reference fit,
#'   after label alignment). A failed refit has `NA` throughout except `seed`.
#' @seealso [enumerate_classes()] to vary the number of classes rather than the
#'   seed, and [parameter_inference()] for uncertainty about the parameters of a
#'   single solution.
#' @references Hipp, J. R., & Bauer, D. J. (2006). Local solutions in the
#'   estimation of growth mixture models. *Psychological Methods, 11*, 36--53.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = c(rnorm(60, -1.5), rnorm(60, 1.5)),
#'   score_b = c(rnorm(60, -1), rnorm(60, 1))
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 3, seed = 1)
#' sensitivity(fit, seeds = 1:3, n_starts = 3)
#'
#' # Growth mixtures are known for local solutions (Hipp and Bauer 2006).
#' growth <- mixture_regression(score ~ wave, data = subset(growth_scores, student <= 60),
#'                              n_classes = 2, id = "student", class_level = "group",
#'                              random = ~ 1, n_starts = 1, seed = 1)
#' sensitivity(growth, seeds = 1:2)
#' @export
sensitivity <- function(x, data = NULL, seeds = 1:10, n_starts = NULL,
                        max_iter = NULL, tol = NULL, tolerance = NULL) {
  # Named families first: a covariate fit is class `multilpa`, so the general
  # path below would accept it and then refit it as the wrong model.
  if (inherits(x, "multilpa_covariates")) {
    stop(errorCondition(paste(
      "Seed sensitivity is not implemented for membership-covariate fits:",
      "refitting this family needs arguments the shared refit does not carry."),
      class = "latents_unsupported_sensitivity", call = NULL))
  }
  stopifnot(
    "`x` must be a fitted `multilpa` model, an `lta()` fit or a growth mixture" =
      inherits(x, c("multilpa", "multilpa_transitions", "multilpa_lta",
                    "latents_growth_mixture")),
    "`n_starts` must be a single positive whole number, or NULL" =
      is.null(n_starts) || (is.numeric(n_starts) && length(n_starts) == 1L &&
        is.finite(n_starts) && n_starts >= 1 &&
        n_starts == floor(n_starts) && n_starts <= .Machine$integer.max),
    # Checked here, not left to each refit: a refit's refusal is caught and
    # recorded as that seed's failure, so a mistyped control would come back
    # as a table of NA rows and a warning instead of an error.
    "`max_iter` must be a single non-negative whole number, or NULL" =
      is.null(max_iter) || (is.numeric(max_iter) && length(max_iter) == 1L &&
        is.finite(max_iter) && max_iter >= 0 && max_iter == floor(max_iter) &&
        max_iter <= .Machine$integer.max),
    "`tol` must be a single positive number, or NULL" =
      is.null(tol) || (is.numeric(tol) && length(tol) == 1L &&
                         is.finite(tol) && tol > 0),
    "`tolerance` must be a single positive number, or NULL" =
      is.null(tolerance) || (is.numeric(tolerance) && length(tolerance) == 1L &&
                               is.finite(tolerance) && tolerance > 0))
  if (length(seeds) < 2L || !all(.multilpa_is_seed(seeds)) ||
      anyDuplicated(seeds) > 0L) {
    stop(errorCondition(sprintf(paste(
      "`seeds` must be at least two distinct whole numbers between -%d and %d,",
      "the seeds set.seed() accepts."),
      .Machine$integer.max, .Machine$integer.max),
      class = "latents_bad_argument", call = NULL))
  }
  family <- if (inherits(x, "latents_growth_mixture")) {
    .sensitivity_growth(x, data, n_starts, max_iter, tol)
  } else if (inherits(x, c("multilpa_transitions", "multilpa_lta"))) {
    .sensitivity_transitions(x, data, n_starts, max_iter, tol)
  } else {
    .sensitivity_profiles(x, data, n_starts, max_iter, tol)
  }

  refits <- lapply(seeds, function(seed) {
    tryCatch(family$refit(seed), error = function(condition) condition)
  })
  failed <- vapply(refits, inherits, logical(1), "condition")
  if (any(failed)) {
    warning(warningCondition(
      sprintf("%d of %d seeds did not produce a fit; their rows are NA. First reason: %s",
              sum(failed), length(seeds),
              conditionMessage(refits[[which(failed)[1L]]])),
      class = "latents_sensitivity_dropped", call = NULL))
  }
  likelihood <- vapply(refits, function(fit)
    if (inherits(fit, "condition")) NA_real_ else fit$log_likelihood, numeric(1))
  converged <- vapply(refits, function(fit)
    if (inherits(fit, "condition")) NA else isTRUE(fit$converged), logical(1))
  iterations <- vapply(refits, function(fit)
    if (inherits(fit, "condition")) NA_integer_ else as.integer(fit$iterations),
    integer(1))
  agreement <- vapply(refits, function(fit) {
    if (inherits(fit, "condition")) NA_real_ else family$agreement(fit)
  }, numeric(1))
  optimum <- .multilpa_optimum_index(likelihood, tolerance %||% family$tolerance)

  data.frame(
    seed = as.integer(seeds),
    log_likelihood = likelihood,
    converged = converged,
    iterations = iterations,
    optimum = optimum,
    best = optimum == 1L,
    agreement = agreement,
    row.names = NULL, stringsAsFactors = FALSE)
}

#' How sensitivity() refits a profile model and compares its assignments
#' @return A list with `refit` (a function of a seed), `agreement` (a
#'   function of a refit) and `tolerance`.
#' @noRd
.sensitivity_profiles <- function(x, data, n_starts, max_iter, tol) {
  if (length(x$fixed %||% character()) > 0L && !isTRUE(x$staged)) {
    stop(errorCondition(paste(
      "Seed sensitivity for a directly fixed fit cannot replay its original",
      "free starting values. Use a joint fit or a fit_staged() result, whose",
      "conditional second-stage search can be repeated."),
      class = "latents_unsupported_sensitivity", call = NULL))
  }
  frame <- .multilpa_resolve_data(x, data)
  .multilpa_check_row_count(x, frame)
  .multilpa_check_alignment(x, frame)
  shared <- .multilpa_refit_arguments(x)
  if (isTRUE(x$staged)) {
    shared$start <- .multilpa_stage_start(x)
    shared$fixed <- x$fixed
  }
  shared$n_starts <- as.integer(n_starts %||% max(nrow(x$starts), 1L))
  shared$max_iter <- max_iter %||% x$max_iter %||% 1000L
  shared$tol <- tol %||% x$tol %||% 1e-8
  list(
    refit = function(seed) do.call(multilpa, c(list(data = frame), shared,
                                               list(seed = seed))),
    agreement = function(fit) {
      mean(.multilpa_align_labels(fit, x)$subject_profiles == x$subject_profiles)
    },
    tolerance = x$replication_tolerance %||% 1e-4)
}

#' How sensitivity() refits a transition model and compares its assignments
#'
#' Each occasion's modal profile is compared, after the refit's profiles are
#' matched to the reference's on the measurement (`.lta_profile_signature()`).
#' @noRd
.sensitivity_transitions <- function(x, data, n_starts, max_iter, tol) {
  if (is.null(data)) {
    stop(errorCondition(paste(
      "Seed sensitivity for a transition model needs `data`, the frame it was",
      "fitted to: each refit keeps the persons' occasions and covariates."),
      class = "latents_bad_argument", call = NULL))
  }
  .lta_check_data(x, data)
  modal <- function(fit) {
    if (inherits(fit, "multilpa_lta")) {
      max.col(fit$expectation$subject_posteriors, ties.method = "first")
    } else fit$subject_profiles
  }
  reference <- modal(x)
  starts <- as.integer(n_starts %||% max(nrow(x$starts), 1L))
  list(
    refit = function(seed) {
      .lta_refit(x, data, starts, max_iter %||% 5000L, tol %||% 1e-10, seed = seed)
    },
    agreement = function(fit) {
      order <- .multilpa_match_order(.lta_profile_signature(x, x),
                                     .lta_profile_signature(fit, x))
      mean(match(modal(fit), order) == reference)
    },
    tolerance = x$replication_tolerance %||% 1e-4)
}

#' How sensitivity() refits a growth mixture and compares its assignments
#'
#' Persons' modal classes are compared after the refit's classes are matched
#' to the reference's by their posterior memberships of the fitted persons.
#' @noRd
.sensitivity_growth <- function(x, data, n_starts, max_iter, tol) {
  if (!is.null(data)) .trajectory_check_data(x, data)
  spec <- x$spec
  frame <- spec$model_data
  weight_column <- NULL
  if (.latents_is_weighted(x)) {
    weight_column <- ".latents_sensitivity_weight"
    frame[[weight_column]] <- spec$sampling_weights[spec$group_index]
  }
  arguments <- c(.growth_refit_arguments(x, weight_column),
                 list(n_starts = as.integer(n_starts %||% x$settings$n_starts),
                      max_iter = max_iter %||% x$settings$max_iter,
                      tol = tol %||% x$settings$tol))
  reference <- max.col(x$expectation$posterior, ties.method = "first")
  list(
    refit = function(seed) do.call(mixture_regression,
                                   c(list(data = frame), arguments, list(seed = seed))),
    agreement = function(fit) {
      matched <- .growth_match_classes(x, fit)
      mean(match(max.col(matched$posterior, ties.method = "first"), matched$order) ==
             reference)
    },
    tolerance = 1e-6 * (1 + abs(x$log_likelihood)))
}

#' Number the distinct maxima a set of log likelihoods reached
#'
#' Sorted from the best downwards, a new optimum begins wherever the drop from
#' the previous one exceeds `tolerance`. Comparing the likelihoods pairwise
#' would not do: `==` on doubles is never right here, and a chain of values each
#' within tolerance of the last is one basin, not several.
#'
#' @param likelihood Maximised log likelihoods, possibly with `NA`.
#' @param tolerance Two values within this count as the same optimum.
#' @return An integer vector the same length as `likelihood`; `1` is the best
#'   maximum, and a failed fit is `NA`.
#' @noRd
.multilpa_optimum_index <- function(likelihood, tolerance) {
  stopifnot(is.numeric(likelihood),
            length(tolerance) == 1L, is.finite(tolerance), tolerance > 0)
  index <- rep(NA_integer_, length(likelihood))
  observed <- which(!is.na(likelihood))
  if (length(observed) == 0L) return(index)
  ranked <- observed[order(likelihood[observed], decreasing = TRUE)]
  values <- likelihood[ranked]
  index[ranked] <- cumsum(c(1L, as.integer(diff(values) < -tolerance)))
  index
}
