utils::globalVariables(c("seed", "status"))

#' Assess sensitivity of a fitted mixture model to initialisation
#'
#' @description
#' `sensitivity()` evaluates how consistently a fitted mixture model is
#' recovered from different random initialisations. Mixture-model estimation
#' can converge to different local maxima, yielding different class assignments
#' for the same data and model specification. Examining this dependence on
#' initialisation helps assess the numerical reliability of the solution
#' being interpreted.
#'
#' The function repeats estimation under several random seeds and compares
#' the resulting likelihood values and class assignments with the reference
#' fit. Class labels are aligned before assignments are compared. By
#' considering likelihood and classification agreement together, the analysis
#' shows whether repeated searches recover a comparable fit and whether they
#' assign the same observations to the same classes.
#'
#' Results include convergence status, log likelihood, and classification
#' agreement for each refit. Similar likelihood values are grouped using a
#' specified tolerance. The data and model specification remain fixed
#' throughout the analysis.
#'
#' @section Interpretation:
#' Consistent results across converged refits support numerical
#' reproducibility, but do not guarantee a global maximum or establish the
#' validity of the model. Similar likelihood values alone do not imply
#' identical parameter estimates. The analysis assesses sensitivity to
#' initialisation, rather than uncertainty arising from sampling or model
#' selection.
#'
#' @details
#' Replication across seeds complements replication across starts within one
#' fit. Each refit uses the reference fit's estimation controls unless
#' alternatives are supplied.
#'
#' A seed whose refit fails gives a row with `NA` estimates, and a
#' `latents_sensitivity_dropped` warning reports how many seeds failed.
#'
#' @param x The reference fit: a [multilpa()] model, an [lta()] model, or a
#'   growth mixture ([mixture_regression()] with `random`). A [fit_staged()]
#'   result varies the second-stage starts while holding its first-stage
#'   measurement solution fixed. A directly fitted model with measurement
#'   blocks held fixed is refused (`latents_unsupported_sensitivity`), because
#'   its original free starting values are not retained. A
#'   membership-covariate fit raises
#'   `latents_unsupported_sensitivity`. For an [lta()] fit, `agreement` is the
#'   share of occasions given the same profile, after matching the refit's
#'   profiles on the measurement; for a growth mixture it is the share of
#'   persons given the same class, after matching the refit's classes by
#'   their posterior memberships of the fitted persons.
#' @param data Optional for a [multilpa()] fit or a growth mixture, which
#'   carry the columns they were built from. Required for an [lta()] fit
#'   (`latents_bad_argument` otherwise). It must be the data the model was
#'   fitted to, in its original row order; shared identifier and indicator
#'   columns are checked against the fit (`latents_bad_inference_data`).
#' @param seeds At least two distinct integer seeds accepted by `set.seed()`.
#'   Default `1:10`. Invalid values raise `latents_bad_argument`. Reusing the
#'   reference seed reproduces its initialisation when the other fitting
#'   controls are unchanged.
#' @param n_starts Starts per refit. Defaults to the number the reference fit
#'   used.
#' @param max_iter,tol Passed to each refit. Default to the reference fit's
#'   (for a fit that does not record them, to the defaults of the function
#'   that fitted it).
#' @param tolerance Maximum log-likelihood difference for grouping refits.
#'   Each group is anchored at its highest log likelihood
#'   and spans no more than `tolerance`. Defaults to the reference fit's
#'   replication tolerance, 1e-6 per observed value.
#' @param cores Number of processor cores (default 1). The seeds are refitted
#'   in parallel, one seed per worker at a time, each refit on one core. Every
#'   refit sets its own seed, so the result is identical to `cores = 1`.
#'   Workers are forked processes on macOS and Linux, and a socket cluster on
#'   Windows and on macOS with Apple's Accelerate BLAS, as in [multilpa()];
#'   socket workers load the installed copy of the package. More cores than
#'   seeds are not used.
#' @return A `data.frame` of class `latents_sensitivity` (which `plot()`
#'   draws), one row per seed, with columns `seed`,
#'   `log_likelihood`, `converged`, `iterations`, `optimum` (1 for the group
#'   containing the highest attained log likelihood, 2 for the next lower
#'   group, and so on),
#'   `best` (whether this seed reached optimum 1) and `agreement` (the
#'   proportion of classified observations given the same profile as the
#'   reference fit, after label matching). For profile and transition models,
#'   this uses observations included in both fits under their `drop_empty`
#'   settings. A failed refit has `NA` in every column except
#'   `seed`. The `optimum` and `best` columns group attained likelihood values
#'   without filtering on convergence; inspect `converged` before interpreting
#'   them as evidence of replicated maxima.
#'
#'   `as.data.frame()` returns a plain `data.frame`. `plot()` accepts `main`
#'   and `subtitle` and returns a ggplot object showing each refit's
#'   log-likelihood difference from the best converged refit and its
#'   classification agreement with the reference fit. Unfinished
#'   refits are distinguished from converged maxima; if none converged, the
#'   highest observed log likelihood is the reference instead.
#' @seealso [enumerate_classes()] to vary the number of classes, and
#'   [parameter_inference()] for uncertainty about the parameters of a single
#'   solution.
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
                        max_iter = NULL, tol = NULL, tolerance = NULL,
                        cores = 1L) {
  cores <- .latents_check_cores(cores)
  # Named families first: a covariate fit is class `multilpa`, so the general
  # path below would accept it and then refit it as the wrong model.
  if (inherits(x, "multilpa_covariates")) {
    stop(errorCondition(paste(
      "Seed sensitivity is not implemented for membership-covariate fits:",
      "refitting this family needs arguments the shared refit does not carry."),
      class = "latents_unsupported_sensitivity", call = NULL))
  }
  .latents_require_fit(x, substitute(x),
                       c("multilpa", "multilpa_transitions", "multilpa_lta",
                         "latents_growth_mixture"),
                       "a fitted profile, class, transition or growth mixture model")
  .latents_require(
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

  # Seeds, not starts, are spread over the workers: every family refits this
  # way (a growth mixture's refit takes no `cores`), each refit seeds itself,
  # and the results come back in the order of `seeds` whoever ran them.
  refits <- .latents_map_candidates(seeds, function(seed) {
    tryCatch(family$refit(seed), error = function(condition) condition)
  }, cores)
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
  tolerance <- tolerance %||% family$tolerance
  optimum <- .multilpa_optimum_index(likelihood, tolerance)

  result <- data.frame(
    seed = as.integer(seeds),
    log_likelihood = likelihood,
    converged = converged,
    iterations = iterations,
    optimum = optimum,
    best = optimum == 1L,
    agreement = agreement,
    row.names = NULL, stringsAsFactors = FALSE)
  # Still a data frame: the class adds plot(), and the tolerance the
  # optimum index used is kept for it.
  attr(result, "tolerance") <- tolerance
  class(result) <- c("latents_sensitivity", "data.frame")
  result
}

#' @noRd
#' @export
as.data.frame.latents_sensitivity <-function(x, row.names = NULL, optional = FALSE, ...) {
  attr(x, "tolerance") <- NULL
  class(x) <- "data.frame"
  x
}

#' @noRd
#' @export
plot.latents_sensitivity <-function(x, main = NULL, subtitle = NULL, ...) {
  .gg_require()
  table <- as.data.frame(x)
  fitted <- is.finite(table$log_likelihood)
  if (!any(fitted)) {
    stop(errorCondition("No refit succeeded, so there is nothing to plot.",
                        class = c("latents_nothing_to_plot", "latents_bad_argument"),
                        call = NULL))
  }
  converged <- fitted & table$converged %in% TRUE
  tolerance <- attr(x, "tolerance") %||% 0
  best <- max(table$log_likelihood[if (any(converged)) converged else fitted])
  reached <- converged & abs(table$log_likelihood - best) <= tolerance
  status <- ifelse(!fitted, "failed",
                   ifelse(!converged, "unfinished",
                           ifelse(reached, "best maximum", "a lower maximum")))
  status <- factor(status, levels = c("best maximum", "a lower maximum", "unfinished", "failed"))
  seed <- factor(table$seed, levels = unique(table$seed))
  likelihood_panel <- if (any(converged)) "Log likelihood below the best converged refit" else
    "Log likelihood below the highest observed value"
  frame <- rbind(
    data.frame(seed = seed, status = status, panel = likelihood_panel,
               value = best - table$log_likelihood),
    data.frame(seed = seed, status = status, panel = "Agreement with the reference fit",
               value = table$agreement))
  frame$panel <- factor(frame$panel, levels = c(likelihood_panel,
                                                "Agreement with the reference fit"))
  frame <- frame[!is.na(frame$value), , drop = FALSE]
  # A failed refit has no value to draw; it keeps its place on the seed axis
  # (the scale does not drop it) and is named in the caption.
  failed_seeds <- table$seed[!fitted]
  unfinished_seeds <- table$seed[fitted & !converged]
  main <- main %||% if (any(converged)) {
    sprintf("%d of %d seeds reached the best converged maximum", sum(reached), nrow(table))
  } else "No seed converged: finite refits are unfinished"
  subtitle <- subtitle %||% paste(
    "Each seed refits the model; the top panel compares its log likelihood",
    "with the reference named on that panel, the bottom its classification agreement.")
  colours <- c("best maximum" = .gg_okabe_ito[[3L]], "a lower maximum" = .gg_okabe_ito[[5L]],
               "unfinished" = .gg_okabe_ito[[2L]],
               "failed" = .gg_reserved_grey)
  shapes <- c("best maximum" = 21L, "a lower maximum" = 24L, "unfinished" = 23L, "failed" = 4L)
  plot <- ggplot2::ggplot(frame, ggplot2::aes(seed, value, colour = status, fill = status,
                                              shape = status)) +
    ggplot2::geom_point(size = 2.6, stroke = 0.6) +
    ggplot2::facet_wrap(~panel, ncol = 1L, scales = "free_y") +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE, name = NULL) +
    ggplot2::scale_fill_manual(values = colours, drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE, name = NULL) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(title = main, subtitle = subtitle, x = "Seed", y = NULL,
                  caption = paste(c(
                    if (any(converged)) sprintf(
                      "Only converged refits within %.3g of the best converged log likelihood count as reaching its maximum.",
                      tolerance),
                    if (length(unfinished_seeds) > 0L) sprintf(
                      "Refit unfinished for seed%s %s; its finite values are shown, not counted as maxima.",
                      if (length(unfinished_seeds) > 1L) "s" else "",
                      paste(unfinished_seeds, collapse = ", ")),
                    if (length(failed_seeds) > 0L) sprintf(
                      "Refit failed (no value to draw) for seed%s %s.",
                      if (length(failed_seeds) > 1L) "s" else "",
                      paste(failed_seeds, collapse = ", "))), collapse = "\n")) +
    .gg_theme() +
    ggplot2::theme(legend.position = "bottom")
  plot
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
  shared$n_starts <- as.integer(n_starts %||% x$n_starts %||% max(nrow(x$starts), 1L))
  shared$n_screen <- .latents_refit_controls(x, n_starts = shared$n_starts)$n_screen
  shared$max_iter <- max_iter %||% x$max_iter %||% 1000L
  shared$tol <- tol %||% x$tol %||% 1e-8
  list(
    refit = function(seed) do.call(multilpa, c(list(data = frame), shared,
                                               list(seed = seed))),
    agreement = function(fit) {
      included <- .multilpa_informative_rows(x) & .multilpa_informative_rows(fit)
      mean(.multilpa_align_labels(fit, x)$subject_profiles[included] ==
             x$subject_profiles[included])
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
  starts <- as.integer(n_starts %||% x$n_starts %||% max(nrow(x$starts), 1L))
  list(
    # The reference fit's own controls, as for a profile fit: with lta()'s
    # defaults hard-coded instead, a fit made with another `tol` or
    # `max_iter` did not reproduce under its own seed (INF-4). lta() fits
    # record both; the defaults of lta() remain for one that does not.
    refit = function(seed) {
      .lta_refit(x, data, starts, max_iter %||% x$max_iter %||% 5000L,
                 tol %||% x$tol %||% 1e-10, seed = seed)
    },
    agreement = function(fit) {
      order <- .multilpa_match_order(.lta_profile_signature(x, x),
                                     .lta_profile_signature(fit, x))
      included <- .multilpa_informative_rows(x) & .multilpa_informative_rows(fit)
      mean(match(modal(fit)[included], order) == reference[included])
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
    tolerance = x$replication_tolerance %||% (1e-6 * (1 + abs(x$log_likelihood))))
}

#' Number the distinct maxima a set of log likelihoods reached
#'
#' Sorted from the best downwards, each optimum is anchored at its best value,
#' and a new optimum begins at the first value more than `tolerance` below the
#' current anchor; that value becomes the next anchor. Every optimum therefore
#' spans at most `tolerance`, so two values more than `tolerance` apart are
#' never one optimum. Chaining from neighbour to neighbour did not have that
#' property: thirty values 0.9e-4 apart, spanning 2.6e-3, were one optimum at a
#' tolerance of 1e-4 (review 2026-10-26 round 2, INF2-5). The rule is the one
#' the equivalence harness uses (`seed_stability()`), read from the top.
#'
#' @param likelihood Maximised log likelihoods, possibly with `NA`.
#' @param tolerance How far below an optimum's best value a log likelihood may
#'   lie and still count as that optimum.
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
  # The running state is (anchor, count); `accumulate` keeps the count each
  # value was given. Reduce() over nothing returns `init` bare, not in a list.
  if (length(values) == 1L) {
    index[ranked] <- 1L
    return(index)
  }
  counts <- Reduce(function(state, value) {
    if (state[[1L]] - value > tolerance) c(value, state[[2L]] + 1) else state
  }, values[-1L], init = c(values[[1L]], 1), accumulate = TRUE)
  index[ranked] <- as.integer(vapply(counts, `[[`, numeric(1), 2L))
  index
}
