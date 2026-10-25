# Latent transition analysis.
#
# The measurement model, its M-step and its conditional Gaussian moments are
# shared verbatim with multilpa(). What changes is the dependence structure: a
# group's occasions are no longer conditionally independent given its class, so
# the sum of per-observation log marginals that multilpa() forms with rowsum()
# becomes a forward-backward pass over the ordered occasions.

#' Lay observations out as one row per group and one column per occasion
#'
#' Groups differ in length and may skip positions, so the recursion runs on a
#' padded rectangle. Two logical masks distinguish the two kinds of absence: an
#' occasion past a group's last one carries neither an emission nor a
#' transition, while a position skipped inside a group's span carries a
#' transition but no emission.
#'
#' @param group_index Integer group indices in first-occurrence order.
#' @param time_values Position of each observation within its group.
#' @param n_groups Number of observed groups.
#' @param occasions `"observed"` numbers each group's own occasions
#'   consecutively; `"grid"` places them on the grid of every position seen in
#'   the data, so that a skipped position still consumes a transition.
#' @return A list with `slot` (groups by occasions, holding the row index or
#'   `NA`), `within` (whether the occasion lies inside the group's span),
#'   `span`, `n_occasions`, `occasions` and `n_observations`.
#' @noRd
.multilpa_sequence_layout <- function(group_index, time_values, n_groups,
                                      occasions = c("observed", "grid")) {
  occasions <- match.arg(occasions)
  stopifnot(
    "`group_index` and `time_values` must have one value per observation" =
      length(group_index) == length(time_values),
    "`n_groups` must be a single positive integer" =
      length(n_groups) == 1L && n_groups >= 1L)
  groups <- factor(group_index, levels = seq_len(n_groups))
  rank_within <- if (identical(occasions, "grid")) {
    match(time_values, sort(unique(time_values)))
  } else {
    sorted <- order(group_index, time_values)
    ranks <- integer(length(group_index))
    ranks[sorted] <- sequence(tabulate(group_index, nbins = n_groups))
    ranks
  }
  first <- as.integer(tapply(rank_within, groups, min))
  last <- as.integer(tapply(rank_within, groups, max))
  span <- last - first + 1L
  n_occasions <- max(span)
  slot <- matrix(NA_integer_, n_groups, n_occasions)
  slot[cbind(group_index, rank_within - first[group_index] + 1L)] <-
    seq_along(group_index)
  list(slot = slot, within = outer(span, seq_len(n_occasions), ">="),
       span = span, n_occasions = n_occasions, occasions = occasions,
       n_observations = length(group_index))
}

#' Emission log densities laid out by occasion, with the shared offset removed
#'
#' The recursion accumulates sums of log densities, so a common offset must be
#' taken out before log transition probabilities are added to them; otherwise a
#' very small density rounds the priors away, exactly as it would in the
#' cross-sectional E-step. The offset cancels from every posterior and is added
#' back once to the log likelihood.
#'
#' @param validate Check inputs for direct calls; disabled inside validated fits.
#' @param log_density Observation-by-profile matrix of measurement log densities.
#' @param layout Sequence layout from `.multilpa_sequence_layout()`.
#' @return A list with `log_density` (one groups-by-profiles matrix per
#'   occasion, offset removed) and `offset` (one total per group).
#' @noRd
.multilpa_sequence_emission <- function(log_density, layout, validate = TRUE) {
  if (validate) stopifnot(is.matrix(log_density), is.list(layout))
  n_groups <- nrow(layout$slot)
  n_profiles <- ncol(log_density)
  blocks <- lapply(seq_len(layout$n_occasions), function(occasion) {
    rows <- layout$slot[, occasion]
    observed <- !is.na(rows)
    # An occasion a group does not occupy contributes no measurement
    # information, which is a log density of zero for every profile alike.
    block <- matrix(0, n_groups, n_profiles)
    block[observed, ] <- log_density[rows[observed], , drop = FALSE]
    block
  })
  offsets <- lapply(blocks, function(block) .multilpa_row_max(block))
  if (any(!is.finite(unlist(offsets, use.names = FALSE)))) {
    stop("All component densities vanished, or a density overflowed.")
  }
  list(log_density = Map(function(block, offset) block - offset, blocks, offsets),
       offset = Reduce(`+`, offsets))
}

#' Forward and backward recursions over every group at once
#'
#' Groups share the occasion grid, so each step of the recursion is a single
#' matrix update rather than one pass per group. A group whose sequence has
#' ended keeps its state unchanged, which leaves its likelihood and posteriors
#' identical to those of the shorter sequence it actually has.
#'
#' @param emission Offset-removed emission list from `.multilpa_sequence_emission()`.
#' @param layout Sequence layout from `.multilpa_sequence_layout()`.
#' @param validate Check inputs for direct calls; disabled inside validated fits.
#' @param initial Initial profile probabilities for one group class.
#' @param transition Profile transition matrix for one group class, rows from.
#' @return A list with `alpha`, `beta` (one groups-by-profiles matrix per
#'   occasion) and `log_scaled`, the offset-removed sequence log likelihoods.
#' @noRd
.multilpa_forward_backward <- function(emission, layout, initial, transition,
                                        validate = TRUE) {
  if (validate) stopifnot(
    "`initial` and `transition` must be conformable positive probabilities" =
      is.list(emission) && is.numeric(initial) && is.matrix(transition) &&
      nrow(transition) == length(initial) && ncol(transition) == length(initial) &&
      all(initial > 0) && all(transition > 0))
  n_groups <- nrow(layout$slot)
  .latents_scaled_forward_backward(
    emission, layout, matrix(initial, n_groups, length(initial), byrow = TRUE),
    transition)
}

#' Posterior occupancy and transition counts for one group class
#'
#' @param validate Check inputs for direct calls; disabled inside validated fits.
#' @param pass Forward and backward recursions from `.multilpa_forward_backward()`.
#' @param emission Offset-removed emission list.
#' @param layout Sequence layout.
#' @param transition Profile transition matrix for this group class.
#' @param weights Posterior probability that each group belongs to this class.
#' @param n_observations Number of input rows.
#' @return A list with `posterior` (observations by profiles, weighted by
#'   `weights`), `initial` (expected first-occasion counts) and `transition`
#'   (expected from-to counts).
#' @noRd
.multilpa_sequence_moments <- function(pass, emission, layout, transition,
                                       weights, n_observations, validate = TRUE) {
  if (validate) stopifnot(
    "`pass` must be the list `.multilpa_forward_backward()` returns" =
      is.list(pass) && all(c("alpha", "beta", "log_scaled") %in% names(pass)),
    "`weights` must give one group-class posterior per group" =
      is.numeric(weights) && length(weights) == nrow(layout$slot),
    "`transition` must be a square matrix of positive probabilities" =
      is.matrix(transition) && is.numeric(transition) &&
      nrow(transition) == ncol(transition) && all(transition > 0))
  n_profiles <- ncol(transition)
  occupancy <- if (!is.null(pass$scaled_alpha)) {
    Map(function(a, b) a * b * weights, pass$scaled_alpha, pass$scaled_beta)
  } else lapply(seq_len(layout$n_occasions), function(occasion) {
    exp(pass$alpha[[occasion]] + pass$beta[[occasion]] - pass$log_scaled) * weights
  })
  posterior <- matrix(0, n_observations, n_profiles)
  placed <- do.call(rbind, lapply(seq_len(layout$n_occasions), function(occasion) {
    rows <- layout$slot[, occasion]
    observed <- !is.na(rows)
    if (!any(observed)) return(NULL)
    cbind(rows[observed], occupancy[[occasion]][observed, , drop = FALSE])
  }))
  # Every input row occupies exactly one cell of the layout, so the placement
  # is an assignment rather than an accumulation.
  posterior[placed[, 1L], ] <- placed[, -1L, drop = FALSE]
  counts <- if (layout$n_occasions < 2L) {
    matrix(0, n_profiles, n_profiles)
  } else if (!is.null(pass$scaled_alpha)) {
    transition * Reduce(`+`, lapply(seq_len(layout$n_occasions)[-1L], function(t) {
      active <- weights * layout$within[, t]
      crossprod(pass$scaled_alpha[[t - 1L]] * active,
                pass$scaled_beta[[t]] * pass$density[[t]] / pass$scale[, t])
    }))
  } else {
    # Log-space callers (including the general engine) retain the complete
    # pair log probability before exponentiating, avoiding Inf * 0.
    n_groups <- nrow(layout$slot)
    origin <- rep(seq_len(n_profiles), times = n_profiles)
    destination <- rep(seq_len(n_profiles), each = n_profiles)
    log_transition <- rep(log(transition), each = n_groups)
    Reduce(`+`, lapply(seq_len(layout$n_occasions)[-1L], function(occasion) {
      active <- weights * layout$within[, occasion]
      log_from <- pass$alpha[[occasion - 1L]] - pass$log_scaled
      log_to <- pass$beta[[occasion]] + emission[[occasion]]
      joint <- log_from[, origin, drop = FALSE] +
        log_to[, destination, drop = FALSE] + log_transition
      matrix(colSums(active * exp(joint)), n_profiles, n_profiles)
    }))
  }
  list(posterior = posterior, initial = colSums(occupancy[[1L]]),
       transition = counts)
}

#' Check fixed transition inputs once before iterating
#' @param x,layout,parameters Inputs to the homogeneous transition E-step.
#' @noRd
.multilpa_transition_validate <- function(x, layout, parameters) {
  n_profiles <- nrow(parameters$means)
  n_types <- length(parameters$group_probabilities)
  stopifnot(
    "`x` must be a numeric matrix without infinite values" =
      is.matrix(x) && is.numeric(x) && !any(is.infinite(x)),
    "`layout` must place every observation in a sequence" =
      is.list(layout) && sum(!is.na(layout$slot)) == nrow(x) &&
      identical(dim(layout$slot), dim(layout$within)),
    "initial and transition probabilities must match the model dimensions" =
      identical(dim(parameters$initial_probabilities), c(n_types, n_profiles)) &&
      identical(dim(parameters$transition_probabilities),
                c(n_profiles, n_profiles, n_types)),
    "transition parameters must be finite with positive normalized probabilities" =
      .multilpa_valid_parameters(parameters))
}


#' Profile prevalence implied by a transition expectation
#'
#' The prevalence of each profile within a group class is a summary of the
#' fitted model rather than a free parameter of it, so it is read off the
#' expectation that is actually reported, exactly as the empty transition rows
#' are. EM leaves the parameters one M-step behind that expectation, and with
#' `max_iter = 0` there is no M-step at all, which previously left the field
#' absent and the assembled fit unusable.
#'
#' @param expectation Expected counts from `.multilpa_transition_expectation()`.
#' @return A group-classes-by-profiles numeric matrix whose rows sum to one.
#' @noRd
.multilpa_transition_prevalence <- function(expectation) {
  stopifnot(
    "`expectation` must carry one joint posterior matrix per group class" =
      is.list(expectation) && is.list(expectation$joint) &&
      length(expectation$joint) >= 1L &&
      all(vapply(expectation$joint, is.matrix, logical(1))))
  counts <- t(vapply(expectation$joint, colSums,
                     numeric(ncol(expectation$joint[[1L]]))))
  totals <- rowSums(counts)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop(errorCondition(
      paste("A group class has no effective membership, so the profile",
            "prevalence within it is not defined."),
      class = "latents_empty_profile", call = NULL))
  }
  counts / totals
}

#' Nested expectation step for a latent transition model
#' @param x Observation-by-indicator matrix of continuous indicators.
#' @param layout Sequence layout.
#' @param parameters Model parameters, including `initial_probabilities` and
#'   `transition_probabilities`.
#' @param codes Optional integer matrix of categorical indicator codes.
#' @param validate Check inputs for direct calls; disabled inside validated fits.
#' @return Log likelihood, posteriors and expected sequence counts.
#' @noRd
.multilpa_transition_expectation <- function(x, layout, parameters, codes = NULL,
                                              validate = TRUE) {
  if (validate) .multilpa_transition_validate(x, layout, parameters)
  n_profiles <- nrow(parameters$means)
  n_types <- length(parameters$group_probabilities)
  n_groups <- nrow(layout$slot)
  measurement <- .latents_measurement_log_density(x, parameters, codes, offset = FALSE)
  log_density <- measurement$log_density
  emission <- .multilpa_sequence_emission(log_density, layout, validate = FALSE)
  passes <- lapply(seq_len(n_types), function(type) {
    .multilpa_forward_backward(emission$log_density, layout,
                               parameters$initial_probabilities[type, ],
                               matrix(parameters$transition_probabilities[, , type],
                                      n_profiles, n_profiles), validate = FALSE)
  })
  group_scores <- matrix(vapply(passes, `[[`, numeric(n_groups), "log_scaled"),
                         n_groups, n_types)
  weighted <- sweep(group_scores, 2L, log(parameters$group_probabilities), "+")
  group_log_likelihood <- .multilpa_log_sum_exp(weighted) + emission$offset
  group_posteriors <- exp(sweep(weighted, 1L, .multilpa_row_max(weighted), "-"))
  group_posteriors <- group_posteriors / rowSums(group_posteriors)
  sequence <- lapply(seq_len(n_types), function(type) {
    .multilpa_sequence_moments(passes[[type]], emission$log_density, layout,
                               matrix(parameters$transition_probabilities[, , type],
                                      n_profiles, n_profiles),
                               group_posteriors[, type], nrow(x), validate = FALSE)
  })
  joint <- lapply(sequence, `[[`, "posterior")
  subject_posteriors <- Reduce(`+`, joint)
  log_likelihood <- sum(group_log_likelihood)
  # The expected counts are checked beside the likelihood rather than left to
  # the M-step. A likelihood formed from the scaled forward pass alone stays
  # finite when the sequence moments do not, so a fit could otherwise converge
  # and still report NaN transition counts.
  sequence_counts <- unlist(lapply(sequence, function(moments) {
    c(moments$initial, moments$transition)
  }), use.names = FALSE)
  if (!is.finite(log_likelihood) || any(!is.finite(subject_posteriors)) ||
      any(!is.finite(sequence_counts)) || any(sequence_counts < 0)) {
    stop("Non-finite likelihood, posterior probabilities or sequence counts.")
  }
  list(log_likelihood = log_likelihood,
       group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors,
       subject_posteriors = subject_posteriors, joint = joint,
       sequence = sequence, gaussian_moments = measurement$moments)
}

#' Maximize the expected complete-data log likelihood of a transition model
#'
#' The measurement parameters are those of the cross-sectional model and are
#' maximized by the same function; only the initial and transition
#' probabilities are new, and both are bounded multinomial shares.
#'
#' @param x Observation-by-indicator matrix.
#' @param expectation Expected counts from `.multilpa_transition_expectation()`.
#' @param variance_model Either varying or equal across profiles.
#' @param min_variance Lower bound on every indicator variance.
#' @param covariance_model Diagonal or full residual covariance.
#' @param codes Optional integer matrix of categorical indicator codes.
#' @param n_categories Integer vector of category counts per categorical indicator.
#' @param min_probability Lower bound on every fitted probability.
#' @return Updated model parameters.
#' @noRd
.multilpa_transition_maximization <- function(x, expectation, variance_model,
                                              min_variance,
                                              covariance_model = "diagonal",
                                              codes = NULL, n_categories = NULL,
                                              min_probability = 1e-10,
                                              structure = NULL, previous = NULL) {
  parameters <- .multilpa_maximization(x, expectation, variance_model,
                                       min_variance, covariance_model, codes,
                                       n_categories, min_probability,
                                       structure = structure, previous = previous)
  n_profiles <- nrow(parameters$means)
  n_types <- length(expectation$sequence)
  bound <- function(counts) {
    .multilpa_bound_probabilities(counts, min_probability = min_probability)
  }
  initial <- t(matrix(vapply(expectation$sequence, function(moments) {
    bound(moments$initial)
  }, numeric(n_profiles)), n_profiles, n_types))
  transition <- array(vapply(expectation$sequence, function(moments) {
    t(apply(moments$transition, 1L, bound))
  }, numeric(n_profiles * n_profiles)),
  c(n_profiles, n_profiles, n_types))
  # A profile that no group occupies before its final occasion leaves its
  # transition row with no information. The bounded solution is then uniform,
  # which is a choice rather than an estimate, so it is reported rather than
  # passed off as fitted.
  empty <- vapply(expectation$sequence, function(moments) {
    rowSums(moments$transition) <= 0
  }, logical(n_profiles))
  # The prevalence of each profile within a group class is a summary of the
  # expectation rather than a free parameter of the model. It is carried here
  # so that a single M-step is self-describing, but the fitted object reports
  # the one `.multilpa_transition_prevalence()` reads off the final
  # expectation, which the last M-step never saw -- and which, when `max_iter`
  # is zero, is the only one there is.
  parameters <- c(parameters[setdiff(names(parameters), "profile_probabilities")],
                  list(initial_probabilities = initial,
                       transition_probabilities = transition,
                       profile_prevalence = parameters$profile_probabilities,
                       empty_transition_rows = matrix(empty, n_profiles, n_types)))
  if (any(!is.finite(unlist(parameters[c("initial_probabilities",
                                         "transition_probabilities")],
                            use.names = FALSE)))) {
    stop("Non-finite transition parameters in the M-step.")
  }
  parameters
}

#' Run EM for a latent transition model
#' @param x Observation-by-indicator matrix.
#' @param layout Sequence layout.
#' @param parameters Initial parameters.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param max_iter Maximum number of updates.
#' @param tol Relative log-likelihood tolerance.
#' @param covariance_model Diagonal or full residual covariance.
#' @param codes Optional categorical codes.
#' @param n_categories Category counts per categorical indicator.
#' @param min_probability Lower bound on every fitted probability.
#' @param acceleration SQUAREM or plain EM.
#' @param pause,resume As for `.latents_em()`.
#' @return Parameters, the final expectation, convergence and history; with
#'   `pause` or `resume`, also the `checkpoint`.
#' @noRd
.multilpa_transition_em <- function(x, layout, parameters, variance_model,
                                    min_variance, max_iter, tol,
                                    covariance_model = "diagonal", codes = NULL,
                                    n_categories = NULL, min_probability = 1e-10,
                                    structure = NULL, pause = Inf, resume = NULL,
                                    acceleration = "squarem") {
  stopifnot(is.matrix(x), is.list(parameters), max_iter >= 0L, tol > 0,
            min_variance > 0, variance_model %in% c("varying", "equal"))
  .multilpa_transition_validate(x, layout, parameters)
  accelerator <- .lta_accelerator(acceleration)
  fit <- .latents_em(
    parameters,
    evaluate = function(point) {
      .multilpa_transition_expectation(x, layout, point, codes, validate = FALSE)
    },
    maximize = function(point, expectation) {
      .multilpa_transition_maximization(
        x, expectation, variance_model, min_variance, covariance_model, codes,
        n_categories, min_probability, structure, point)
    },
    max_iter = max_iter, tol = tol, pause = pause, resume = resume,
    settings = .latents_em_settings(accelerate = accelerator))
  result <- list(parameters = fit$state, expectation = fit$expectation,
                 converged = fit$converged, iterations = fit$iterations,
                 history = fit$history)
  if (!is.null(fit$checkpoint)) result$checkpoint <- fit$checkpoint
  result
}

#' Initialize a latent transition model
#'
#' The measurement model and the group-class split reuse the cross-sectional
#' initializer. Transition matrices start persistent, because a mixture whose
#' states are freely exchangeable at every occasion carries no sequence
#' information to separate them; the persistence is randomized across starts.
#'
#' @param x Observation-by-indicator matrix.
#' @param group_index Integer group indices.
#' @param n_profiles Number of profiles.
#' @param n_types Number of group classes.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param start_index Restart number.
#' @param covariance_model Diagonal or full residual covariance.
#' @param codes Optional categorical codes.
#' @param n_categories Category counts per categorical indicator.
#' @param min_probability Lower bound on every fitted probability.
#' @param kind The kind of start, as for `.multilpa_initialize()`.
#' @param row_weights Optional sampling weights in observation order, as for
#'   `.multilpa_initialize()`.
#' @return Strictly positive initial parameters for EM.
#' @noRd
.multilpa_transition_initialize <- function(x, group_index, n_profiles, n_types,
                                            variance_model, min_variance,
                                            start_index,
                                            covariance_model = "diagonal",
                                            codes = NULL, n_categories = NULL,
                                            min_probability = 1e-10, extra = NULL,
                                            kind = "kmeans", row_weights = NULL) {
  # The first start is Ward's hierarchical clustering, as in multilpa().
  parameters <- .multilpa_initialize(x, group_index, n_profiles, n_types,
                                     variance_model, min_variance, start_index,
                                     covariance_model, codes, n_categories,
                                     min_probability,
                                     hierarchical = start_index == 1L, extra = extra,
                                     kind = kind, row_weights = row_weights)
  initial <- parameters$profile_probabilities
  persistence <- if (start_index == 1L) 0.7 else stats::runif(1L, 0.4, 0.9)
  transition <- array(vapply(seq_len(n_types), function(type) {
    row <- matrix(initial[type, ], n_profiles, n_profiles, byrow = TRUE)
    bounded <- persistence * diag(n_profiles) + (1 - persistence) * row
    t(apply(bounded, 1L, .multilpa_bound_probabilities,
            min_probability = min_probability))
  }, numeric(n_profiles * n_profiles)), c(n_profiles, n_profiles, n_types))
  parameters$initial_probabilities <- t(matrix(
    vapply(seq_len(n_types), function(type) {
      .multilpa_bound_probabilities(initial[type, ], min_probability)
    }, numeric(n_profiles)), n_profiles, n_types))
  parameters$transition_probabilities <- transition
  parameters$profile_probabilities <- NULL
  parameters
}

#' Refuse numbers written as text whose alphabetical order is not numeric
#'
#' Occasions are ordered by `sort()`, which orders character strings
#' alphabetically. Positions written as `"1"`, ..., `"10"` would put `"10"`
#' second and fit transitions between the wrong occasions without a word.
#' Strings whose alphabetical order is already right (zero-padded numbers,
#' ISO dates, labels) are accepted.
#'
#' @param values The time values.
#' @param time The time column name, for the message.
#' @return `NULL`, invisibly; raises `latents_bad_time` (also `latents_bad_data`).
#' @noRd
.lta_check_time_order <- function(values, time) {
  if (!is.character(values)) return(invisible(NULL))
  # type.convert() turns the strings into numbers only when all of them parse.
  numbers <- utils::type.convert(values, as.is = TRUE)
  if (!is.numeric(numbers) || !all(is.finite(numbers))) return(invisible(NULL))
  lexical <- sort(unique(values))
  if (!is.unsorted(numbers[match(lexical, values)])) return(invisible(NULL))
  stop(errorCondition(sprintf(paste(
    "`%s` holds numbers stored as text, which sort alphabetically (\"10\"",
    "before \"2\") and would put the occasions out of order. Convert it with",
    "as.numeric() first."), time),
    class = c("latents_bad_time", "latents_bad_data"), call = NULL))
}

#' Count the free parameters of a latent transition model
#' @param n_profiles Number of profiles.
#' @param n_group_classes Number of group classes.
#' @param n_continuous Number of continuous indicators.
#' @param n_categories Category counts per categorical indicator, or `NULL`.
#' @param variance_model Variance constraint.
#' @param covariance_model Diagonal or full residual covariance.
#' @return Number of free parameters.
#' @noRd
.multilpa_count_transition_parameters <- function(n_profiles, n_group_classes,
                                                  n_continuous, n_categories,
                                                  variance_model,
                                                  covariance_model,
                                                  structure = NULL) {
  # The cross-sectional count carries the measurement model and the group-class
  # split. Its per-class profile prevalences are replaced here by a per-class
  # initial distribution and a per-class transition matrix.
  cross_sectional <- .multilpa_count_parameters(
    n_profiles, n_group_classes, n_continuous, n_categories, variance_model,
    covariance_model, structure)
  cross_sectional + n_group_classes * n_profiles * (n_profiles - 1L)
}

#' Fit a latent transition model
#'
#' Fits latent profiles to repeated observations of the same group, estimates
#' the probabilities of moving between profiles from one occasion to the next,
#' and returns a `multilpa_transitions` (or `multilpa_lta`) object. The
#' measurement model is that of [multilpa()] (Gaussian, categorical or mixed
#' indicators) and is invariant across occasions by default, so a profile has
#' the same meaning at every occasion. With more than one group class, each
#' class has its own initial distribution and its own transition matrix, so
#' the group classes are a mixture over transition patterns. [multilpa()] is
#' the cross-sectional counterpart.
#'
#' @inheritParams multilpa
#' @param search How the starts are found: `"kmeans"` (the default) or
#'   `"emEM"`, as described for [multilpa()]. Unlike [multilpa()], `lta()`
#'   does not choose its default by the covariance structure.
#' @param cores Number of processor cores for the random starts (default 1).
#'   The starting values are drawn in order in the current process and only
#'   the EM runs are distributed, so with `seed` given the result is identical
#'   to `cores = 1`. Under `search = "emEM"` both screening and the final runs
#'   are distributed. Workers are forked processes on macOS and Linux, and a
#'   socket cluster on Windows and on macOS with Apple's Accelerate BLAS, as
#'   in [enumerate_classes()]. Lost worker results are retried once with half
#'   as many workers (rounded up), with one `latents_worker_lost` warning
#'   reporting losses and recovery and suggesting fewer cores. A start still
#'   lost is recorded as failed; the fit uses the remaining starts. If all
#'   starts fail, the usual `latents_all_starts_failed` error is raised.
#'   Socket workers use bounded chunks; completed tasks within a failed chunk
#'   may also be retried. To fit many models, give `cores` to
#'   [enumerate_classes()], which runs starts across candidates in one pool.
#' @param id Name of the observed group identifier column. Character, factor
#'   and numeric identifiers are supported; missing identifiers are not.
#'   Required; `id = NULL` is not accepted.
#' @param n_profiles Integer number of individual profiles, at least two
#'   (`latents_bad_transition` otherwise).
#' @param time Name of the column giving each observation's occasion within
#'   its group. Required. Values must be complete and unique within each
#'   group. They are ordered as [sort()] orders them: numbers numerically, a
#'   factor by its levels, and character strings alphabetically. Numbers
#'   stored as text whose alphabetical order is not their numeric order
#'   (`"2"`, `"10"`) raise `latents_bad_time` (also `latents_bad_data`);
#'   zero-padded numbers and ISO dates are accepted.
#' @param n_group_classes Positive integer number of latent group classes
#'   (default 1). One gives an ordinary latent transition model; more than one
#'   fits a mixture of transition models over groups.
#' @param select_start Which start the fit reports. `"likelihood"` (the
#'   default) takes the highest log likelihood across starts, preferring a
#'   converged start only among those within a relative `1e-10` of it.
#'   `"converged"` takes the highest log likelihood among the converged starts
#'   when at least one converged, and among all starts otherwise. The two
#'   differ when an unconverged start exceeds a converged one by a negligible
#'   amount, as happens when a membership logit drifts along a flat ridge of
#'   the likelihood; only the converged start is a maximum on which standard
#'   errors can be based. `summary()` and `get_results(fit, "starts")` show
#'   every start.
#' @param occasions How a group's occasions are placed on the transition grid.
#'   `"observed"` (the default) numbers each group's own observations
#'   consecutively, so a transition always links two adjacent observations.
#'   `"grid"` places them on the grid of every position seen in the data, so a
#'   group that skips a position still makes a transition across the gap and
#'   contributes no measurement information at it. The two agree when every
#'   group is observed at every position. Under either, a group's first
#'   observation is its first occasion: occasion-specific transitions and
#'   measurement count occasions from each group's own start, not from the
#'   first position in the data.
#' @param transitions `"homogeneous"` (the default) uses one transition matrix
#'   for every occasion; `"occasion"` estimates a separate one for each move
#'   (occasion 1 to 2, 2 to 3, ...).
#' @param transition_covariates Names of columns whose values shift the
#'   transition probabilities through a multinomial logit per origin profile:
#'   the log odds of moving to each other profile rather than staying. A
#'   column may change over occasions; the value at the destination occasion
#'   is used. Factors are expanded to indicator columns. With `occasions =
#'   "grid"`, a move into a position the group skipped has no row to read the
#'   covariates from, and the fit raises `latents_bad_data`; use
#'   `occasions = "observed"`, or add the skipped rows with their covariate
#'   values (their indicators may be missing under `missing = "fiml"`).
#' @param initial_covariates Names of columns that shift the initial profile
#'   distribution (a multinomial logit, reference: the last profile), read at
#'   each group's first occasion.
#' @param measurement `"invariant"` (the default) keeps the profiles' means,
#'   variances and response probabilities equal across occasions;
#'   `"occasion"` estimates them separately at each occasion.
#' @param tol Positive relative log-likelihood tolerance (default `1e-10`,
#'   tighter than the `1e-8` of [multilpa()]). Convergence requires an
#'   absolute change no greater than `tol * (1 + abs(previous log likelihood))`.
#'   The group-class shares are re-estimated from posteriors one EM step
#'   behind them; at `1e-8` the two differed by up to about `2e-4`, at `1e-10`
#'   by about `2e-5`, for roughly 1.5 times the iterations.
#' @param acceleration `"squarem"` (the default) accelerates EM with SQUAREM
#'   (Varadhan & Roland, 2008), as [multilpa()] does: each cycle extrapolates
#'   along the last two EM steps, including the initial and transition
#'   probabilities, and ends with an ordinary EM step. An extrapolated point
#'   that would lower the likelihood, or that leaves a probability
#'   non-positive or a row not summing to one, is replaced by the plain EM
#'   point, so the iteration stays monotone. Convergence is assessed over a
#'   whole cycle, so it is never looser than plain EM's. Under the strict
#'   default `tol` plain EM can creep along a flat likelihood for thousands of
#'   steps; SQUAREM usually stops in far fewer, at the same or a slightly
#'   higher maximum. `"none"` runs plain EM, which reproduces fits from
#'   versions before this argument existed. `iterations` and `max_iter` count
#'   EM steps under either setting.
#' @param max_iter Nonnegative integer maximum number of EM updates per start
#'   (default 5000, higher than the 1000 of [multilpa()] to match the tighter
#'   `tol`). A fit that stops short warns with `latents_unconverged`.
#'   `max_iter = 0` updates nothing and evaluates the starting values.
#' @param order `1` (the default) or `2`. With `2`, the profile at each
#'   occasion from the third on depends on the profiles at the two previous
#'   occasions (a second-order Markov chain); the first move stays first
#'   order.
#' @param mover_stayer `FALSE` (the default), or `TRUE` to add a class of
#'   stayers: groups that never change profile (identity transitions) and have
#'   their own initial profile distribution (Goodman's mover-stayer model).
#'   `n_group_classes` then counts the mover classes; the stayer class is
#'   reported as `"stayers"`.
#' @param weights `NULL`, or the name of a column of `data` holding one
#'   sampling weight per `id` (constant over its occasions). As in
#'   [multilpa()]: pseudo maximum likelihood with the weights scaled to sum to
#'   the number of `id` units, sandwich standard errors by default, and
#'   integer weights equal to repeating each sequence. A weighted fit uses the
#'   general transition engine and returns a `multilpa_lta`, whose
#'   model-implied `"transitions"` and `"initial"` tables average the groups
#'   with their weights.
#' @param ordinal,count,count_model,count_dispersion Ordinal and count
#'   indicators, as in [multilpa()]; with `measurement = "occasion"` their
#'   parameters are estimated per occasion. Read them with
#'   `get_results(fit, "ordinal")` and `get_results(fit, "count_means")`. A
#'   fit with them uses the general transition engine.
#' @param model A covariance structure as `multilpa()` takes it (`"EEI"`,
#'   `"VVI"`, `"VEI"`, ... ; mclust's codes). Standard errors are given for
#'   every structure with the extensions (covariate, occasion-varying or
#'   second-order transitions, occasion-specific measurement), and for EEI,
#'   VVI, EEE and VVV in the homogeneous model, which refuses the other ten
#'   with `latents_unsupported_inference`.
#' @return A `multilpa_transitions` object containing `means`, `variances`,
#'   optional `covariances` and `response_probabilities`,
#'   `initial_probabilities` (group classes by profiles),
#'   `transition_probabilities` (profiles by profiles by group classes, rows
#'   indexing the profile moved from), `group_probabilities`, posterior
#'   matrices, classifications, log likelihood, information criteria, restart
#'   diagnostics and convergence history. [get_results()] returns the tidy
#'   transition table and the other tables, [parameter_inference()] Wald or
#'   bootstrap standard errors, and [bootstrap_lrt()] compares nested fits.
#'   Global optimality is not guaranteed.
#'
#'   With `transitions = "occasion"`, `transition_covariates`,
#'   `initial_covariates`, `measurement = "occasion"`, `order = 2`,
#'   `mover_stayer = TRUE`, `weights`, `ordinal` or `count`, the result is a
#'   `multilpa_lta` object, read with [get_results()]: transition and initial
#'   logit coefficients with Wald standard errors (observed, robust or OPG,
#'   from analytic scores), model-implied transition probabilities per
#'   occasion, and profiles per occasion. Missing indicators
#'   (`missing = "fiml"`) and covariance structures are supported; diagonal
#'   fits are finished by a quasi-Newton search on the exact likelihood,
#'   others by EM. Agreement with Mplus (User's Guide examples 8.13, 8.14;
#'   occasion-specific thresholds), depmixS4 and LMest is recorded in
#'   `equivalence/lta-extensions/` of the source repository.
#' @details By default transitions are first order and homogeneous over
#'   occasions: the probability of moving from one profile to another does not
#'   depend on the occasion or on earlier profiles. Measurement parameters are
#'   shared across occasions and across group classes. A profile that no group
#'   occupies before its final occasion leaves its transition row without
#'   information; the row is then uniform by construction, not estimated. Such
#'   a row raises a warning when the model is fitted, is flagged by the
#'   `estimated` column of the transition table and is listed by
#'   `get_results(fit, "transitions", estimated = FALSE)`. Profile labels are
#'   arbitrary and are not comparable across fits without alignment.
#'
#'   `max_iter = 0` updates nothing: the starting values are evaluated and
#'   returned with the expectation they imply, and the fit reports
#'   `converged = FALSE` after zero iterations. Every reported quantity,
#'   including `expected_count` and the implied profile prevalence, is then
#'   read off that single expectation.
#' @references Collins, L. M., & Lanza, S. T. (2010). Latent class and latent
#'   transition analysis. Wiley.
#'
#'   Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods
#'   for accelerating the convergence of any EM algorithm. Scandinavian
#'   Journal of Statistics, 35, 335--353.
#'
#'   Vermunt, J. K. (2003). Multilevel latent class models. Sociological
#'   Methodology, 33, 213--239. doi:10.1111/j.0081-1750.2003.t01-1-00131.x.
#' @seealso [get_results()] for the fitted transition probabilities and for the
#'   assignments in order, and [multilpa()] for the
#'   cross-sectional model this one shares its measurement parameters with.
#' @examples
#' # Students take their courses in their own order, which `sequence` records,
#' # so the same enrolments that fit a cross-sectional model fit a transition
#' # one. Engagement mostly persists from one course to the next.
#' fit <- lta(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, time = "sequence", n_starts = 2, seed = 1
#' )
#' get_results(fit, what = "transitions")
#' summary(fit)
#' @export
lta <- function(data, vars, id, n_profiles, time,
                            n_group_classes = 1L,
                            variance_model = c("varying", "equal"),
                            n_starts = 10L, max_iter = 5000L, tol = 1e-10,
                            min_variance = 1e-6, seed = NULL,
                            missing = c("error", "fiml"),
                            covariance_model = c("diagonal", "full"),
                            categorical = character(), min_probability = 1e-10,
                            occasions = c("observed", "grid"),
                            select_start = c("likelihood", "converged"),
                            transitions = c("homogeneous", "occasion"),
                            transition_covariates = character(),
                            initial_covariates = character(),
                            measurement = c("invariant", "occasion"),
                            order = 1L, model = NULL, mover_stayer = FALSE,
                            weights = NULL, ordinal = character(),
                            count = character(),
                            count_model = c("poisson", "negative_binomial"),
                            count_dispersion = c("varying", "equal"),
                            drop_empty = TRUE, search = c("kmeans", "emEM"),
                            n_screen = NULL,
                            acceleration = c("squarem", "none"), cores = 1L) {
  # Setup forces the public arguments and draws every start in its seed frame.
  plan <- do.call(.lta_setup, c(as.list(environment()), list(call = match.call())),
                  quote = TRUE)
  .lta_run_plan(plan, cores)
}

#' Prepare a latent transition fit without running its starts
#' @inheritParams lta
#' @param call The public fitting call to retain on the result.
#' @return A prepared start plan with a runner and a parent-side finaliser.
#' @noRd
.lta_setup <- function(data, vars, id, n_profiles, time,
                            n_group_classes = 1L,
                            variance_model = c("varying", "equal"),
                            n_starts = 10L, max_iter = 5000L, tol = 1e-10,
                            min_variance = 1e-6, seed = NULL,
                            missing = c("error", "fiml"),
                            covariance_model = c("diagonal", "full"),
                            categorical = character(), min_probability = 1e-10,
                            occasions = c("observed", "grid"),
                            select_start = c("likelihood", "converged"),
                            transitions = c("homogeneous", "occasion"),
                            transition_covariates = character(),
                            initial_covariates = character(),
                            measurement = c("invariant", "occasion"),
                            order = 1L, model = NULL, mover_stayer = FALSE,
                            weights = NULL, ordinal = character(),
                            count = character(),
                            count_model = c("poisson", "negative_binomial"),
                            count_dispersion = c("varying", "equal"),
                            drop_empty = TRUE, search = c("kmeans", "emEM"),
                            n_screen = NULL,
                            acceleration = c("squarem", "none"), cores = 1L, call = NULL) {
  cores <- .latents_check_cores(cores)
  acceleration <- match.arg(acceleration)
  count_model <- match.arg(count_model)
  count_dispersion <- match.arg(count_dispersion)
  .multilpa_check_drop_empty(drop_empty)
  select_start <- match.arg(select_start)
  search <- match.arg(search)
  transitions <- match.arg(transitions)
  measurement_model <- match.arg(measurement)
  if (!is.numeric(order) || length(order) != 1L || !order %in% c(1, 2)) {
    stop(errorCondition("`order` must be 1 or 2.", class = "latents_bad_argument",
                        call = NULL))
  }
  stopifnot("`mover_stayer` must be TRUE or FALSE" = isTRUE(mover_stayer) ||
              isFALSE(mover_stayer))
  general <- mover_stayer || order == 2 || !identical(transitions, "homogeneous") ||
    length(transition_covariates) > 0L || length(initial_covariates) > 0L ||
    !identical(measurement_model, "invariant") || !is.null(weights) ||
    length(ordinal) > 0L || length(count) > 0L
  .latents_check_extra_arguments(vars, categorical, ordinal, count)
  stopifnot(is.data.frame(data), is.character(vars), is.character(id),
            "`categorical` must be a character vector of indicator names" =
              is.character(categorical) && !anyNA(categorical),
            "`min_probability` must be a single number in (0, 1)" =
              is.numeric(min_probability) && length(min_probability) == 1L &&
              is.finite(min_probability) && min_probability > 0 &&
              min_probability < 1,
            "`time` must be supplied; a transition model is defined by the ordering" =
              !is.null(time))
  if (is.null(call)) call <- match.call()
  variance_model <- match.arg(variance_model)
  covariance_model <- match.arg(covariance_model)
  missing <- match.arg(missing)
  occasions <- match.arg(occasions)
  # `model` names the covariance structure in mclust's three letters, as in
  # multilpa(); NULL keeps the `variance_model` x `covariance_model` choice.
  structure <- NULL
  if (!is.null(model)) {
    if (!is.character(model) || length(model) != 1L ||
        !model %in% .multilpa_structures()) {
      stop(errorCondition(sprintf("`model` must be one of %s.",
        paste(sprintf("\"%s\"", .multilpa_structures()), collapse = ", ")),
        class = "latents_bad_argument", call = NULL))
    }
    pieces <- .multilpa_structure_arguments(model)
    structure <- .multilpa_resolve_structure(variance_model, covariance_model,
                                             pieces$volume, pieces$shape,
                                             pieces$orientation)
    if (.multilpa_is_ellipsoidal(structure)) covariance_model <- "full"
    if (structure %in% c("EEI", "EEE")) variance_model <- "equal"
    if (structure %in% c("VVI", "VVV")) variance_model <- "varying"
    .multilpa_reset_structure_log()
  }
  # The data contract is checked BEFORE the ordering column. Reversed, a missing
  # group column surfaced as an internal `split()` failure ("group length is 0
  # but data length > 0") instead of this package's own classed error.
  .multilpa_check_arguments(data, vars, id, n_profiles, n_group_classes,
                            n_starts, max_iter, tol, min_variance,
                            min_probability, seed, categorical)
  time_values <- .multilpa_time_values(data, time, id, vars)
  if (general) {
    # The general engine: transitions per occasion or covariate, covariate
    # initial distribution, or measurement per occasion.
    if (n_profiles < 2L) {
      stop(errorCondition(
        "`n_profiles` must be at least two; a single profile has nothing to move between.",
        class = "latents_bad_transition", call = NULL))
    }
    return(.lta_setup_general(data, vars, id, time, n_profiles, n_group_classes,
                            variance_model, n_starts, max_iter, tol, min_variance,
                            seed, categorical, min_probability, occasions,
                            transitions, transition_covariates, initial_covariates,
                            measurement_model, select_start, call, as.integer(order),
                            mover_stayer, missing, covariance_model, structure,
                            weights, ordinal, count, count_model, count_dispersion,
                            drop_empty = drop_empty, search = search,
                            n_screen = n_screen, acceleration = acceleration,
                            cores = cores))
  }
  if (n_profiles < 2L) {
    stop(errorCondition(
      "`n_profiles` must be at least two; a single profile has nothing to move between.",
      class = "latents_bad_transition", call = NULL))
  }
  measurement <- .multilpa_prepare_indicators(data, vars, categorical,
                                              missing, min_probability)
  continuous <- measurement$continuous
  codes <- measurement$codes
  n_categories <- measurement$n_categories
  x <- measurement$x
  groups <- .multilpa_prepare_groups(data[[id]])
  layout <- .multilpa_sequence_layout(groups$index, time_values, groups$n, occasions)
  if (layout$n_occasions < 2L) {
    stop(errorCondition(
      "No group is observed at two occasions, so no transition can be estimated.",
      class = "latents_bad_transition", call = NULL))
  }
  distinct_rows <- if (is.null(codes)) nrow(unique(x)) else
    nrow(unique(cbind(x, codes)))
  if (n_profiles > distinct_rows) {
    stop(errorCondition("n_profiles cannot exceed the number of distinct observed indicator rows.",
        class = "latents_unidentified", call = NULL))
  }
  if (n_group_classes > groups$n) {
    stop(errorCondition("n_group_classes cannot exceed the number of groups.",
        class = "latents_unidentified", call = NULL))
  }
  if (n_group_classes > 1L && all(groups$sizes == 1L)) {
    stop(errorCondition("Multiple group classes are not identifiable with only singleton groups.",
        class = "latents_unidentified", call = NULL))
  }
  .latents_local_seed(seed)
  # Centering is a translation of the indicators, with no density Jacobian and
  # no change of scale; the fitted means are shifted back before they are read.
  centers <- if (ncol(x) > 0L) colMeans(x, na.rm = TRUE) else numeric(0)
  x <- sweep(x, 2L, centers, "-")
  if (any(!is.finite(x[!is.na(x)]^2))) {
    stop(errorCondition("Indicator scales overflow squared residuals; rescale the data.",
        class = "latents_bad_data", call = NULL))
  }
  # Each start runs to completion or is recorded as failed; the best is
  # chosen by the shared rule (see .latents_run_starts()).
  # Every start is prepared before any is run; EM draws no random numbers, so
  # the draws are those of preparing and running each start in turn.
  screen <- .latents_search_screen(search, n_starts, n_screen, max_iter)
  prepare <- function(start_index) {
    initial <- .multilpa_transition_initialize(
      x, groups$index, n_profiles, n_group_classes, variance_model,
      min_variance, start_index, covariance_model, codes, n_categories,
      min_probability, kind = .latents_start_kind(start_index, n_starts))
    if (!is.null(structure)) {
      initial <- .multilpa_project_start(initial, structure, nrow(x), min_variance)
    }
    initial
  }
  initial <- .latents_prepare_starts(
    if (is.null(screen)) n_starts else screen$n_pool, prepare)
  context <- as.list(environment())
  context[c("prepare", "initial")] <- NULL
  .lta_plan(context, initial, .lta_run_start, .lta_finalise,
             c("x", "layout", "variance_model", "min_variance", "max_iter", "tol",
               "covariance_model", "codes", "n_categories", "min_probability",
               "structure", "acceleration"))
}

#' Run one homogeneous transition start
#' @param start_index Original start index.
#' @param initial Starting values, or NULL when resuming.
#' @param pause,resume Screening budget and checkpoint.
#' @param context Forced EM inputs from setup.
#' @return One EM attempt.
#' @noRd
.lta_run_start <- function(start_index, initial, pause = Inf, resume = NULL, context) {
  with(context, {
    .multilpa_transition_em(x, layout, initial %||% resume$state, variance_model,
                            min_variance, max_iter, tol, covariance_model, codes,
                            n_categories, min_probability, structure,
                            pause = pause, resume = resume,
                            acceleration = acceleration)
  })
}

#' Finalise a homogeneous transition fit from its selected starts
#' @param starts The shared start-selection result.
#' @param context Forced setup inputs.
#' @return A fitted multilpa_transitions object.
#' @noRd
.lta_finalise <- function(starts, context) {
  with(context, {
    found <- starts
    attempts <- starts$attempts
    valid <- starts$valid
    scores <- starts$scores
    best_start <- starts$best_start
    best <- attempts[[best_start]]
    parameters <- best$parameters
    parameters$means <- sweep(parameters$means, 2L, centers, "+")
    profile_names <- paste0("profile_", seq_len(n_profiles))
    type_names <- paste0("group_class_", seq_len(n_group_classes))
    parameters <- .multilpa_label_parameters(parameters, profile_names, continuous,
                                             categorical, measurement$encoded$levels)
    dimnames(parameters$initial_probabilities) <- list(type_names, profile_names)
    # Taken from the expectation that is reported rather than from the M-step,
    # which saw the previous one -- and which never runs at all when `max_iter`
    # is zero, where reading the M-step's field left the assembly setting
    # dimnames on NULL.
    parameters$profile_prevalence <- matrix(
      .multilpa_transition_prevalence(best$expectation),
      n_group_classes, n_profiles, dimnames = list(type_names, profile_names))
    dimnames(parameters$transition_probabilities) <-
      list(profile_names, profile_names, type_names)
    names(parameters$group_probabilities) <- type_names
    subject_posteriors <- best$expectation$subject_posteriors
    group_posteriors <- best$expectation$group_posteriors
    dimnames(subject_posteriors) <- list(rownames(data), profile_names)
    dimnames(group_posteriors) <- list(groups$ids, type_names)
    n_parameters <- .multilpa_count_transition_parameters(
      n_profiles, n_group_classes, ncol(x), n_categories, variance_model,
      covariance_model, structure)
    log_likelihood <- best$expectation$log_likelihood
    transition_counts <- array(
      vapply(best$expectation$sequence, `[[`, numeric(n_profiles * n_profiles),
             "transition"),
      c(n_profiles, n_profiles, n_group_classes),
      dimnames = list(profile_names, profile_names, type_names))
    # EM leaves the parameters one M-step behind the final expectation, so the
    # empty rows are re-read from the counts that are actually reported rather
    # than from the ones the last update happened to see.
    parameters$empty_transition_rows <- matrix(
      apply(transition_counts, c(1L, 3L), sum) <= 0, n_profiles, n_group_classes,
      dimnames = list(profile_names, type_names))
    starts <- do.call(rbind, lapply(seq_along(attempts), function(start_index) {
      attempt <- attempts[[start_index]]
      data.frame(start = start_index, log_likelihood = scores[start_index],
                 converged = if (valid[start_index]) attempt$converged else FALSE,
                 iterations = if (valid[start_index]) attempt$iterations else NA_integer_,
                 error = if (valid[start_index]) NA_character_ else attempt$error)
    }))
    starts <- .latents_search_columns(starts, found, n_starts)
    observed_per_row <- rowSums(!is.na(x)) +
      if (is.null(codes)) 0L else rowSums(!is.na(codes))
    n_informative <- sum(observed_per_row > 0L)
    n_informative_groups <- length(unique(groups$index[observed_per_row > 0L]))
    # Under `drop_empty`, sizes and counts are of the classified rows and
    # groups, as Mplus counts them; otherwise of every row and group (0.9.37).
    counted_rows <- if (isTRUE(drop_empty)) observed_per_row > 0L else rep(TRUE, nrow(x))
    counted_groups <- if (isTRUE(drop_empty))
      seq_len(groups$n) %in% groups$index[observed_per_row > 0L] else rep(TRUE, groups$n)
    group_sample_size <- if (isTRUE(drop_empty)) n_informative_groups else groups$n
    boundary <- .multilpa_covariance_boundary(parameters, min_variance)
    empty_rows <- any(parameters$empty_transition_rows)
    small_classes <- any(colSums(subject_posteriors) < 1) ||
      any(colSums(group_posteriors) < 1)
    result <- c(parameters, list(
      call = call, vars = vars, continuous = continuous,
      continuous_types = vapply(data[continuous], typeof, character(1)),
      categorical = categorical,
      categorical_levels = measurement$encoded$levels,
      categorical_values = .multilpa_categorical_values(
        data, categorical, measurement$encoded),
      min_probability = min_probability,
      indicator_data = as.matrix(measurement$frame), categorical_data = codes,
      id = id, group_ids = groups$ids, time = time, time_values = time_values,
      group_values = groups$values, group_index = groups$index,
      group_sizes = setNames(groups$sizes, groups$ids),
      n_observations = nrow(x), n_informative = n_informative,
      informative = observed_per_row > 0L,
      n_informative_groups = n_informative_groups, drop_empty = drop_empty,
      n_groups = groups$n, n_profiles = as.integer(n_profiles),
      n_group_classes = as.integer(n_group_classes),
      n_occasions = layout$n_occasions, occasions = occasions,
      sequence_lengths = setNames(layout$span, groups$ids),
      balanced = all(layout$span == layout$n_occasions) &&
        !anyNA(layout$slot),
      variance_model = variance_model, covariance_model = covariance_model,
      covariance_structure = structure,
      missing = missing, min_variance = min_variance,
      standard_deviations = sqrt(parameters$variances),
      measurement_model = if (is.null(codes)) "gaussian" else
        if (ncol(x) == 0L) "categorical" else "mixed",
      transition_counts = transition_counts,
      subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
      subject_profiles = max.col(subject_posteriors, ties.method = "first"),
      group_classes = max.col(group_posteriors, ties.method = "first"),
      log_likelihood = log_likelihood,
      group_log_likelihood = setNames(best$expectation$group_log_likelihood,
                                      groups$ids),
      n_parameters = n_parameters, aic = -2 * log_likelihood + 2 * n_parameters,
      bic = -2 * log_likelihood + log(group_sample_size) * n_parameters,
      bic_groups = -2 * log_likelihood + log(group_sample_size) * n_parameters,
      bic_individual = -2 * log_likelihood + log(n_informative) * n_parameters,
      converged = best$converged, iterations = best$iterations,
      log_likelihood_history = best$history, starts = starts,
      best_start = best_start, n_failed_starts = sum(!valid), boundary = boundary,
      small_classes = small_classes,
      effective_profile_counts = colSums(subject_posteriors[counted_rows, , drop = FALSE]),
      effective_group_counts = colSums(group_posteriors[counted_groups, , drop = FALSE]),
      n_best_replicated = sum(valid & .latents_final_stage(found) &
                                abs(scores - log_likelihood) <=
                                1e-6 * (1 + abs(log_likelihood))),
      replication_tolerance = 1e-6 * (1 + abs(log_likelihood))),
      .latents_search_record(search, n_starts, screen, select_start, "allow"),
      list(max_iter = as.integer(max_iter), tol = tol, acceleration = acceleration))
    class(result) <- "multilpa_transitions"
    # Every qualification a caller might branch on carries its catalogued class,
    # so a simulation loop can muffle the ones it expects and let the rest through.
    if (any(!valid)) {
      warning(warningCondition(sprintf(
        paste("%d of %d starts failed; read their messages with",
              "get_results(fit, \"starts\")."),
        sum(!valid), length(attempts)),
        class = "latents_failed_starts", call = NULL))
    }
    if (!best$converged && max_iter > 0L) {
      warning(warningCondition(
        "The best start did not converge; increase max_iter and inspect starts.",
        class = "latents_unconverged", call = NULL))
    }
    if (boundary) {
      warning(warningCondition(if (covariance_model == "full")
        "A covariance eigenvalue reached min_variance; this is a bound-active constrained fit." else
        "A variance reached min_variance; this is a bound-active constrained fit.",
        class = "latents_boundary", call = NULL))
    }
    if (empty_rows) {
      warning(warningCondition(paste(
        "A profile is never occupied before a final occasion; its",
        "transition row is uniform by construction, not estimated.",
        "List the affected rows with",
        "get_results(fit, \"transitions\", estimated = FALSE)."),
        class = "latents_empty_transition_row", call = NULL))
    }
    if (small_classes) {
      warning(warningCondition(
        "A profile or group class has effective membership below one.",
        class = "latents_small_classes", call = NULL))
    }
    .latents_warn_thin_profiles(
      # Empty and categorical-only rows do not inform Gaussian measurement,
      # even when drop_empty = FALSE includes them in classification counts.
      colSums(subject_posteriors[rowSums(!is.na(x)) > 0L, , drop = FALSE]),
      .latents_profile_parameters(length(continuous), variance_model, covariance_model,
                                  structure))
    result$gaussian_diagnostics <- .latents_gaussian_diagnostics(
      parameters, x, result$subject_posteriors,
      min_variance = min_variance)
    result
  })
}

#' First-order transition probabilities between profiles
#'
#' Documented on `?get_results`, which is where a caller reaches this table from.
#'
#' @param x A fitted `multilpa_transitions` model.
#' @param estimated,stable `TRUE` or `FALSE` restricts the table; `NULL` keeps
#'   every row.
#' @return A base `data.frame`, one row per group class and ordered pair of
#'   profiles.
#' @noRd
.multilpa_transitions_table <- function(x, estimated = NULL, stable = NULL) {
  stopifnot(
    "`object` must be a fitted `multilpa_transitions` model" =
      inherits(x, "multilpa_transitions"),
    "`estimated` must be NULL, TRUE or FALSE" =
      .multilpa_optional_flag(estimated),
    "`stable` must be NULL, TRUE or FALSE" = .multilpa_optional_flag(stable))
  .multilpa_restrict_transitions(.multilpa_transition_frame(x),
                                 estimated = estimated, stable = stable)
}

#' Is a value absent, or a single non-missing logical?
#'
#' @param value The argument to check.
#' @return A single logical.
#' @noRd
.multilpa_optional_flag <- function(value) {
  is.null(value) ||
    (is.logical(value) && length(value) == 1L && !is.na(value))
}

#' Transition probabilities of a fit or of its summary
#'
#' Both a fitted model and its summary carry the arrays this table is built
#' from under the same names, so both get the same table from one definition.
#'
#' @param x A fitted `multilpa_transitions` model or its summary.
#' @return One row per group class and ordered pair of profiles.
#' @noRd
.multilpa_transition_frame <- function(x) {
  n_profiles <- x$n_profiles
  n_types <- x$n_group_classes
  cells <- n_profiles * n_profiles
  data.frame(
    group_class = rep(seq_len(n_types), each = cells),
    from = rep(rep(seq_len(n_profiles), each = n_profiles), times = n_types),
    to = rep(seq_len(n_profiles), times = n_profiles * n_types),
    probability = as.vector(aperm(x$transition_probabilities, c(2L, 1L, 3L))),
    expected_count = as.vector(aperm(x$transition_counts, c(2L, 1L, 3L))),
    stable = rep(as.vector(t(diag(TRUE, n_profiles))), times = n_types),
    estimated = rep(!as.vector(x$empty_transition_rows), each = n_profiles),
    group_class_probability = rep(unname(x$group_probabilities), each = cells),
    row.names = NULL)
}

#' Keep the transition rows a caller asked for
#'
#' @param frame A transition table, as built by `.multilpa_transition_frame()`.
#' @param estimated,stable `NULL` to keep every row, or a single logical to keep
#'   the rows whose column equals it.
#' @return The table, with row names renumbered.
#' @noRd
.multilpa_restrict_transitions <- function(frame, estimated, stable) {
  keep <- rep(TRUE, nrow(frame))
  if (!is.null(estimated)) keep <- keep & frame$estimated == estimated
  if (!is.null(stable)) keep <- keep & frame$stable == stable
  frame <- frame[keep, , drop = FALSE]
  row.names(frame) <- NULL
  frame
}

#' Initial profile probabilities as a tidy table
#' @param x A fitted `multilpa_transitions` model.
#' @return One row per group class and profile.
#' @noRd
.multilpa_initial_frame <- function(x) {
  n_profiles <- x$n_profiles
  n_types <- x$n_group_classes
  data.frame(
    group_class = rep(seq_len(n_types), each = n_profiles),
    profile = rep(seq_len(n_profiles), times = n_types),
    probability = as.vector(t(x$initial_probabilities)),
    prevalence = as.vector(t(x$profile_prevalence)),
    group_class_probability = rep(unname(x$group_probabilities),
                                  each = n_profiles),
    row.names = NULL)
}

#' @noRd
#' @export
as.data.frame.multilpa_transitions <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `multilpa_transitions`" = inherits(x, "multilpa_transitions"))
  .multilpa_coerce(x, row.names, list(...))
}

#' @noRd
#' @export
print.multilpa_transitions <- function(x, rows = 20L, ...) {
  stopifnot(inherits(x, "multilpa_transitions"))
  cat(sprintf("Latent transition model: %d profiles, %d group %s\n",
              x$n_profiles, x$n_group_classes,
              if (x$n_group_classes == 1L) "class" else "classes"))
  cat(sprintf("%d observations in %d groups, up to %d occasions (%s)\n",
              x$n_observations, x$n_groups, x$n_occasions,
              if (isTRUE(x$balanced)) "balanced" else
                sprintf("unbalanced, %s grid", x$occasions)))
  cat(sprintf("Log likelihood: %.6f; %d parameters; BIC (groups): %.4f\n",
              x$log_likelihood, x$n_parameters, x$bic_groups))
  cat(sprintf("Converged: %s after %d iterations; best of %d starts\n",
              x$converged, x$iterations, nrow(x$starts)))
  .latents_print_gaussian_diagnostics(x$gaussian_diagnostics)
  .multilpa_print_primary(x, rows = rows)
  invisible(x)
}

#' @rdname latents-model-methods
#' @export
logLik.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' @rdname latents-model-methods
#' @export
nobs.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  object$n_groups
}

#' @noRd
#' @export
summary.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  fields <- c("call", "n_observations", "n_groups", "n_profiles",
              "n_group_classes", "n_occasions", "occasions", "balanced",
              "variance_model", "covariance_model", "missing", "means",
              "variances", "standard_deviations", "covariances",
              "response_probabilities", "initial_probabilities",
              "transition_probabilities", "transition_counts",
              "profile_prevalence",
              "empty_transition_rows", "group_probabilities", "log_likelihood",
              "n_parameters", "aic", "bic", "bic_individual", "converged",
              "iterations", "boundary", "min_variance", "small_classes",
              "starts", "best_start", "n_best_replicated",
              "replication_tolerance", "effective_profile_counts",
              "effective_group_counts")
  # `covariances` exists only under a full covariance model and
  # `response_probabilities` only with categorical indicators. Every other
  # field is mandatory: subsetting a list by an absent name yields a NULL
  # element keyed NA rather than an error, so an absent field is caught here
  # instead of surfacing later as an `NA`-named hole in the summary.
  optional <- c("covariances", "response_probabilities")
  absent <- setdiff(setdiff(fields, optional), names(object))
  if (length(absent) > 0L) {
    stop(errorCondition(
      sprintf("The fit is missing fields the summary requires: %s.",
              paste(absent, collapse = ", ")),
      class = "latents_incomplete_fit", call = NULL))
  }
  result <- object[intersect(fields, names(object))]
  # The starts run to the end; under emEM the table also holds the screened.
  result$n_starts <- object$n_starts %||% nrow(object$starts)
  # Every table the fit can produce, built once here, so `get_results()` on the
  # summary serves the same tables the fit would and `print()` can show them
  # all without recomputing anything.
  result$tables <- get_results(object, "all")
  class(result) <- "summary_multilpa_transitions"
  result
}

#' @noRd
#' @export
as.data.frame.summary_multilpa_transitions <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa_transitions`" = inherits(x, "summary_multilpa_transitions"))
  .multilpa_coerce(x, row.names, list(...))
}

#' @noRd
#' @export
print.summary_multilpa_transitions <- function(x, digits = 4L, rows = 10L, ...) {
  stopifnot("`x` must be a `summary_multilpa_transitions` object" =
              inherits(x, "summary_multilpa_transitions"))
  .multilpa_check_print_arguments(digits, rows)
  cat(sprintf("Latent transition model: %d profiles and %d group %s\n",
              x$n_profiles, x$n_group_classes,
              if (x$n_group_classes == 1L) "class" else "classes"))
  cat(sprintf("Observations: %d; groups: %d; up to %d occasions (%s)\n",
              x$n_observations, x$n_groups, x$n_occasions,
              if (isTRUE(x$balanced)) "balanced" else
                sprintf("unbalanced, %s grid", x$occasions)))
  cat(sprintf("Parameters: %d; converged: %s\n", x$n_parameters, x$converged))
  cat(sprintf("Log likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              x$log_likelihood, x$aic, x$bic, x$bic_individual))
  cat(sprintf("Best likelihood replicated in %d/%d starts (absolute tolerance %.3g).\n",
              x$n_best_replicated, x$n_starts, x$replication_tolerance))
  if (!x$converged) cat("WARNING: best start did not converge.\n")
  if (x$boundary) cat("WARNING: at least one variance is at min_variance.\n")
  if (any(x$empty_transition_rows)) {
    cat("WARNING: a transition row is uniform by construction, not estimated.\n")
  }
  if (x$small_classes) cat("WARNING: an effective class membership is below one.\n")
  .multilpa_print_tables(x$tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(x$tables)
  invisible(x)
}

#' @rdname latents-model-methods
#' @export
coef.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  # One encoding and one set of labels for coef(), vcov() and
  # parameter_inference(), so their names agree entry for entry.
  estimates <- .multilpa_transition_encode(object, "natural")
  stats::setNames(unname(estimates), .multilpa_parameter_names(
    .multilpa_transition_labels(names(estimates),
                                .multilpa_transition_view(object))))
}

#' @rdname latents-plots
#' @export
plot.multilpa_transitions <- function(x, what = c("transitions", "profiles",
                                                  "bars", "heatmap", "responses",
                                                  "sequences", "sizes", "entropy",
                                                  "posteriors", "avepp", "all"),
                                      data = NULL,
                                      scale = c("raw", "standardized"),
                                      category = "last", labels = TRUE,
                                      cell_labels = TRUE, main = NULL,
                                      subtitle = NULL, ..., intervals = TRUE) {
  stopifnot("`x` must be a fitted `multilpa_transitions` model" =
              inherits(x, "multilpa_transitions"))
  .multilpa_check_plot_arguments(list(...), labels, intervals, cell_labels,
                                 category)
  what <- match.arg(what)
  scale <- match.arg(scale)
  .gg_require()
  if (identical(what, "all")) {
    return(.gg_every_view(x, match.call(), parent.frame()))
  }
  .multilpa_draw_view(x, what, data = data, scale = scale,
                      category = category, labels = labels, intervals = intervals,
                      cell_labels = cell_labels, main = main, subtitle = subtitle)
}
