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
#' @param log_density Observation-by-profile matrix of measurement log densities.
#' @param layout Sequence layout from [.multilpa_sequence_layout()].
#' @return A list with `log_density` (one groups-by-profiles matrix per
#'   occasion, offset removed) and `offset` (one total per group).
#' @noRd
.multilpa_sequence_emission <- function(log_density, layout) {
  stopifnot(is.matrix(log_density), is.list(layout))
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
  offsets <- lapply(blocks, function(block) apply(block, 1L, max))
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
#' @param emission Offset-removed emission list from [.multilpa_sequence_emission()].
#' @param layout Sequence layout from [.multilpa_sequence_layout()].
#' @param initial Initial profile probabilities for one group class.
#' @param transition Profile transition matrix for one group class, rows from.
#' @return A list with `alpha`, `beta` (one groups-by-profiles matrix per
#'   occasion) and `log_scaled`, the offset-removed sequence log likelihoods.
#' @noRd
.multilpa_forward_backward <- function(emission, layout, initial, transition) {
  stopifnot(is.list(emission), is.numeric(initial), is.matrix(transition),
            nrow(transition) == length(initial),
            ncol(transition) == length(initial),
            all(initial > 0), all(transition > 0))
  n_groups <- nrow(layout$slot)
  n_profiles <- length(initial)
  n_occasions <- layout$n_occasions
  log_transition <- log(transition)
  later_occasions <- seq_len(n_occasions)[-1L]
  advance <- function(state, log_weights) {
    matrix(vapply(seq_len(n_profiles), function(profile) {
      .multilpa_log_sum_exp(sweep(state, 2L, log_weights[, profile], "+"))
    }, numeric(n_groups)), n_groups, n_profiles)
  }
  hold <- function(updated, previous, occasion) {
    outside <- !layout$within[, occasion]
    if (any(outside)) updated[outside, ] <- previous[outside, , drop = FALSE]
    updated
  }
  alpha <- Reduce(function(previous, occasion) {
    hold(advance(previous, log_transition) + emission[[occasion]],
         previous, occasion)
  }, later_occasions, init = sweep(emission[[1L]], 2L, log(initial), "+"),
  accumulate = TRUE)
  beta <- Reduce(function(occasion, later) {
    hold(advance(later + emission[[occasion]], t(log_transition)),
         later, occasion)
  }, later_occasions, init = matrix(0, n_groups, n_profiles),
  right = TRUE, accumulate = TRUE)
  if (!is.list(alpha)) alpha <- list(alpha)
  if (!is.list(beta)) beta <- list(beta)
  list(alpha = alpha, beta = beta,
       log_scaled = .multilpa_log_sum_exp(alpha[[n_occasions]]))
}

#' Posterior occupancy and transition counts for one group class
#'
#' @param pass Forward and backward recursions from [.multilpa_forward_backward()].
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
                                       weights, n_observations) {
  stopifnot(is.list(pass), is.numeric(weights),
            length(weights) == nrow(layout$slot))
  n_profiles <- ncol(transition)
  occupancy <- lapply(seq_len(layout$n_occasions), function(occasion) {
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
  } else Reduce(`+`, lapply(seq_len(layout$n_occasions)[-1L], function(occasion) {
    active <- weights * layout$within[, occasion]
    from <- exp(pass$alpha[[occasion - 1L]] - pass$log_scaled)
    to <- exp(pass$beta[[occasion]] + emission[[occasion]])
    transition * crossprod(from * active, to)
  }))
  list(posterior = posterior, initial = colSums(occupancy[[1L]]),
       transition = counts)
}

#' Nested expectation step for a latent transition model
#' @param x Observation-by-indicator matrix of continuous indicators.
#' @param layout Sequence layout.
#' @param parameters Model parameters, including `initial_probabilities` and
#'   `transition_probabilities`.
#' @param codes Optional integer matrix of categorical indicator codes.
#' @return Log likelihood, posteriors and expected sequence counts.
#' @noRd
.multilpa_transition_expectation <- function(x, layout, parameters, codes = NULL) {
  stopifnot(is.matrix(x), is.numeric(x), !any(is.infinite(x)), is.list(parameters))
  n_profiles <- nrow(parameters$means)
  n_types <- length(parameters$group_probabilities)
  n_groups <- nrow(layout$slot)
  gaussian <- if (ncol(x) > 0L && (anyNA(x) || !is.null(parameters$covariances))) {
    .multilpa_gaussian_moments(x, parameters)
  } else NULL
  log_density <- if (!is.null(gaussian)) gaussian$log_density else {
    vapply(seq_len(n_profiles), function(profile) {
      residuals <- sweep(x, 2L, parameters$means[profile, ], "-")
      -0.5 * rowSums(sweep(residuals^2, 2L, parameters$variances[profile, ], "/") +
        matrix(log(2 * pi) + log(parameters$variances[profile, ]),
               nrow(x), ncol(x), byrow = TRUE))
    }, numeric(nrow(x)))
  }
  if (!is.null(codes)) {
    log_density <- log_density +
      .multilpa_categorical_log_density(codes, parameters$response_probabilities)
  }
  emission <- .multilpa_sequence_emission(log_density, layout)
  passes <- lapply(seq_len(n_types), function(type) {
    .multilpa_forward_backward(emission$log_density, layout,
                               parameters$initial_probabilities[type, ],
                               matrix(parameters$transition_probabilities[, , type],
                                      n_profiles, n_profiles))
  })
  group_scores <- matrix(vapply(passes, `[[`, numeric(n_groups), "log_scaled"),
                         n_groups, n_types)
  weighted <- sweep(group_scores, 2L, log(parameters$group_probabilities), "+")
  group_log_likelihood <- .multilpa_log_sum_exp(weighted) + emission$offset
  group_posteriors <- exp(sweep(weighted, 1L, apply(weighted, 1L, max), "-"))
  group_posteriors <- group_posteriors / rowSums(group_posteriors)
  sequence <- lapply(seq_len(n_types), function(type) {
    .multilpa_sequence_moments(passes[[type]], emission$log_density, layout,
                               matrix(parameters$transition_probabilities[, , type],
                                      n_profiles, n_profiles),
                               group_posteriors[, type], nrow(x))
  })
  joint <- lapply(sequence, `[[`, "posterior")
  subject_posteriors <- Reduce(`+`, joint)
  log_likelihood <- sum(group_log_likelihood)
  if (!is.finite(log_likelihood) || any(!is.finite(subject_posteriors))) {
    stop("Non-finite likelihood or posterior probabilities.")
  }
  list(log_likelihood = log_likelihood,
       group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors,
       subject_posteriors = subject_posteriors, joint = joint,
       sequence = sequence, gaussian_moments = gaussian$moments)
}

#' Maximize the expected complete-data log likelihood of a transition model
#'
#' The measurement parameters are those of the cross-sectional model and are
#' maximized by the same function; only the initial and transition
#' probabilities are new, and both are bounded multinomial shares.
#'
#' @param x Observation-by-indicator matrix.
#' @param expectation Expected counts from [.multilpa_transition_expectation()].
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
                                              min_probability = 1e-10) {
  parameters <- .multilpa_maximization(x, expectation, variance_model,
                                       min_variance, covariance_model, codes,
                                       n_categories, min_probability)
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
  # The prevalence of each profile within a group class is implied by the
  # initial and transition probabilities; it is carried as a summary, not as a
  # free parameter of the model.
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
#' @return Parameters, the final expectation, convergence and history.
#' @noRd
.multilpa_transition_em <- function(x, layout, parameters, variance_model,
                                    min_variance, max_iter, tol,
                                    covariance_model = "diagonal", codes = NULL,
                                    n_categories = NULL, min_probability = 1e-10) {
  stopifnot(is.matrix(x), is.list(parameters), max_iter >= 0L, tol > 0,
            min_variance > 0, variance_model %in% c("varying", "equal"))
  expectation <- .multilpa_transition_expectation(x, layout, parameters, codes)
  history <- expectation$log_likelihood
  converged <- FALSE
  iteration <- 0L
  while (iteration < max_iter && !converged) {
    parameters <- .multilpa_transition_maximization(
      x, expectation, variance_model, min_variance, covariance_model, codes,
      n_categories, min_probability)
    updated <- .multilpa_transition_expectation(x, layout, parameters, codes)
    improvement <- updated$log_likelihood - expectation$log_likelihood
    if (improvement < -1e-10 * (1 + abs(expectation$log_likelihood))) {
      stop("EM likelihood decreased beyond numerical roundoff.")
    }
    converged <- abs(improvement) <= tol * (1 + abs(expectation$log_likelihood))
    iteration <- iteration + 1L
    history <- c(history, updated$log_likelihood)
    expectation <- updated
  }
  list(parameters = parameters, expectation = expectation,
       converged = converged, iterations = iteration, history = history)
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
#' @return Strictly positive initial parameters for EM.
#' @noRd
.multilpa_transition_initialize <- function(x, group_index, n_profiles, n_types,
                                            variance_model, min_variance,
                                            start_index,
                                            covariance_model = "diagonal",
                                            codes = NULL, n_categories = NULL,
                                            min_probability = 1e-10) {
  parameters <- .multilpa_initialize(x, group_index, n_profiles, n_types,
                                     variance_model, min_variance, start_index,
                                     covariance_model, codes, n_categories,
                                     min_probability)
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
                                                  covariance_model) {
  # The cross-sectional count carries the measurement model and the group-class
  # split. Its per-class profile prevalences are replaced here by a per-class
  # initial distribution and a per-class transition matrix.
  cross_sectional <- .multilpa_count_parameters(
    n_profiles, n_group_classes, n_continuous, n_categories, variance_model,
    covariance_model)
  cross_sectional + n_group_classes * n_profiles * (n_profiles - 1L)
}

#' Fit a latent transition model
#'
#' Fits latent profiles to repeated observations of the same group and
#' estimates the probabilities of moving between them from one occasion to the
#' next. The measurement model is that of [multilpa()] -- Gaussian, categorical
#' or mixed indicators, invariant across occasions -- so profiles keep the same
#' meaning at every occasion and a change of profile is a change of state
#' rather than a change of definition. With more than one group class, each
#' class carries its own initial distribution and its own transition matrix,
#' which is how groups that differ in their dynamics are separated from groups
#' that differ only in where they start.
#'
#' @inheritParams multilpa
#' @param n_profiles Integer number of individual profiles, at least two.
#' @param time Name of the column giving each observation's occasion within its
#'   group. Required: this model is defined by the ordering. Values must be
#'   complete and unique within each group.
#' @param n_group_classes Positive integer number of latent group classes. One
#'   gives an ordinary single-level latent transition model; more than one fits
#'   a mixture of transition models over groups.
#' @param occasions How a group's occasions are placed on the transition grid.
#'   `"observed"` numbers each group's own observations consecutively, so a
#'   transition always links two adjacent observations. `"grid"` places them on
#'   the grid of every position seen in the data, so a group that skips a
#'   position still consumes a transition across the gap and contributes no
#'   measurement information at it. The two agree whenever every group is
#'   observed at every position.
#' @return A `multilpa_transitions` object containing `means`, `variances`,
#'   optional `covariances` and `response_probabilities`,
#'   `initial_probabilities` (group classes by profiles),
#'   `transition_probabilities` (profiles by profiles by group classes, rows
#'   indexing the profile moved from), `group_probabilities`, posterior
#'   matrices, classifications, log likelihood, information criteria, restart
#'   diagnostics and convergence history. Use [transitions()] for the tidy
#'   transition table and [as.data.frame.multilpa_transitions()] for the other
#'   tables. No standard
#'   errors, likelihood-ratio tests or guarantees of global optimality are
#'   given for this model family.
#' @details Transitions are first order and homogeneous over occasions: the
#'   probability of moving from one profile to another does not depend on the
#'   occasion or on earlier profiles. Measurement parameters are shared across
#'   occasions and across group classes. A profile that no group occupies before
#'   its final occasion leaves its transition row without information; the row
#'   is then uniform by construction rather than estimated. Such a row is warned
#'   about when the model is fitted, is flagged by the `estimated` column of
#'   [transitions()] and is listed by `transitions(fit, estimated = FALSE)`.
#'   Profile labels are arbitrary and
#'   are not comparable across fits without alignment.
#' @references Collins, L. M., & Lanza, S. T. (2010). Latent class and latent
#'   transition analysis. Wiley.
#'
#'   Vermunt, J. K. (2003). Multilevel latent class models. Sociological
#'   Methodology, 33, 213--239. doi:10.1111/j.0081-1750.2003.t01-1-00131.x.
#' @seealso [transitions()] for the fitted transition probabilities,
#'   [sequences()] for the assignments in order, and [multilpa()] for the
#'   cross-sectional model this one shares its measurement parameters with.
#' @examples
#' # Forty people measured on six occasions, mostly staying where they are.
#' set.seed(7)
#' example_data <- do.call(rbind, lapply(seq_len(40), function(person) {
#'   path <- Reduce(function(previous, step) {
#'     if (stats::runif(1) < 0.8) previous else 3L - previous
#'   }, seq_len(5), init = 1L + person %% 2L, accumulate = TRUE)
#'   data.frame(person = person, wave = seq_len(6), state = unlist(path))
#' }))
#' centre <- ifelse(example_data$state == 1L, -2, 2)
#' example_data$score_a <- stats::rnorm(nrow(example_data), centre, 0.7)
#' example_data$score_b <- stats::rnorm(nrow(example_data), 0.75 * centre, 0.7)
#'
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' transitions(fit)
#' summary(fit)
#' @export
fit_transitions <- function(data, indicators, group, n_profiles, time,
                            n_group_classes = 1L,
                            variance_model = c("varying", "equal"),
                            n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                            min_variance = 1e-6, seed = NULL,
                            missing = c("error", "fiml"),
                            covariance_model = c("diagonal", "full"),
                            categorical = character(), min_probability = 1e-10,
                            occasions = c("observed", "grid")) {
  stopifnot(is.data.frame(data), is.character(indicators), is.character(group),
            "`categorical` must be a character vector of indicator names" =
              is.character(categorical) && !anyNA(categorical),
            "`min_probability` must be a single number in (0, 1)" =
              is.numeric(min_probability) && length(min_probability) == 1L &&
              is.finite(min_probability) && min_probability > 0 &&
              min_probability < 1,
            "`time` must be supplied; a transition model is defined by the ordering" =
              !is.null(time))
  call <- match.call()
  variance_model <- match.arg(variance_model)
  covariance_model <- match.arg(covariance_model)
  missing <- match.arg(missing)
  occasions <- match.arg(occasions)
  # The data contract is checked BEFORE the ordering column. Reversed, a missing
  # group column surfaced as an internal `split()` failure ("group length is 0
  # but data length > 0") instead of this package's own classed error.
  .multilpa_check_arguments(data, indicators, group, n_profiles, n_group_classes,
                            n_starts, max_iter, tol, min_variance,
                            min_probability, seed, categorical)
  time_values <- .multilpa_time_values(data, time, group, indicators)
  if (n_profiles < 2L) {
    stop(errorCondition(
      "`n_profiles` must be at least two; a single profile has nothing to move between.",
      class = "multilpa_bad_transition", call = NULL))
  }
  measurement <- .multilpa_prepare_indicators(data, indicators, categorical,
                                              missing, min_probability)
  continuous <- measurement$continuous
  codes <- measurement$codes
  n_categories <- measurement$n_categories
  x <- measurement$x
  groups <- .multilpa_prepare_groups(data[[group]])
  layout <- .multilpa_sequence_layout(groups$index, time_values, groups$n, occasions)
  if (layout$n_occasions < 2L) {
    stop(errorCondition(
      "No group is observed at two occasions, so no transition can be estimated.",
      class = "multilpa_bad_transition", call = NULL))
  }
  distinct_rows <- if (is.null(codes)) nrow(unique(x)) else
    nrow(unique(cbind(x, codes)))
  if (n_profiles > distinct_rows) {
    stop("n_profiles cannot exceed the number of distinct observed indicator rows.")
  }
  if (n_group_classes > groups$n) {
    stop("n_group_classes cannot exceed the number of groups.")
  }
  if (n_group_classes > 1L && all(groups$sizes == 1L)) {
    stop("Multiple group classes are not identifiable with only singleton groups.")
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  # Centering is a translation of the indicators, with no density Jacobian and
  # no change of scale; the fitted means are shifted back before they are read.
  centers <- if (ncol(x) > 0L) colMeans(x, na.rm = TRUE) else numeric(0)
  x <- sweep(x, 2L, centers, "-")
  if (any(!is.finite(x[!is.na(x)]^2))) {
    stop("Indicator scales overflow squared residuals; rescale the data.")
  }
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch({
      initial <- .multilpa_transition_initialize(
        x, groups$index, n_profiles, n_group_classes, variance_model,
        min_variance, start_index, covariance_model, codes, n_categories,
        min_probability)
      .multilpa_transition_em(x, layout, initial, variance_model, min_variance,
                              max_iter, tol, covariance_model, codes,
                              n_categories, min_probability)
    }, error = function(error) list(error = conditionMessage(error)))
  })
  valid <- vapply(attempts, function(attempt) is.null(attempt$error), logical(1))
  if (!any(valid)) {
    stop(sprintf("All %d starts failed: %s", n_starts,
                 paste(unique(vapply(attempts, `[[`, character(1), "error")),
                       collapse = "; ")))
  }
  scores <- vapply(attempts, function(attempt) {
    if (is.null(attempt$error)) attempt$expectation$log_likelihood else -Inf
  }, numeric(1))
  best_start <- which.max(scores)
  best <- attempts[[best_start]]
  parameters <- best$parameters
  parameters$means <- sweep(parameters$means, 2L, centers, "+")
  profile_names <- paste0("profile_", seq_len(n_profiles))
  type_names <- paste0("group_class_", seq_len(n_group_classes))
  parameters <- .multilpa_label_parameters(parameters, profile_names, continuous,
                                           categorical, measurement$encoded$levels)
  dimnames(parameters$initial_probabilities) <- list(type_names, profile_names)
  dimnames(parameters$profile_prevalence) <- list(type_names, profile_names)
  dimnames(parameters$transition_probabilities) <-
    list(profile_names, profile_names, type_names)
  names(parameters$group_probabilities) <- type_names
  subject_posteriors <- best$expectation$subject_posteriors
  group_posteriors <- best$expectation$group_posteriors
  dimnames(subject_posteriors) <- list(rownames(data), profile_names)
  dimnames(group_posteriors) <- list(groups$ids, type_names)
  n_parameters <- .multilpa_count_transition_parameters(
    n_profiles, n_group_classes, ncol(x), n_categories, variance_model,
    covariance_model)
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
  observed_per_row <- rowSums(!is.na(x)) +
    if (is.null(codes)) 0L else rowSums(!is.na(codes))
  n_informative <- sum(observed_per_row > 0L)
  boundary <- .multilpa_covariance_boundary(parameters, min_variance)
  empty_rows <- any(parameters$empty_transition_rows)
  small_classes <- any(colSums(subject_posteriors) < 1) ||
    any(colSums(group_posteriors) < 1)
  result <- c(parameters, list(
    call = call, indicators = indicators, continuous = continuous,
    categorical = categorical,
    categorical_levels = measurement$encoded$levels,
    min_probability = min_probability,
    indicator_data = as.matrix(measurement$frame), categorical_data = codes,
    group = group, group_ids = groups$ids, time = time, time_values = time_values,
    group_values = groups$values, group_index = groups$index,
    group_sizes = setNames(groups$sizes, groups$ids),
    n_observations = nrow(x), n_informative = n_informative,
    n_groups = groups$n, n_profiles = as.integer(n_profiles),
    n_group_classes = as.integer(n_group_classes),
    n_occasions = layout$n_occasions, occasions = occasions,
    sequence_lengths = setNames(layout$span, groups$ids),
    balanced = all(layout$span == layout$n_occasions) &&
      !anyNA(layout$slot),
    variance_model = variance_model, covariance_model = covariance_model,
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
    bic = -2 * log_likelihood + log(groups$n) * n_parameters,
    bic_groups = -2 * log_likelihood + log(groups$n) * n_parameters,
    bic_individual = -2 * log_likelihood + log(n_informative) * n_parameters,
    converged = best$converged, iterations = best$iterations,
    log_likelihood_history = best$history, starts = starts,
    best_start = best_start, n_failed_starts = sum(!valid), boundary = boundary,
    small_classes = small_classes,
    effective_profile_counts = colSums(subject_posteriors),
    effective_group_counts = colSums(group_posteriors),
    n_best_replicated = sum(valid & abs(scores - log_likelihood) <=
                              1e-6 * (1 + abs(log_likelihood))),
    replication_tolerance = 1e-6 * (1 + abs(log_likelihood))))
  class(result) <- "multilpa_transitions"
  if (any(!valid)) {
    warning(sprintf(paste("%d of %d starts failed; read their messages with",
                          "as.data.frame(fit, what = \"starts\")."),
                    sum(!valid), n_starts), call. = FALSE)
  }
  if (!best$converged && max_iter > 0L) {
    warning("The best start did not converge; increase max_iter and inspect starts.",
            call. = FALSE)
  }
  if (boundary) {
    warning(if (covariance_model == "full")
      "A covariance eigenvalue reached min_variance; this is a bound-active constrained fit." else
      "A variance reached min_variance; this is a bound-active constrained fit.",
      call. = FALSE)
  }
  if (empty_rows) {
    warning(paste("A profile is never occupied before a final occasion; its",
                  "transition row is uniform by construction, not estimated.",
                  "List the affected rows with",
                  "transitions(fit, estimated = FALSE)."),
            call. = FALSE)
  }
  if (small_classes) {
    warning("A profile or group class has effective membership below one.",
            call. = FALSE)
  }
  result
}

#' Fitted transition probabilities
#'
#' Returns the estimated probability of moving from each profile to each
#' profile between consecutive occasions, together with the expected number of
#' transitions each probability was estimated from. The diagonal is the
#' probability of staying, so a model whose profiles are stable and one whose
#' profiles are churning are told apart by reading this table, not the profile
#' prevalences.
#'
#' @param object A fitted `multilpa_transitions` model.
#' @param estimated Keep only the rows the data could estimate (`TRUE`), only
#'   the rows that are uniform by construction (`FALSE`), or all of them
#'   (`NULL`, the default).
#' @param stable Keep only the probabilities of staying in a profile (`TRUE`),
#'   only the probabilities of moving to a different one (`FALSE`), or all of
#'   them (`NULL`, the default). The two arguments combine, so
#'   `transitions(fit, estimated = TRUE, stable = FALSE)` is every estimated
#'   move out of a profile.
#' @return A base `data.frame` with one row per group class and ordered pair of
#'   profiles, and the columns `group_class`, `from`, `to`, `probability`, the
#'   `expected_count` of such transitions, `stable` (`TRUE` where `from` and
#'   `to` are the same profile, so the probability is one of staying),
#'   `estimated` (`FALSE` where the profile moved from is never occupied before
#'   a final occasion, so nothing was ever observed leaving it and the whole row
#'   is uniform by construction rather than estimated) and
#'   `group_class_probability`. `stable` and `estimated` are logical and never
#'   `NA`. Rows are ordered by group class, then by the profile moved from, then
#'   by the profile moved to; `estimated` and `stable` drop rows without
#'   reordering the ones they keep.
#'
#'   `group_class_probability` is the fitted share of groups in `group_class`,
#'   so it is constant within a group class and repeats down its rows. It is
#'   carried here because the overall transition matrix is the class-weighted
#'   average of the per-class ones, and that weight would otherwise have to be
#'   fetched out of the fitted object by hand.
#' @seealso [fit_transitions()] to fit the model, and [sequences()] for the
#'   assignments the transitions are estimated from.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' transitions(fit)
#' transitions(fit, stable = FALSE)
#' @export
transitions <- function(object, estimated = NULL, stable = NULL) {
  stopifnot(
    "`object` must be a fitted `multilpa_transitions` model" =
      inherits(object, "multilpa_transitions"),
    "`estimated` must be NULL, TRUE or FALSE" =
      .multilpa_optional_flag(estimated),
    "`stable` must be NULL, TRUE or FALSE" = .multilpa_optional_flag(stable))
  .multilpa_restrict_transitions(.multilpa_transition_frame(object),
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

#' Tidy a fitted latent transition model
#'
#' Returns any part of a fitted transition model as a tidy `data.frame`, so
#' that results can be printed, joined and written out without reaching into
#' the fitted object.
#'
#' @param x A fitted `multilpa_transitions` model.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"transitions"` gives exactly
#'   [transitions()], `"initial"` the initial profile probabilities and implied
#'   prevalences, `"profiles"` the Gaussian measurement model, `"responses"`
#'   the categorical measurement model, `"posteriors"` the individual
#'   posteriors, `"group_posteriors"` the group posteriors, `"sequence_lengths"`
#'   how many occasions each group contributes, `"starts"` the restart
#'   diagnostics, and `"information_criteria"`, `"classification"` and
#'   `"entropy"` exactly [information_criteria()], [classification_table()] and
#'   [entropy_table()].
#' @param format For the posterior tables, `"long"` (the default) gives one row
#'   per unit and class, and `"wide"` gives one row per unit with one column per
#'   class. The wide column set grows with the number of classes, so no caller
#'   can address it generically; it stays available because joining posteriors
#'   back onto the fitting data wants one row per unit. Supplying it with any
#'   other `what` is an error rather than a silent no-op.
#' @param ... Passed to the underlying accessor; `"transitions"` accepts
#'   `estimated` and `stable`, `"classification"` accepts `level`, and
#'   `"information_criteria"` accepts `definitions`. The per-assigned-class
#'   breakdown that `classification_table(detail = TRUE)` used to return is now
#'   its own verb, [average_posteriors()].
#' @return A base `data.frame` whose columns depend on `what`. `"initial"` has
#'   one row per group class and profile, with columns `group_class`,
#'   `profile`, `probability`, `prevalence` and `group_class_probability`.
#'   `"sequence_lengths"` has one row per group, with columns `group`,
#'   `group_class`, `occasions` (the span the transitions run over, which under
#'   `occasions = "grid"` counts positions the group skipped), `observations`
#'   (how many rows the group actually contributes) and `complete`. The
#'   remaining values match the tables named under `what`.
#' @seealso [transitions()], [fit_transitions()].
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' as.data.frame(fit)
#' as.data.frame(fit, what = "initial")
#' @export
as.data.frame.multilpa_transitions <- function(x, row.names = NULL,
                                               optional = FALSE,
                                               what = c("transitions", "initial",
                                                        "profiles", "responses",
                                                        "posteriors",
                                                        "group_posteriors",
                                                        "sequence_lengths",
                                                        "starts",
                                                        "information_criteria",
                                                        "classification",
                                                        "entropy"),
                                               format = c("long", "wide"), ...) {
  stopifnot("`x` must be a `multilpa_transitions` fit" =
              inherits(x, "multilpa_transitions"))
  what <- match.arg(what)
  format <- match.arg(format)
  stopifnot("`format` applies only to the posterior tables" =
              identical(format, "long") ||
              what %in% c("posteriors", "group_posteriors"))
  result <- switch(what,
    transitions = transitions(x, ...),
    initial = .multilpa_initial_frame(x),
    profiles = .multilpa_profile_frame(x, ...),
    responses = .multilpa_response_frame(x, ...),
    posteriors = .multilpa_posterior_frame(x, format),
    group_posteriors = .multilpa_group_posterior_frame(x, format),
    sequence_lengths = data.frame(
      group = x$group_values, group_class = x$group_classes,
      occasions = unname(x$sequence_lengths),
      observations = unname(tabulate(x$group_index, nbins = x$n_groups)),
      complete = unname(x$sequence_lengths) == x$n_occasions,
      row.names = NULL),
    starts = x$starts,
    information_criteria = information_criteria(x, ...),
    classification = classification_table(x, ...),
    entropy = entropy_table(x))
  row.names(result) <- row.names
  result
}

#' Summarize a fitted latent transition model
#'
#' Reports the model's size, how the occasions are laid out, its fit and its
#' convergence, and names the verbs that return the fitted quantities.
#'
#' @param x A fitted `multilpa_transitions` model.
#' @param ... Ignored.
#' @return `x`, invisibly; called for the summary it prints.
#' @seealso [transitions()], [as.data.frame.multilpa_transitions()].
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' print(fit)
#' @export
print.multilpa_transitions <- function(x, ...) {
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
  cat("Transition probabilities: transitions(x); other tables: as.data.frame(x)\n")
  invisible(x)
}

#' Log likelihood of a fitted latent transition model
#'
#' @param object A fitted `multilpa_transitions` model.
#' @param ... Ignored.
#' @return A `logLik` object carrying the maximized observed-data log
#'   likelihood, the free parameter count as `df`, and the number of groups as
#'   `nobs`. Groups are the independent units, because a group's occasions are
#'   dependent by construction in this model.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' logLik(fit)
#' @export
logLik.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' Number of independent units in a fitted latent transition model
#'
#' @param object A fitted `multilpa_transitions` model.
#' @param ... Ignored.
#' @return The number of observed groups. Occasions within a group are
#'   dependent by construction, so they are not independent observations.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' nobs(fit)
#' @export
nobs.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  object$n_groups
}

#' Summarize a fitted latent transition model
#'
#' @param object A fitted `multilpa_transitions` model.
#' @param ... Reserved for compatibility with `summary()`.
#' @return A `summary_multilpa_transitions` object carrying the measurement
#'   parameters, the initial and transition probabilities and counts, the
#'   effective class counts, the fit statistics and the restart diagnostics.
#'   Its tables are read with [as.data.frame.summary_multilpa_transitions()];
#'   `print()` reports the whole model.
#' @seealso [transitions()] for the transition probabilities as a tidy table.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' summary(fit)
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
      class = "multilpa_incomplete_fit", call = NULL))
  }
  result <- object[intersect(fields, names(object))]
  class(result) <- "summary_multilpa_transitions"
  result
}

#' Tables of a latent transition summary
#'
#' Returns the summary's own estimates as tidy `data.frame`s, so that a summary
#' can be printed, joined and written out without reaching into it.
#'
#' @param x A `summary_multilpa_transitions` object.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"transitions"` (the default) gives the
#'   same table as [transitions()] on the fit it summarizes, `"initial"` the
#'   initial profile probabilities and implied prevalences, and `"starts"` the
#'   restart diagnostics.
#' @param estimated,stable Passed to the same arguments of [transitions()] when
#'   `what = "transitions"`, and ignored otherwise.
#' @param ... Ignored.
#' @return A base `data.frame` whose columns depend on `what`: one row per group
#'   class and ordered pair of profiles for `"transitions"`, one row per group
#'   class and profile for `"initial"`, and one row per start for `"starts"`.
#' @seealso [summary.multilpa_transitions()], [transitions()].
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' as.data.frame(summary(fit))
#' as.data.frame(summary(fit), what = "initial")
#' @export
as.data.frame.summary_multilpa_transitions <- function(
    x, row.names = NULL, optional = FALSE,
    what = c("transitions", "initial", "starts"),
    estimated = NULL, stable = NULL, ...) {
  stopifnot(
    "`x` must be a `summary_multilpa_transitions` object" =
      inherits(x, "summary_multilpa_transitions"),
    "`estimated` must be NULL, TRUE or FALSE" =
      .multilpa_optional_flag(estimated),
    "`stable` must be NULL, TRUE or FALSE" = .multilpa_optional_flag(stable))
  what <- match.arg(what)
  result <- switch(what,
    transitions = .multilpa_restrict_transitions(
      .multilpa_transition_frame(x), estimated = estimated, stable = stable),
    initial = .multilpa_initial_frame(x),
    starts = x$starts)
  row.names(result) <- row.names
  result
}

#' Print a latent transition summary
#'
#' @param x A `summary_multilpa_transitions` object.
#' @param digits Number of printed significant digits.
#' @param ... Additional arguments passed to matrix printing.
#' @return The summary, invisibly; called for what it prints.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' print(summary(fit), digits = 3)
#' @export
print.summary_multilpa_transitions <- function(x, digits = 4L, ...) {
  stopifnot(inherits(x, "summary_multilpa_transitions"), is.numeric(digits),
            length(digits) == 1L, is.finite(digits), digits >= 1, digits <= 22)
  cat(sprintf("Latent transition model: %d profiles and %d group %s\n",
              x$n_profiles, x$n_group_classes,
              if (x$n_group_classes == 1L) "class" else "classes"))
  cat(sprintf("Observations: %d; groups: %d; up to %d occasions (%s)\n",
              x$n_observations, x$n_groups, x$n_occasions,
              if (isTRUE(x$balanced)) "balanced" else
                sprintf("unbalanced, %s grid", x$occasions)))
  cat(sprintf("Parameters: %d; converged: %s\n", x$n_parameters, x$converged))
  if (ncol(x$means) > 0L) {
    cat("\nProfile means:\n")
    print(x$means, digits = digits, ...)
    cat("\nProfile standard deviations:\n")
    print(x$standard_deviations, digits = digits, ...)
  }
  if (length(x$response_probabilities) > 0L) {
    cat("\nCategorical response probabilities:\n")
    print(x$response_probabilities, digits = digits, ...)
  }
  cat("\nInitial profile probabilities within each group class:\n")
  print(x$initial_probabilities, digits = digits, ...)
  cat("\nTransition probabilities (rows: profile moved from):\n")
  print(x$transition_probabilities, digits = digits, ...)
  cat("\nGroup-class probabilities:\n")
  print(x$group_probabilities, digits = digits, ...)
  cat("\nEffective individual memberships:\n")
  print(x$effective_profile_counts, digits = digits, ...)
  cat(sprintf("\nLog likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              x$log_likelihood, x$aic, x$bic, x$bic_individual))
  cat(sprintf("Best likelihood replicated in %d/%d starts (absolute tolerance %.3g).\n",
              x$n_best_replicated, nrow(x$starts), x$replication_tolerance))
  if (!x$converged) cat("WARNING: best start did not converge.\n")
  if (x$boundary) cat("WARNING: at least one variance is at min_variance.\n")
  if (any(x$empty_transition_rows)) {
    cat("WARNING: a transition row is uniform by construction, not estimated.\n")
  }
  if (x$small_classes) cat("WARNING: an effective class membership is below one.\n")
  cat("\nEstimates as tidy tables: as.data.frame(summary(fit))\n")
  invisible(x)
}

#' Fitted parameters of a latent transition model
#'
#' @param object A fitted `multilpa_transitions` model.
#' @param ... Reserved for compatibility with `coef()`.
#' @return A named numeric vector of every free parameter on its natural scale:
#'   profile means and variances, categorical response probabilities where
#'   present, each group class's initial profile probabilities, its transition
#'   probabilities, and the group-class probabilities.
#'
#'   Names follow the package-wide `level.parameter.outcome.term` grammar shared
#'   with [coef.multilpa()], so one pattern matches across fit classes: for
#'   example `measurement.mean.profile_1.reading` and
#'   `profile.transition_probability.profile_2.group_class_1:profile_1`, whose
#'   term names the origin the move is from. `level`, `parameter` and `outcome`
#'   never contain a dot, so everything after the third dot is the term and the
#'   name parses back into the columns [parameter_inference()] reports even when
#'   an indicator name itself contains a dot.
#'
#'   No standard errors accompany these: [vcov.multilpa_transitions()] refuses
#'   rather than returning an invalid matrix. [transitions()] and
#'   [as.data.frame.multilpa_transitions()] give the same quantities as tidy
#'   tables, which is the form to prefer.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' coef(fit)
#' @export
coef.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  # Names follow the package-wide `level.parameter.outcome.term` grammar, which is
  # parameter_inference()'s tidy decomposition written on one line. level, parameter
  # and outcome never contain a dot, so everything after the third dot is the term
  # and the name parses back even when an indicator name contains a dot.
  measurement <- c(
    stats::setNames(as.vector(t(object$means)),
                    sprintf("measurement.mean.%s.%s",
                            rep(rownames(object$means), each = ncol(object$means)),
                            rep(colnames(object$means), nrow(object$means)))),
    stats::setNames(as.vector(t(object$variances)),
                    sprintf("measurement.variance.%s.%s",
                            rep(rownames(object$variances), each = ncol(object$variances)),
                            rep(colnames(object$variances), nrow(object$variances)))))
  responses <- unlist(lapply(names(object$response_probabilities), function(indicator) {
    block <- object$response_probabilities[[indicator]]
    stats::setNames(as.vector(t(block)),
                    sprintf("measurement.response.%s.%s:%s",
                            rep(rownames(block), each = ncol(block)), indicator,
                            rep(colnames(block), nrow(block))))
  }), use.names = TRUE)
  moves <- transitions(object)
  c(measurement, responses,
    stats::setNames(as.vector(t(object$initial_probabilities)),
                    sprintf("profile.initial_probability.%s.%s",
                            rep(colnames(object$initial_probabilities),
                                object$n_group_classes),
                            rep(rownames(object$initial_probabilities),
                                each = object$n_profiles))),
    stats::setNames(moves$probability,
                    sprintf("profile.transition_probability.profile_%d.group_class_%d:profile_%d",
                            moves$to, moves$group_class, moves$from)),
    stats::setNames(object$group_probabilities,
                    sprintf("group.probability.group_class_%d",
                            seq_along(object$group_probabilities))))
}

#' Standard errors are not available for a latent transition model
#'
#' @param object A fitted `multilpa_transitions` model.
#' @param ... Ignored.
#' @return Nothing; `vcov()`, `confint()` and `parameter_inference()` all raise
#'   a `multilpa_no_inference` condition on a latent transition fit. The
#'   analytic score and Jacobian this package uses for observed-information and
#'   sandwich standard errors do not yet cover the initial and transition
#'   multinomial logits, so no interval is reported rather than an invalid one.
#'   `confint()` refuses explicitly instead of letting `confint.default()`
#'   reach `vcov()` and refuse by accident.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
#' )
#' example_data$score_a <- stats::rnorm(nrow(example_data))
#' example_data$score_b <- stats::rnorm(nrow(example_data))
#' fit <- fit_transitions(example_data, c("score_a", "score_b"), "person",
#'                        n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
#' tryCatch(confint(fit), multilpa_no_inference = function(condition) {
#'   conditionMessage(condition)
#' })
#' @export
#' @importFrom stats vcov
vcov.multilpa_transitions <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  .multilpa_refuse_transition_inference()
}

#' @rdname vcov.multilpa_transitions
#' @param data Ignored; present for compatibility with the generic.
#' @export
parameter_inference.multilpa_transitions <- function(object, data, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  .multilpa_refuse_transition_inference()
}

#' @rdname vcov.multilpa_transitions
#' @param parm Ignored; present for compatibility with the generic.
#' @param level Ignored; present for compatibility with the generic.
#' @export
#' @importFrom stats confint
confint.multilpa_transitions <- function(object, parm, level = 0.95, ...) {
  stopifnot(inherits(object, "multilpa_transitions"))
  .multilpa_refuse_transition_inference()
}

#' Refuse to report standard errors for a latent transition model
#'
#' One definition so that `vcov()`, `confint()` and `parameter_inference()`
#' cannot drift apart in what they say or in the class they raise.
#'
#' @return Nothing; always raises `multilpa_no_inference`.
#' @noRd
.multilpa_refuse_transition_inference <- function() {
  stop(errorCondition(
    paste("Standard errors are not available for a latent transition model.",
          "Read the estimates with transitions() and as.data.frame()."),
    class = "multilpa_no_inference", call = NULL))
}

#' No plot method is defined for a latent transition model
#'
#' @param x A fitted `multilpa_transitions` model.
#' @param ... Ignored.
#' @return Nothing; always raises a `multilpa_no_plot` condition, rather than
#'   letting the fit fall through to the default method and fail obscurely.
#' @examples
#' # plot() on a transition fit raises multilpa_no_plot by design.
#' @export
plot.multilpa_transitions <- function(x, ...) {
  stopifnot(inherits(x, "multilpa_transitions"))
  stop(errorCondition(
    paste("No plot method is defined for a latent transition model.",
          "Use transitions() for the transition probabilities, and",
          "plot(what = \"sequences\") on a multilpa() fit for profile paths."),
    class = "multilpa_no_plot", call = NULL))
}
