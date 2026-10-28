# The mixed latent Markov model, written from its definition in base R.
#
# Used by suites/lta-mixed.R. A mixture of first-order latent Markov chains
# (Van de Pol & Langeheine 1990; Bartolucci, Farcomeni & Pennoni 2013, sec.
# 6.2): each sequence i belongs to one of K1 group classes u, with
# probability la_u; within class u its states follow a homogeneous chain with
# initial distribution piv_u and transition matrix Pi_u; each response depends
# only on the current state, through state-conditional category probabilities
# psi_j shared by every class and occasion:
#
#   f(y_i) = sum_u la_u sum_{s_1..s_T} piv_u(s_1) prod_{t>1} Pi_u(s_{t-1}, s_t)
#            prod_t prod_j psi_j(s_t, y_ijt).
#
# Each class's chain is evaluated by the scaled forward-backward recursion
# (Baum et al. 1970), over all sequences at once and one occasion at a time.
# Nothing here calls the package; every quantity the reference program (LMest)
# also reports is checked against it by the suite before anything else is
# trusted.
#
# Van de Pol, F. & Langeheine, R. (1990). Mixed Markov latent class models.
# Sociological Methodology, 20, 213-247.
# Bartolucci, F., Farcomeni, A. & Pennoni, F. (2013). Latent Markov Models for
# Longitudinal Data. Chapman & Hall/CRC.

#' Category codes of balanced sequences, by sequence and occasion
#'
#' @param frame A data frame ordered by sequence, then occasion, every
#'   sequence observed at every one of `occasions` occasions.
#' @param levels Named list of category labels per indicator, in the order the
#'   probability matrices' columns use.
#' @param occasions Number of occasions.
#' @return A list of integer matrices, one per indicator: rows are rows of
#'   `frame`, holding the 1-based category; and the layout (`occasions`,
#'   `n_sequences`).
.lta_mixed_codes <- function(frame, levels, occasions) {
  stopifnot("the frame must hold whole sequences of `occasions` rows" =
              nrow(frame) %% occasions == 0L)
  codes <- lapply(stats::setNames(nm = names(levels)), function(item) {
    code <- match(as.character(frame[[item]]), levels[[item]])
    stopifnot("every response must be one of the given categories" = !anyNA(code))
    code
  })
  list(codes = codes, occasions = occasions, n_sequences = nrow(frame) %/% occasions)
}

#' Emission densities of every occasion, scaled per row
#'
#' @param responses Named list of states-by-categories probability matrices.
#' @param layout From [.lta_mixed_codes()].
#' @return A list with `density`, one sequences-by-states matrix per occasion
#'   (`exp(log f(y_it | s) - peak_it)`), and `peak`, the sum over occasions of
#'   each sequence's `peak_it`, which the log-likelihood adds back.
.lta_mixed_emission <- function(responses, layout) {
  log_emission <- Reduce(`+`, lapply(names(layout$codes), function(item) {
    t(log(responses[[item]]))[layout$codes[[item]], , drop = FALSE]
  }))
  rows_at <- function(occasion) (seq_len(layout$n_sequences) - 1L) * layout$occasions + occasion
  by_occasion <- lapply(seq_len(layout$occasions), function(occasion) {
    block <- log_emission[rows_at(occasion), , drop = FALSE]
    peak <- apply(block, 1L, max)
    list(density = exp(block - peak), peak = peak)
  })
  list(density = lapply(by_occasion, `[[`, "density"),
       peak = Reduce(`+`, lapply(by_occasion, `[[`, "peak")))
}

#' The scaled forward-backward recursion of one class's chain
#'
#' `alpha_t` is normalised to sum to one by `scale_t`; `beta_t = (E_{t+1}
#' beta_{t+1}) Pi' / scale_{t+1}`, so that `alpha_t * beta_t` is the smoothed
#' state posterior and `sum_t log scale_t` the sequence log-likelihood (less
#' the emission peaks).
#'
#' @param initial Initial state distribution of the class.
#' @param transition Its transition matrix, rows the state moved from.
#' @param density The `density` of [.lta_mixed_emission()].
#' @param backward Whether to run the backward pass too.
#' @return A list with `alpha`, `scale` and `log_scale_sum` (one value per
#'   sequence), and `beta` when asked for.
.lta_mixed_chain <- function(initial, transition, density, backward = TRUE) {
  first <- sweep(density[[1L]], 2L, initial, "*")
  forward <- Reduce(function(previous, occasion) {
    alpha <- (previous$alpha %*% transition) * density[[occasion]]
    scale <- rowSums(alpha)
    list(alpha = alpha / scale, scale = scale)
  }, seq_along(density)[-1L], init = list(alpha = first / rowSums(first), scale = rowSums(first)),
  accumulate = TRUE)
  out <- list(alpha = lapply(forward, `[[`, "alpha"), scale = lapply(forward, `[[`, "scale"))
  out$log_scale_sum <- Reduce(`+`, lapply(out$scale, log))
  if (backward) {
    out$beta <- Reduce(function(occasion, following) {
      ((density[[occasion + 1L]] * following) %*% t(transition)) / out$scale[[occasion + 1L]]
    }, seq_len(length(density) - 1L), init = matrix(1, nrow(first), ncol(first)),
    right = TRUE, accumulate = TRUE)
  }
  out
}

#' Evaluate the mixed latent Markov model at given parameters
#'
#' @param parameters A list: `la` (group-class probabilities), `initial`
#'   (group classes by states), `transition` (states by states by group
#'   classes, rows the state moved from) and `responses` (named list of
#'   states-by-categories matrices).
#' @param layout From [.lta_mixed_codes()].
#' @return A list: `log_likelihood`; `sequence_log_likelihood`, log f(y_i);
#'   `group_posteriors`, sequences by group classes, P(u | y_i);
#'   `posteriors`, rows (in data order) by states, the marginal smoothed
#'   posterior `sum_u P(u | y_i) P(s_it | u, y_i)`; `counts`, states by states
#'   by group classes, the expected transitions `sum_i P(u | y_i) sum_t P(s_t
#'   = a, s_{t+1} = b | u, y_i)`; `initial_counts`, group classes by states,
#'   `sum_i P(u | y_i) P(s_1 | u, y_i)`; `prevalence`, group classes by
#'   states, each class's expected share of occasions in each state;
#'   `origin_occupancy`, states by group classes, the expected number of
#'   occasions before the last spent in each state, `sum_i P(u | y_i)
#'   sum_{t < T} P(s_t = a | u, y_i)`, which is the information a row of
#'   transitions out of `a` rests on; and
#'   `response_counts`, per indicator the states-by-categories expected counts
#'   under the marginal posteriors.
.lta_mixed_evaluate <- function(parameters, layout) {
  emission <- .lta_mixed_emission(parameters$responses, layout)
  group_classes <- seq_along(parameters$la)
  chains <- lapply(group_classes, function(u) {
    .lta_mixed_chain(parameters$initial[u, ], parameters$transition[, , u], emission$density)
  })
  # log la_u + log f(y_i | u), and the log-sum-exp over classes.
  joint <- vapply(group_classes, function(u) {
    log(parameters$la[[u]]) + chains[[u]]$log_scale_sum + emission$peak
  }, numeric(layout$n_sequences))
  joint <- matrix(joint, layout$n_sequences)
  top <- apply(joint, 1L, max)
  sequence_log_likelihood <- top + log(rowSums(exp(joint - top)))
  z <- exp(joint - sequence_log_likelihood)
  occasions <- seq_len(layout$occasions)
  smoothed <- lapply(group_classes, function(u) {
    lapply(occasions, function(t) chains[[u]]$alpha[[t]] * chains[[u]]$beta[[t]])
  })
  marginal <- lapply(occasions, function(t) {
    Reduce(`+`, lapply(group_classes, function(u) smoothed[[u]][[t]] * z[, u]))
  })
  n_states <- ncol(parameters$initial)
  # Stacked by occasion, row (t - 1) N + i holds sequence i at occasion t;
  # the data are sequence-major, row (i - 1) T + t.
  data_row <- as.vector(outer(seq_len(layout$n_sequences), occasions,
                              function(i, t) (i - 1L) * layout$occasions + t))
  posteriors <- matrix(0, layout$n_sequences * layout$occasions, n_states)
  posteriors[data_row, ] <- do.call(rbind, marginal)
  counts <- vapply(group_classes, function(u) {
    chain <- chains[[u]]
    pairs <- Reduce(`+`, lapply(occasions[-length(occasions)], function(t) {
      crossprod(chain$alpha[[t]] * z[, u],
                emission$density[[t + 1L]] * chain$beta[[t + 1L]] / chain$scale[[t + 1L]])
    }))
    pairs * parameters$transition[, , u]
  }, matrix(0, n_states, n_states))
  counts <- array(counts, c(n_states, n_states, length(group_classes)))
  initial_counts <- t(vapply(group_classes, function(u) colSums(smoothed[[u]][[1L]] * z[, u]),
                             numeric(n_states)))
  occupied <- t(vapply(group_classes, function(u) {
    Reduce(`+`, lapply(occasions, function(t) colSums(smoothed[[u]][[t]] * z[, u])))
  }, numeric(n_states)))
  origin_occupancy <- vapply(group_classes, function(u) {
    Reduce(`+`, lapply(occasions[-length(occasions)], function(t) {
      colSums(smoothed[[u]][[t]] * z[, u])
    }), numeric(n_states))
  }, numeric(n_states))
  response_counts <- lapply(stats::setNames(nm = names(layout$codes)), function(item) {
    categories <- ncol(parameters$responses[[item]])
    vapply(seq_len(categories), function(category) {
      colSums(posteriors[layout$codes[[item]] == category, , drop = FALSE])
    }, numeric(n_states))
  })
  list(log_likelihood = sum(sequence_log_likelihood),
       sequence_log_likelihood = sequence_log_likelihood, group_posteriors = z,
       posteriors = posteriors, counts = counts,
       initial_counts = matrix(initial_counts, length(group_classes)),
       prevalence = matrix(occupied / rowSums(occupied), length(group_classes)),
       origin_occupancy = matrix(origin_occupancy, n_states),
       response_counts = lapply(response_counts, function(m) matrix(m, n_states)))
}

#' The log-likelihood alone, by the forward pass
#' @inheritParams .lta_mixed_evaluate
#' @return A single number.
.lta_mixed_log_likelihood <- function(parameters, layout) {
  emission <- .lta_mixed_emission(parameters$responses, layout)
  joint <- vapply(seq_along(parameters$la), function(u) {
    log(parameters$la[[u]]) + emission$peak +
      .lta_mixed_chain(parameters$initial[u, ], parameters$transition[, , u],
                       emission$density, backward = FALSE)$log_scale_sum
  }, numeric(layout$n_sequences))
  joint <- matrix(joint, layout$n_sequences)
  top <- apply(joint, 1L, max)
  sum(top + log(rowSums(exp(joint - top))))
}

#' One EM step of the mixed latent Markov model
#'
#' The M-step of Van de Pol & Langeheine (1990) with equal sequence weights:
#' every distribution is its expected counts, normalised.
#'
#' @inheritParams .lta_mixed_evaluate
#' @param evaluated `.lta_mixed_evaluate(parameters, layout)`, when already
#'   computed.
#' @return The updated parameters.
.lta_mixed_em_step <- function(parameters, layout,
                               evaluated = .lta_mixed_evaluate(parameters, layout)) {
  normalise_rows <- function(m) m / rowSums(m)
  transition <- parameters$transition
  transition[] <- vapply(seq_along(parameters$la), function(u) {
    normalise_rows(evaluated$counts[, , u])
  }, matrix(0, nrow(transition), ncol(transition)))
  responses <- lapply(evaluated$response_counts, normalise_rows)
  responses <- Map(function(updated, old) {
    dimnames(updated) <- dimnames(old)
    updated
  }, responses, parameters$responses[names(responses)])
  list(la = colMeans(evaluated$group_posteriors),
       initial = normalise_rows(evaluated$initial_counts), transition = transition,
       responses = responses)
}

#' Every parameter in one vector, in a fixed order
#' @param parameters As for [.lta_mixed_evaluate()].
#' @return A numeric vector.
.lta_mixed_flatten <- function(parameters) {
  c(parameters$la, as.vector(parameters$initial), as.vector(parameters$transition),
    unlist(lapply(parameters$responses, as.vector), use.names = FALSE))
}

#' Iterate the EM map to its fixed point
#'
#' Steps are taken until no parameter moves by more than `threshold`; the
#' rule is on the parameters, because on a flat likelihood the log-likelihood
#' stops changing well before they do. The map is iterated in SQUAREM cycles
#' (Varadhan & Roland 2008, scheme S3) -- two EM steps, an extrapolation
#' along them and an EM step from the extrapolated point, kept only when every
#' probability stays positive and the log-likelihood does not fall -- because
#' plain EM contracts slowly along the group-class directions. The
#' extrapolation is an affine combination of parameter sets, so every
#' distribution still sums to one. The fixed point is accepted only when a
#' plain EM step from it moves nothing by more than `threshold`.
#'
#' @inheritParams .lta_mixed_evaluate
#' @param threshold Largest parameter change that counts as no move.
#' @param max_cycles Bound on the SQUAREM cycles.
#' @return A list with `parameters` and `cycles`. A classed warning is raised
#'   when the bound is reached.
.lta_mixed_fixed_point <- function(parameters, layout, threshold = 1e-12,
                                   max_cycles = 2000L) {
  # The inverse of .lta_mixed_flatten(), on the shapes of `like`.
  shaped <- function(vector, like) {
    sizes <- c(length(like$la), length(like$initial), length(like$transition),
               vapply(like$responses, length, integer(1)))
    pieces <- split(vector, rep(seq_along(sizes), sizes))
    like$la[] <- pieces[[1L]]
    like$initial[] <- pieces[[2L]]
    like$transition[] <- pieces[[3L]]
    like$responses <- Map(function(block, values) {
      block[] <- values
      block
    }, like$responses, pieces[-seq_len(3L)])
    like
  }
  step <- function(p) .lta_mixed_em_step(p, layout)
  state <- Reduce(function(state, cycle) {
    if (state$settled) return(state)
    p0 <- state$parameters
    one <- step(p0)
    v0 <- .lta_mixed_flatten(p0)
    v1 <- .lta_mixed_flatten(one)
    if (max(abs(v1 - v0)) < threshold) {
      return(list(parameters = p0, settled = TRUE, cycles = cycle))
    }
    two <- step(one)
    v2 <- .lta_mixed_flatten(two)
    r <- v1 - v0
    v <- v2 - 2 * v1 + v0
    alpha <- min(-sqrt(sum(r^2) / sum(v^2)), -1)
    jumped <- v0 - 2 * alpha * r + alpha^2 * v
    following <- two
    if (is.finite(alpha) && all(jumped > 0)) {
      candidate <- step(shaped(jumped, p0))
      if (.lta_mixed_log_likelihood(candidate, layout) >=
            .lta_mixed_log_likelihood(two, layout)) following <- candidate
    }
    list(parameters = following, settled = FALSE, cycles = cycle)
  }, seq_len(max_cycles), init = list(parameters = parameters, settled = FALSE, cycles = 0L))
  if (!state$settled) {
    warning(warningCondition(sprintf(
      "the reference EM was still moving after %d SQUAREM cycles", max_cycles),
      class = "multilpa_harness_no_converge"))
  }
  state[c("parameters", "cycles")]
}

#' Observed-information standard errors of every probability
#'
#' The log-likelihood is written in baseline-category logits of every
#' distribution (last entry the baseline), which have no boundary, and its
#' Hessian taken by Richardson extrapolation (numDeriv). The probabilities'
#' covariance follows by the delta method, `J V J'` with `J = diag(p) - p p'`
#' over the free logits of each distribution; at a maximum this is the inverse
#' observed information of the probabilities whatever free coordinates are
#' used, so it is the estimand of the package's observed-information errors.
#'
#' @inheritParams .lta_mixed_evaluate
#' @return A list shaped like `parameters` (`la`, `initial`, `transition`,
#'   `responses`) holding standard errors.
.lta_mixed_standard_errors <- function(parameters, layout) {
  # The distributions, each a vector of probabilities, in one fixed order.
  k1 <- length(parameters$la)
  k2 <- ncol(parameters$initial)
  distributions <- c(
    list(la = parameters$la),
    lapply(seq_len(k1), function(u) parameters$initial[u, ]),
    unlist(lapply(seq_len(k1), function(u) {
      lapply(seq_len(k2), function(from) parameters$transition[from, , u])
    }), recursive = FALSE),
    unlist(lapply(parameters$responses, function(block) {
      lapply(seq_len(k2), function(state) block[state, ])
    }), recursive = FALSE))
  widths <- vapply(distributions, length, integer(1))
  free <- widths - 1L
  ends <- cumsum(free)
  logits <- unlist(lapply(distributions, function(p) log(p[-length(p)] / p[[length(p)]])))
  softmax <- function(eta) {
    e <- exp(c(eta, 0) - max(c(eta, 0)))
    e / sum(e)
  }
  rebuild <- function(theta) {
    blocks <- lapply(seq_along(distributions), function(d) {
      softmax(theta[ends[[d]] - free[[d]] + seq_len(free[[d]])])
    })
    initial <- do.call(rbind, blocks[1L + seq_len(k1)])
    transition <- array(0, c(k2, k2, k1))
    transition[] <- vapply(seq_len(k1), function(u) {
      do.call(rbind, blocks[1L + k1 + (u - 1L) * k2 + seq_len(k2)])
    }, matrix(0, k2, k2))
    offset <- 1L + k1 + k1 * k2
    responses <- Map(function(block, item_index) {
      out <- do.call(rbind, blocks[offset + (item_index - 1L) * k2 + seq_len(k2)])
      dimnames(out) <- dimnames(block)
      out
    }, parameters$responses, seq_along(parameters$responses))
    list(la = blocks[[1L]], initial = initial, transition = transition,
         responses = responses)
  }
  stopifnot("the logit coordinates must reproduce the parameters" =
              max(abs(.lta_mixed_flatten(rebuild(logits)) -
                        .lta_mixed_flatten(parameters))) < 1e-12)
  hessian <- numDeriv::hessian(function(theta) .lta_mixed_log_likelihood(rebuild(theta), layout),
                               logits)
  covariance <- solve(-hessian)
  errors <- lapply(seq_along(distributions), function(d) {
    p <- distributions[[d]]
    jacobian <- (diag(p, length(p)) - tcrossprod(p))[, seq_len(free[[d]]), drop = FALSE]
    at <- ends[[d]] - free[[d]] + seq_len(free[[d]])
    sqrt(diag(jacobian %*% covariance[at, at, drop = FALSE] %*% t(jacobian)))
  })
  rebuilt <- rebuild(logits)
  rebuilt$la <- errors[[1L]]
  rebuilt$initial[] <- do.call(rbind, errors[1L + seq_len(k1)])
  rebuilt$transition[] <- vapply(seq_len(k1), function(u) {
    do.call(rbind, errors[1L + k1 + (u - 1L) * k2 + seq_len(k2)])
  }, matrix(0, k2, k2))
  offset <- 1L + k1 + k1 * k2
  rebuilt$responses <- Map(function(block, item_index) {
    block[] <- do.call(rbind, errors[offset + (item_index - 1L) * k2 + seq_len(k2)])
    block
  }, rebuilt$responses, seq_along(rebuilt$responses))
  rebuilt
}
