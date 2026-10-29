# Independent pre-scaling log-space implementation, retained only for tests.
.lta_log_reference <- function(emission, layout, log_initial, log_transition) {
  n_groups <- nrow(layout$slot)
  n_profiles <- ncol(log_initial)
  later <- seq_len(layout$n_occasions)[-1L]
  hold <- function(updated, previous, occasion) {
    outside <- !layout$within[, occasion]
    if (any(outside)) updated[outside, ] <- previous[outside, , drop = FALSE]
    updated
  }
  forward_step <- function(state, transition) {
    matrix(vapply(seq_len(n_profiles), function(to) {
      .latents_log_sum_exp(state + transition[, , to])
    }, numeric(n_groups)), n_groups, n_profiles)
  }
  backward_step <- function(state, transition) {
    matrix(vapply(seq_len(n_profiles), function(from) {
      .latents_log_sum_exp(state + matrix(transition[, from, ], n_groups, n_profiles))
    }, numeric(n_groups)), n_groups, n_profiles)
  }
  alpha <- Reduce(function(previous, t) {
    hold(forward_step(previous, log_transition[[t - 1L]]) + emission[[t]],
         previous, t)
  }, later, init = log_initial + emission[[1L]], accumulate = TRUE)
  beta <- Reduce(function(t, following) {
    hold(backward_step(following + emission[[t]], log_transition[[t - 1L]]),
         following, t)
  }, later, init = matrix(0, n_groups, n_profiles), right = TRUE, accumulate = TRUE)
  if (!is.list(alpha)) alpha <- list(alpha)
  if (!is.list(beta)) beta <- list(beta)
  list(alpha = alpha, beta = beta,
       log_scaled = .latents_log_sum_exp(alpha[[layout$n_occasions]]))
}


.lta_moments_reference <- function(pass, emission, layout, transition,
                                       weights, n_observations) {
  stopifnot(
    "`pass` must be the list `.multilpa_forward_backward()` returns" =
      is.list(pass) && all(c("alpha", "beta", "log_scaled") %in% names(pass)),
    "`weights` must give one group-class posterior per group" =
      is.numeric(weights) && length(weights) == nrow(layout$slot),
    "`transition` must be a square matrix of positive probabilities" =
      is.matrix(transition) && is.numeric(transition) &&
      nrow(transition) == ncol(transition) && all(transition > 0))
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
  # The expected count of a move is a sum of pair posteriors, and a pair
  # posterior is a probability: its logarithm is never positive, so it
  # exponentiates without overflow. Its three factors are not probabilities.
  # On a long sequence the forward term is scaled far above one exactly where
  # the backward term is far below it, so exponentiating them apart -- as this
  # did until the complete pair log probability was formed first -- gave
  # `Inf * 0` and returned NaN counts from a fit whose log likelihood, read off
  # the scaled forward pass alone, was correct. No per-occasion maximum is
  # subtracted because none is needed: the exponent is already bounded above by
  # zero, and any rescaling that kept the three factors apart would carry the
  # reciprocal of a transition probability and could overflow on its own.
  n_groups <- nrow(layout$slot)
  # One column per ordered pair of profiles, laid out column-major exactly as
  # the transition matrix is, so that the groups are summed out in one pass.
  origin <- rep(seq_len(n_profiles), times = n_profiles)
  destination <- rep(seq_len(n_profiles), each = n_profiles)
  log_transition <- rep(log(transition), each = n_groups)
  counts <- if (layout$n_occasions < 2L) {
    matrix(0, n_profiles, n_profiles)
  } else Reduce(`+`, lapply(seq_len(layout$n_occasions)[-1L], function(occasion) {
    active <- weights * layout$within[, occasion]
    log_from <- pass$alpha[[occasion - 1L]] - pass$log_scaled
    log_to <- pass$beta[[occasion]] + emission[[occasion]]
    joint <- log_from[, origin, drop = FALSE] +
      log_to[, destination, drop = FALSE] + log_transition
    matrix(colSums(active * exp(joint)), n_profiles, n_profiles)
  }))
  list(posterior = posterior, initial = colSums(occupancy[[1L]]),
       transition = counts)
}
