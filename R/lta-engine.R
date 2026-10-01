# General latent transition engine. Every transition from occasion t - 1 to t
# of group j has its own multinomial logit, driven by a design row Z[j, t, ]:
# an intercept alone gives the homogeneous model; occasion indicators give a
# transition matrix per occasion; covariate columns (fixed or changing over
# occasions) give covariate-dependent transitions. The initial profile
# distribution is a multinomial logit on a group-level design. The reference
# outcome of a transition logit is staying in the origin profile, so a
# coefficient is the log odds of moving to that profile rather than staying.
# Measurement is the multilpa() model, invariant across occasions or with its
# own means and variances at each occasion.

#' Log-sum-exp over the last dimension that tolerates rows of -Inf
#' @noRd
.lta_log_sum_exp <- function(values) {
  top <- .multilpa_row_max(values)
  finite <- is.finite(top)
  out <- rep(-Inf, nrow(values))
  if (any(finite)) {
    shifted <- exp(values[finite, , drop = FALSE] - top[finite])
    out[finite] <- top[finite] + log(rowSums(shifted))
  }
  out
}

#' Covariate columns as a numeric matrix, factors expanded to indicators
#' @return A matrix with one row per data row and named columns, or `NULL`.
#' @noRd
.lta_covariate_matrix <- function(data, covariates) {
  if (length(covariates) == 0L) return(NULL)
  missing_columns <- setdiff(covariates, names(data))
  if (length(missing_columns) > 0L) {
    stop(errorCondition(sprintf("Covariates %s are not columns of `data`.",
                                paste(sprintf("`%s`", missing_columns), collapse = ", ")),
                        class = "latents_bad_data", call = NULL))
  }
  frame <- data[, covariates, drop = FALSE]
  if (anyNA(frame)) {
    stop(errorCondition("Transition and initial covariates must be complete.",
                        class = "latents_bad_data", call = NULL))
  }
  design <- stats::model.matrix(stats::reformulate(covariates), frame)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (!all(is.finite(design))) {
    stop(errorCondition("Covariates must be finite.", class = "latents_bad_data",
                        call = NULL))
  }
  design
}

#' Initial and transition designs on the sequence layout
#'
#' @param layout Sequence layout.
#' @param covariate_rows Row-level covariate matrix for transitions, or `NULL`.
#' @param initial_rows Row-level covariate matrix for the initial
#'   distribution, or `NULL`; read at each group's first occasion.
#' @param transitions `"homogeneous"` (one intercept) or `"occasion"` (one
#'   intercept per occasion from the second on).
#' @return A list with `initial` (groups x p0), `transition` (a list over
#'   occasions 2..T of groups x p matrices; transition covariates are read at
#'   the destination occasion, zero where it is not observed) and the column
#'   names.
#' @noRd
.lta_designs <- function(layout, covariate_rows, initial_rows, transitions) {
  n_groups <- nrow(layout$slot)
  first_row <- layout$slot[cbind(seq_len(n_groups), 1L)]
  initial <- cbind(matrix(1, n_groups, 1L, dimnames = list(NULL, "(Intercept)")),
                   if (!is.null(initial_rows)) initial_rows[first_row, , drop = FALSE])
  occasion_names <- paste0("occasion_", seq_len(layout$n_occasions)[-1L])
  transition <- lapply(seq_len(layout$n_occasions)[-1L], function(t) {
    intercepts <- if (identical(transitions, "occasion")) {
      matrix(as.numeric(seq_len(layout$n_occasions)[-1L] == t), n_groups,
             length(occasion_names), byrow = TRUE,
             dimnames = list(NULL, occasion_names))
    } else matrix(1, n_groups, 1L, dimnames = list(NULL, "(Intercept)"))
    if (is.null(covariate_rows)) return(intercepts)
    rows <- layout$slot[, t]
    values <- matrix(0, n_groups, ncol(covariate_rows),
                     dimnames = list(NULL, colnames(covariate_rows)))
    observed <- !is.na(rows)
    values[observed, ] <- covariate_rows[rows[observed], , drop = FALSE]
    cbind(intercepts, values)
  })
  list(initial = initial, transition = transition)
}

#' Log transition probabilities of one class at one occasion
#'
#' @param design Groups x p design at that occasion.
#' @param coefficients p x (K - 1) x K array for this class: for origin k, the
#'   logits of the K - 1 destinations other than k, in profile order.
#' @return A groups x K x K array of log probabilities, origin by destination.
#' @noRd
.lta_log_transition <- function(design, coefficients) {
  n_profiles <- dim(coefficients)[3L]
  n_groups <- nrow(design)
  out <- array(0, c(n_groups, n_profiles, n_profiles))
  invisible(lapply(seq_len(n_profiles), function(k) {
    others <- setdiff(seq_len(n_profiles), k)
    # Staying is the reference category: .multilpa_log_softmax() puts it last.
    log_p <- .multilpa_log_softmax(design, matrix(coefficients[, , k], ncol(design)))
    out[, k, others] <<- log_p[, seq_along(others), drop = FALSE]
    out[, k, k] <<- log_p[, length(others) + 1L]
  }))
  out
}

#' Log initial probabilities of one class: groups x K (reference: last profile)
#' @noRd
.lta_log_initial <- function(design, coefficients) {
  .multilpa_log_softmax(design, coefficients)
}

#' Forward-backward with group- and occasion-specific transitions
#'
#' @param emission Offset-removed emission list (one groups x K matrix per
#'   occasion).
#' @param layout Sequence layout.
#' @param log_initial Groups x K log initial probabilities.
#' @param log_transition List over occasions 2..T of groups x K x K arrays.
#' @return A list with `alpha`, `beta` and `log_scaled`.
#' @noRd
.lta_forward_backward <- function(emission, layout, log_initial, log_transition) {
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
      .lta_log_sum_exp(state + transition[, , to])
    }, numeric(n_groups)), n_groups, n_profiles)
  }
  backward_step <- function(state, transition) {
    matrix(vapply(seq_len(n_profiles), function(from) {
      .lta_log_sum_exp(state + matrix(transition[, from, ], n_groups, n_profiles))
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
       log_scaled = .lta_log_sum_exp(alpha[[layout$n_occasions]]))
}

#' Posterior occupancy and pair probabilities of one class, per group and
#' occasion, weighted by the group's class posterior
#' @return A list with `posterior` (rows x K), `initial` (groups x K) and
#'   `pairs` (list over occasions 2..T of groups x K x K arrays).
#' @noRd
.lta_moments <- function(pass, emission, layout, log_transition, weights,
                         n_observations) {
  n_profiles <- ncol(pass$alpha[[1L]])
  n_groups <- nrow(layout$slot)
  occupancy <- lapply(seq_len(layout$n_occasions), function(t) {
    exp(pass$alpha[[t]] + pass$beta[[t]] - pass$log_scaled) * weights
  })
  posterior <- matrix(0, n_observations, n_profiles)
  invisible(lapply(seq_len(layout$n_occasions), function(t) {
    rows <- layout$slot[, t]
    observed <- !is.na(rows)
    if (any(observed)) posterior[rows[observed], ] <<- occupancy[[t]][observed, , drop = FALSE]
  }))
  pairs <- lapply(seq_len(layout$n_occasions)[-1L], function(t) {
    active <- weights * layout$within[, t]
    from <- pass$alpha[[t - 1L]] - pass$log_scaled
    to <- pass$beta[[t]] + emission[[t]]
    # The pair log probability is formed whole before exponentiating, so the
    # scaled forward and backward factors cannot overflow apart.
    joint <- array(from, c(n_groups, n_profiles, n_profiles)) +
      aperm(array(to, c(n_groups, n_profiles, n_profiles)), c(1L, 3L, 2L)) +
      log_transition[[t - 1L]]
    exp(joint) * active
  })
  list(posterior = posterior, initial = occupancy[[1L]], pairs = pairs)
}

#' Pair-state moments back to profile moments
#'
#' Occupancy and initial counts sum over the previous profile; the first move
#' becomes first-order counts (origin j, destination l); every later move is
#' kept by origin pair (i, j) and destination l.
#' @noRd
.lta_reduce_moments <- function(moments, n_states) {
  n_pairs <- n_states * n_states
  current <- rep(seq_len(n_states), times = n_states)
  collapse <- function(values) {
    # Sum the columns of each current profile.
    vapply(seq_len(n_states), function(j) rowSums(values[, current == j, drop = FALSE]),
           numeric(nrow(values)))
  }
  pair_to <- function(i, j) (j - 1L) * n_states + seq_len(n_states)
  first <- moments$pairs[[1L]]
  first_order <- array(0, c(dim(first)[1L], n_states, n_states))
  invisible(lapply(seq_len(n_pairs), function(s) {
    i <- (s - 1L) %/% n_states + 1L
    j <- (s - 1L) %% n_states + 1L
    first_order[, j, ] <<- first_order[, j, ] + first[, s, pair_to(i, j)]
  }))
  later <- lapply(moments$pairs[-1L], function(pair) {
    out <- array(0, c(dim(pair)[1L], n_pairs, n_states))
    invisible(lapply(seq_len(n_pairs), function(s) {
      i <- (s - 1L) %/% n_states + 1L
      j <- (s - 1L) %% n_states + 1L
      out[, s, ] <<- pair[, s, pair_to(i, j)]
    }))
    out
  })
  list(posterior = matrix(collapse(moments$posterior), nrow(moments$posterior)),
       initial = matrix(collapse(moments$initial), nrow(moments$initial)),
       pairs = list(first_order), pairs2 = later)
}

#' Measurement log densities, rows x K, invariant or per occasion
#' @noRd
.lta_log_density <- function(x, codes, parameters, occasion_of_row) {
  gaussian_density <- function(rows, means, variances) {
    vapply(seq_len(nrow(means)), function(k) {
      residuals <- sweep(x[rows, , drop = FALSE], 2L, means[k, ], "-")
      -0.5 * rowSums(sweep(residuals^2, 2L, variances[k, ], "/") +
                       matrix(log(2 * pi) + log(variances[k, ]), length(rows),
                              ncol(x), byrow = TRUE))
    }, numeric(length(rows)))
  }
  n_profiles <- nrow(parameters$measurement[[1L]]$means)
  out <- matrix(0, nrow(x), n_profiles)
  invisible(lapply(seq_along(parameters$measurement), function(t) {
    rows <- if (length(parameters$measurement) == 1L) seq_len(nrow(x)) else
      which(occasion_of_row == t)
    if (length(rows) == 0L) return(NULL)
    block <- parameters$measurement[[t]]
    if (ncol(x) > 0L) {
      out[rows, ] <<- matrix(gaussian_density(rows, block$means, block$variances),
                             length(rows), n_profiles)
    }
    if (!is.null(codes)) {
      out[rows, ] <<- out[rows, , drop = FALSE] + .multilpa_categorical_log_density(
        codes[rows, , drop = FALSE], block$response_probabilities)
    }
  }))
  out
}

#' Second-order transitions as a first-order chain on profile pairs
#'
#' State (i, j) at occasion t holds the profile at t - 1 (i) and at t (j);
#' index (i - 1) K + j. Occasion 1 occupies only the diagonal pairs (j, j)
#' with the initial probabilities; the move to occasion 2 uses the
#' first-order logits; every later move (i, j) -> (j, l) uses the
#' second-order logits of origin pair (i, j). Emissions read the current
#' profile j.
#' @return A list with the augmented `log_initial` (groups x K^2) and
#'   `log_transition` (list of groups x K^2 x K^2 arrays, -Inf where a move is
#'   impossible).
#' @noRd
.lta_augment <- function(log_initial, first_order, second_order) {
  n_groups <- nrow(log_initial)
  n_states <- ncol(log_initial)
  n_pairs <- n_states * n_states
  pair <- function(i, j) (i - 1L) * n_states + j
  initial <- matrix(-Inf, n_groups, n_pairs)
  initial[, pair(seq_len(n_states), seq_len(n_states))] <- log_initial
  expand <- function(per_origin) {
    # per_origin(i, j) gives the groups x K log probabilities of the next
    # profile l from pair (i, j).
    out <- array(-Inf, c(n_groups, n_pairs, n_pairs))
    cells <- expand.grid(i = seq_len(n_states), j = seq_len(n_states))
    invisible(lapply(seq_len(nrow(cells)), function(r) {
      i <- cells$i[r]
      j <- cells$j[r]
      out[, pair(i, j), pair(j, seq_len(n_states))] <<- per_origin(i, j)
    }))
    out
  }
  transition <- c(
    list(expand(function(i, j) matrix(first_order[[1L]][, j, ], n_groups))),
    lapply(second_order, function(log_p) {
      expand(function(i, j) matrix(log_p[, pair(i, j), ], n_groups))
    }))
  list(log_initial = initial, log_transition = transition)
}

#' Log second-order transition probabilities of one class at one occasion:
#' groups x K^2 origin pairs x K destinations (reference: staying in j)
#' @noRd
.lta_log_transition2 <- function(design, coefficients) {
  n_states <- dim(coefficients)[2L] + 1L
  n_groups <- nrow(design)
  out <- array(0, c(n_groups, n_states * n_states, n_states))
  invisible(lapply(seq_len(n_states * n_states), function(s) {
    j <- (s - 1L) %% n_states + 1L
    others <- setdiff(seq_len(n_states), j)
    log_p <- .multilpa_log_softmax(design, matrix(coefficients[, , s], ncol(design)))
    out[, s, others] <<- log_p[, seq_along(others), drop = FALSE]
    out[, s, j] <<- log_p[, length(others) + 1L]
  }))
  out
}

#' E-step of the general transition model
#' @noRd
.lta_expectation <- function(x, codes, layout, designs, parameters, occasion_of_row) {
  n_types <- length(parameters$group_probabilities)
  n_groups <- nrow(layout$slot)
  log_density <- .lta_log_density(x, codes, parameters, occasion_of_row)
  emission <- .multilpa_sequence_emission(log_density, layout)
  second <- !is.null(parameters$transition2)
  n_states <- nrow(parameters$measurement[[1L]]$means)
  # Emissions on the augmented pair states read the current profile.
  augmented_emission <- if (second) lapply(emission$log_density, function(block) {
    block[, rep(seq_len(n_states), times = n_states), drop = FALSE]
  }) else NULL
  per_class <- lapply(seq_len(n_types), function(h) {
    log_initial <- .lta_log_initial(designs$initial,
                                    matrix(parameters$initial[, , h], ncol(designs$initial)))
    log_transition <- lapply(designs$transition, function(design) {
      .lta_log_transition(design, array(parameters$transition[, , , h],
                                        dim(parameters$transition)[1:3]))
    })
    if (!second) {
      return(list(log_initial = log_initial, log_transition = log_transition,
                  pass = .lta_forward_backward(emission$log_density, layout,
                                               log_initial, log_transition)))
    }
    log_transition2 <- lapply(designs$transition[-1L], function(design) {
      .lta_log_transition2(design, array(parameters$transition2[, , , h],
                                         dim(parameters$transition2)[1:3]))
    })
    augmented <- .lta_augment(log_initial, log_transition, log_transition2)
    list(log_initial = log_initial, log_transition = log_transition,
         log_transition2 = log_transition2, augmented = augmented,
         pass = .lta_forward_backward(augmented_emission, layout,
                                      augmented$log_initial, augmented$log_transition))
  })
  scores <- matrix(vapply(per_class, function(c) c$pass$log_scaled, numeric(n_groups)),
                   n_groups, n_types)
  weighted <- sweep(scores, 2L, log(parameters$group_probabilities), "+")
  group_log_likelihood <- .multilpa_log_sum_exp(weighted) + emission$offset
  group_posteriors <- exp(weighted - .multilpa_row_max(weighted))
  group_posteriors <- group_posteriors / rowSums(group_posteriors)
  moments <- lapply(seq_len(n_types), function(h) {
    if (!second) {
      return(.lta_moments(per_class[[h]]$pass, emission$log_density, layout,
                          per_class[[h]]$log_transition, group_posteriors[, h],
                          nrow(x)))
    }
    .lta_reduce_moments(.lta_moments(per_class[[h]]$pass, augmented_emission, layout,
                                     per_class[[h]]$augmented$log_transition,
                                     group_posteriors[, h], nrow(x)), n_states)
  })
  joint <- lapply(moments, `[[`, "posterior")
  subject_posteriors <- Reduce(`+`, joint)
  log_likelihood <- sum(group_log_likelihood)
  if (!is.finite(log_likelihood) || any(!is.finite(subject_posteriors))) {
    stop("Non-finite likelihood or posterior probabilities.")
  }
  list(log_likelihood = log_likelihood, group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors, subject_posteriors = subject_posteriors,
       joint = joint, moments = moments, per_class = per_class)
}

#' M-step of the general transition model
#' @noRd
.lta_maximization <- function(x, codes, layout, designs, expectation, parameters,
                              occasion_of_row, variance_model, min_variance,
                              n_categories, min_probability) {
  n_types <- length(parameters$group_probabilities)
  n_profiles <- ncol(expectation$subject_posteriors)
  # Measurement: the cross-sectional M-step, on all rows or occasion by occasion.
  measurement <- lapply(seq_along(parameters$measurement), function(t) {
    rows <- if (length(parameters$measurement) == 1L) seq_len(nrow(x)) else
      which(occasion_of_row == t)
    sub <- list(subject_posteriors = expectation$subject_posteriors[rows, , drop = FALSE],
                group_posteriors = expectation$group_posteriors,
                joint = lapply(expectation$joint, function(j) j[rows, , drop = FALSE]))
    updated <- .multilpa_maximization(x[rows, , drop = FALSE], sub, variance_model,
                                      min_variance, "diagonal",
                                      if (is.null(codes)) NULL else codes[rows, , drop = FALSE],
                                      n_categories, min_probability)
    updated[intersect(c("means", "variances", "response_probabilities"), names(updated))]
  })
  initial <- array(vapply(seq_len(n_types), function(h) {
    counts <- expectation$moments[[h]]$initial
    keep <- rowSums(counts) > 0
    .multilpa_weighted_logits(designs$initial[keep, , drop = FALSE],
                              counts[keep, , drop = FALSE],
                              matrix(parameters$initial[, , h],
                                     ncol(designs$initial)))$coefficients
  }, numeric(ncol(designs$initial) * (n_profiles - 1L))),
  c(ncol(designs$initial), n_profiles - 1L, n_types))
  p <- ncol(designs$transition[[1L]])
  second <- !is.null(parameters$transition2)
  # With second-order transitions the first-order logits govern only the move
  # into occasion 2.
  stacked_design <- if (second) designs$transition[[1L]] else
    do.call(rbind, designs$transition)
  transition <- array(vapply(seq_len(n_types), function(h) {
    vapply(seq_len(n_profiles), function(k) {
      others <- setdiff(seq_len(n_profiles), k)
      # Counts of moves out of k, destination columns with staying last.
      counts <- do.call(rbind, lapply(expectation$moments[[h]]$pairs, function(pair) {
        matrix(pair[, k, c(others, k)], nrow(pair))
      }))
      keep <- rowSums(counts) > 0
      if (!any(keep)) return(as.vector(parameters$transition[, , k, h]))
      as.vector(.multilpa_weighted_logits(
        stacked_design[keep, , drop = FALSE], counts[keep, , drop = FALSE],
        matrix(parameters$transition[, , k, h], p))$coefficients)
    }, numeric(p * (n_profiles - 1L)))
  }, numeric(p * (n_profiles - 1L) * n_profiles)),
  c(p, n_profiles - 1L, n_profiles, n_types))
  transition2 <- if (!second) NULL else {
    later_design <- do.call(rbind, designs$transition[-1L])
    array(vapply(seq_len(n_types), function(h) {
      vapply(seq_len(n_profiles^2), function(pair) {
        j <- (pair - 1L) %% n_profiles + 1L
        others <- setdiff(seq_len(n_profiles), j)
        counts <- do.call(rbind, lapply(expectation$moments[[h]]$pairs2, function(moves) {
          matrix(moves[, pair, c(others, j)], nrow(moves))
        }))
        keep <- rowSums(counts) > 0
        if (!any(keep)) return(as.vector(parameters$transition2[, , pair, h]))
        as.vector(.multilpa_weighted_logits(
          later_design[keep, , drop = FALSE], counts[keep, , drop = FALSE],
          matrix(parameters$transition2[, , pair, h], p))$coefficients)
      }, numeric(p * (n_profiles - 1L)))
    }, numeric(p * (n_profiles - 1L) * n_profiles^2)),
    c(p, n_profiles - 1L, n_profiles^2, n_types))
  }
  c(list(measurement = measurement, initial = initial, transition = transition),
    if (second) list(transition2 = transition2),
    list(group_probabilities = colMeans(expectation$group_posteriors)))
}

#' EM for the general transition model from one start
#' @noRd
.lta_em <- function(x, codes, layout, designs, parameters, occasion_of_row,
                    variance_model, min_variance, n_categories, min_probability,
                    max_iter, tol) {
  expectation <- .lta_expectation(x, codes, layout, designs, parameters,
                                  occasion_of_row)
  history <- expectation$log_likelihood
  converged <- FALSE
  iteration <- 0L
  # Sequential by nature: each step starts from the last point.
  while (iteration < max_iter && !converged) {
    updated <- .lta_maximization(x, codes, layout, designs, expectation, parameters,
                                 occasion_of_row, variance_model, min_variance,
                                 n_categories, min_probability)
    updated_expectation <- .lta_expectation(x, codes, layout, designs, updated,
                                            occasion_of_row)
    gain <- updated_expectation$log_likelihood - expectation$log_likelihood
    if (gain < -1e-8 * (1 + abs(expectation$log_likelihood))) {
      stop("EM likelihood decreased beyond numerical roundoff.")
    }
    converged <- abs(gain) <= tol * (1 + abs(expectation$log_likelihood))
    iteration <- iteration + 1L
    history <- c(history, updated_expectation$log_likelihood)
    parameters <- updated
    expectation <- updated_expectation
  }
  list(parameters = parameters, expectation = expectation, converged = converged,
       iterations = iteration, history = history)
}

#' Starting values: the homogeneous initializer, converted to logit
#' coefficients with zero covariate effects
#' @noRd
.lta_initialize <- function(x, group_index, layout, designs, n_profiles, n_types,
                            variance_model, min_variance, start_index, codes,
                            n_categories, min_probability, n_measurement,
                            order = 1L) {
  start <- .multilpa_transition_initialize(x, group_index, n_profiles, n_types,
                                           variance_model, min_variance, start_index,
                                           "diagonal", codes, n_categories,
                                           min_probability)
  block <- start[intersect(c("means", "variances", "response_probabilities"),
                           names(start))]
  p0 <- ncol(designs$initial)
  p <- ncol(designs$transition[[1L]])
  initial <- array(0, c(p0, n_profiles - 1L, n_types))
  transition <- array(0, c(p, n_profiles - 1L, n_profiles, n_types))
  # Intercept columns carry the starting logits; covariate columns start at 0.
  intercepts_initial <- which(colnames(designs$initial) == "(Intercept)")
  intercepts_transition <- which(colnames(designs$transition[[1L]]) == "(Intercept)" |
                                   startsWith(colnames(designs$transition[[1L]]), "occasion_"))
  invisible(lapply(seq_len(n_types), function(h) {
    shares <- start$initial_probabilities[h, ]
    initial[intercepts_initial, , h] <<- log(shares[-n_profiles] / shares[n_profiles])
    tp <- matrix(start$transition_probabilities[, , h], n_profiles)
    lapply(seq_len(n_profiles), function(k) {
      others <- setdiff(seq_len(n_profiles), k)
      logits <- log(tp[k, others] / tp[k, k])
      transition[intercepts_transition, , k, h] <<-
        matrix(logits, length(intercepts_transition), length(others), byrow = TRUE)
    })
  }))
  # Second order starts as the first-order chain: each origin pair (i, j)
  # takes the first-order logits of its current profile j.
  transition2 <- if (order < 2L) NULL else {
    array(vapply(seq_len(n_types), function(h) {
      vapply(seq_len(n_profiles^2), function(pair) {
        as.vector(transition[, , (pair - 1L) %% n_profiles + 1L, h])
      }, numeric(p * (n_profiles - 1L)))
    }, numeric(p * (n_profiles - 1L) * n_profiles^2)),
    c(p, n_profiles - 1L, n_profiles^2, n_types))
  }
  c(list(measurement = rep(list(block), n_measurement), initial = initial,
         transition = transition),
    if (order >= 2L) list(transition2 = transition2),
    list(group_probabilities = start$group_probabilities))
}

#' Fit the general transition model
#'
#' Called by `lta()` when transitions vary by occasion or covariate, the
#' initial distribution depends on covariates, or measurement varies by
#' occasion.
#' @return An object of class `multilpa_lta`.
#' @noRd
.lta_fit_general <- function(data, vars, id, time, n_profiles, n_group_classes,
                             variance_model, n_starts, max_iter, tol, min_variance,
                             seed, categorical, min_probability, occasions,
                             transitions, transition_covariates, initial_covariates,
                             measurement_model, select_start, call, order = 1L) {
  measurement <- .multilpa_prepare_indicators(data, vars, categorical, "error",
                                              min_probability)
  x <- measurement$x
  codes <- measurement$codes
  n_categories <- measurement$n_categories
  time_values <- .multilpa_time_values(data, time, id, vars)
  groups <- .multilpa_prepare_groups(data[[id]])
  layout <- .multilpa_sequence_layout(groups$index, time_values, groups$n, occasions)
  if (layout$n_occasions < 2L) {
    stop(errorCondition(
      "No group is observed at two occasions, so no transition can be estimated.",
      class = "latents_bad_transition", call = NULL))
  }
  occasion_of_row <- integer(nrow(x))
  occasion_of_row[layout$slot[!is.na(layout$slot)]] <- col(layout$slot)[!is.na(layout$slot)]
  designs <- .lta_designs(layout, .lta_covariate_matrix(data, transition_covariates),
                          .lta_covariate_matrix(data, initial_covariates), transitions)
  n_measurement <- if (identical(measurement_model, "occasion")) layout$n_occasions else 1L
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)  # nolint: object_name_linter. R's name for the RNG state.
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  centers <- if (ncol(x) > 0L) colMeans(x) else numeric(0)
  x <- sweep(x, 2L, centers, "-")
  categorical_names <- setdiff(vars, measurement$continuous)
  if (!is.null(n_categories)) names(n_categories) <- categorical_names
  # What the parameter maps and the scores read, before a fit exists.
  spec <- list(n_profiles = as.integer(n_profiles),
               n_group_classes = as.integer(n_group_classes),
               continuous = measurement$continuous, categorical = categorical_names,
               n_categories = n_categories, variance_model = variance_model,
               measurement = vector("list", n_measurement), designs = designs,
               layout = layout, occasion_of_row = occasion_of_row, x = x,
               codes = codes, group_index = groups$index, n_groups = groups$n,
               group_ids = groups$ids, order = as.integer(order))
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch({
      start <- .lta_initialize(x, groups$index, layout, designs, n_profiles,
                               n_group_classes, variance_model, min_variance,
                               start_index, codes, n_categories, min_probability,
                               n_measurement, order)
      # EM to a loose tolerance, then a quasi-Newton finish on the exact
      # likelihood: EM alone crawls when profiles are weakly separated.
      em <- .lta_em(x, codes, layout, designs, start, occasion_of_row, variance_model,
                    min_variance, n_categories, min_probability, max_iter,
                    max(tol, 1e-6))
      if (max_iter == 0L) return(em)
      finish <- .lta_quasi_newton(spec, em$parameters, min_variance, tol)
      expectation <- .lta_expectation(x, codes, layout, designs, finish$parameters,
                                      occasion_of_row)
      list(parameters = finish$parameters, expectation = expectation,
           converged = isTRUE(finish$converged), iterations = em$iterations,
           history = c(em$history, expectation$log_likelihood),
           max_gradient = finish$gradient)
    }, error = function(error) list(error = conditionMessage(error)))
  })
  valid <- vapply(attempts, function(a) is.null(a$error), logical(1))
  if (!any(valid)) {
    stop(errorCondition(sprintf("All %d starts failed: %s", n_starts,
      paste(unique(vapply(attempts, `[[`, character(1), "error")), collapse = "; ")),
      class = "latents_all_starts_failed", call = NULL))
  }
  scores <- vapply(attempts, function(a) if (is.null(a$error))
    a$expectation$log_likelihood else -Inf, numeric(1))
  converged <- vapply(attempts, function(a) isTRUE(a$converged), logical(1))
  best_start <- .multilpa_select_start(scores, converged, select_start)
  best <- attempts[[best_start]]
  .lta_result(best, attempts, scores, valid, converged, best_start, centers, data,
              vars, id, time, groups, layout, designs, x, codes, measurement,
              occasion_of_row, n_profiles, n_group_classes, variance_model,
              transitions, transition_covariates, initial_covariates,
              measurement_model, min_variance, min_probability, call, order)
}

#' Free parameters of the general transition model
#' @noRd
.lta_count_parameters <- function(n_profiles, n_types, d, n_categories,
                                  variance_model, n_measurement, p0, p, order = 1L) {
  per_measurement <- n_profiles * d +
    (if (identical(variance_model, "equal")) 1L else n_profiles) * d +
    if (is.null(n_categories)) 0L else n_profiles * sum(n_categories - 1L)
  as.integer(n_measurement * per_measurement + (n_types - 1L) +
               n_types * p0 * (n_profiles - 1L) +
               n_types * p * (n_profiles - 1L) * n_profiles +
               if (order >= 2L) n_types * p * (n_profiles - 1L) * n_profiles^2 else 0L)
}

#' Assemble a `multilpa_lta` fit
#' @noRd
.lta_result <- function(best, attempts, scores, valid, converged, best_start,
                        centers, data, vars, id, time, groups, layout, designs, x,
                        codes, measurement, occasion_of_row, n_profiles,
                        n_group_classes, variance_model, transitions,
                        transition_covariates, initial_covariates,
                        measurement_model, min_variance, min_probability, call,
                        order = 1L) {
  profile_names <- paste0("profile_", seq_len(n_profiles))
  class_names <- paste0("group_class_", seq_len(n_group_classes))
  parameters <- best$parameters
  parameters$measurement <- lapply(parameters$measurement, function(block) {
    if (!is.null(block$means)) {
      block$means <- sweep(block$means, 2L, centers, "+")
      dimnames(block$means) <- dimnames(block$variances) <-
        list(profile_names, measurement$continuous)
    }
    block
  })
  expectation <- best$expectation
  n_parameters <- .lta_count_parameters(n_profiles, n_group_classes, ncol(x),
                                        measurement$n_categories, variance_model,
                                        length(parameters$measurement),
                                        ncol(designs$initial),
                                        ncol(designs$transition[[1L]]), order)
  log_likelihood <- expectation$log_likelihood
  subject_posteriors <- expectation$subject_posteriors
  group_posteriors <- expectation$group_posteriors
  dimnames(subject_posteriors) <- list(NULL, profile_names)
  dimnames(group_posteriors) <- list(groups$ids, class_names)
  names(parameters$group_probabilities) <- class_names
  dimnames(parameters$initial) <- list(colnames(designs$initial),
                                       profile_names[-n_profiles], class_names)
  # Dimension 2 indexes the destinations other than the origin, in profile
  # order; the transition table names them.
  dimnames(parameters$transition) <- list(colnames(designs$transition[[1L]]),
                                          paste0("destination_", seq_len(n_profiles - 1L)),
                                          profile_names, class_names)
  if (!is.null(parameters$transition2)) {
    pairs <- expand.grid(current = profile_names, previous = profile_names,
                         stringsAsFactors = FALSE)
    dimnames(parameters$transition2) <- list(
      colnames(designs$transition[[1L]]),
      paste0("destination_", seq_len(n_profiles - 1L)),
      paste0(pairs$previous, "->", pairs$current), class_names)
  }
  boundary <- any(vapply(parameters$measurement, function(block) {
    any(block$variances <= min_variance * (1 + 1e-8))
  }, logical(1)))
  result <- list(
    call = call, vars = vars, continuous = measurement$continuous,
    categorical = setdiff(vars, measurement$continuous),
    id = id, time = time, group_ids = groups$ids, group_index = groups$index,
    n_profiles = as.integer(n_profiles), n_group_classes = as.integer(n_group_classes),
    n_groups = groups$n, n_observations = nrow(x), n_occasions = layout$n_occasions,
    transitions = transitions, transition_covariates = transition_covariates,
    initial_covariates = initial_covariates, measurement_model = measurement_model,
    variance_model = variance_model,
    measurement = parameters$measurement, initial_coefficients = parameters$initial,
    transition_coefficients = parameters$transition,
    second_order_coefficients = parameters$transition2, order = as.integer(order),
    group_probabilities = parameters$group_probabilities,
    designs = designs, layout = layout, occasion_of_row = occasion_of_row,
    centers = centers, x = x, codes = codes,
    n_categories = if (is.null(measurement$n_categories)) NULL else
      stats::setNames(measurement$n_categories, setdiff(vars, measurement$continuous)),
    subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
    expectation = expectation,
    log_likelihood = log_likelihood, n_parameters = n_parameters,
    aic = -2 * log_likelihood + 2 * n_parameters,
    bic = -2 * log_likelihood + log(groups$n) * n_parameters,
    bic_individual = -2 * log_likelihood + log(nrow(x)) * n_parameters,
    converged = isTRUE(best$converged), iterations = best$iterations,
    log_likelihood_history = best$history, boundary = boundary,
    starts = data.frame(start = seq_along(attempts),
                        log_likelihood = ifelse(valid, scores, NA_real_),
                        converged = converged,
                        error = vapply(attempts, function(a) a$error %||% NA_character_,
                                       character(1)),
                        selected = seq_along(attempts) == best_start),
    n_best_replicated = sum(valid & abs(scores - log_likelihood) <=
                              1e-6 * (1 + abs(log_likelihood))),
    min_variance = min_variance, min_probability = min_probability)
  class(result) <- "multilpa_lta"
  if (!result$converged) {
    warning(warningCondition("The best start did not converge; increase max_iter.",
                             class = "latents_unconverged", call = NULL))
  }
  if (boundary) {
    warning(warningCondition("A variance reached min_variance; this is a bound-active fit.",
                             class = "latents_boundary", call = NULL))
  }
  result
}
