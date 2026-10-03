# Latent structures: how classes are distributed over the units whose log
# densities a measurement block supplies.
#
# A block (the regression of an outcome, a person's growth trajectory with
# its random effects integrated out) turns parameters into a units-by-classes
# matrix of log densities. A structure turns that matrix into the likelihood
# and the posteriors, and owns the class-membership parameters:
#
#   "observation"  every unit has its own class, with multinomial membership
#                  logits of unit covariates (a row of a regression mixture,
#                  a person of a growth mixture);
#   "group"        the rows of a unit share one class, so their log densities
#                  add before the class prior applies (latent class growth
#                  analysis: the rows of one person);
#   "two-level"    units have their own class, and a group class of their
#                  cluster shifts how probable each class is (Vermunt 2003).
#
# A structure is a list naming its `nesting`, the number of classes
# `n_classes`, the unit membership design `w` (and `intercept_only`), and
# for "group" and "two-level" the unit-to-cluster map `group_index`,
# `n_groups`, `n` and, for "two-level", `n_group_classes`, the cluster
# design `v` and `group_intercept_only`. Sampling weights, one per
# independent unit, are `sampling_weights`, with `row_weights` the weight of
# each row of the density matrix. A regression specification is its own
# structure; a growth model builds one over persons.
#
# Reference: Vermunt, J. K. (2003). Multilevel latent class models.
# Sociological Methodology, 33, 213-239.

#' Row-wise log-sum-exp of a matrix
#' @param log_values Numeric matrix.
#' @return Numeric vector, one value per row.
#' @noRd
.mixture_row_lse <- function(log_values) {
  # pmax() over the columns is the row maximum without a per-row R call.
  row_max <- do.call(pmax, lapply(seq_len(ncol(log_values)), function(k) {
    log_values[, k]
  }))
  finite <- is.finite(row_max)
  if (all(finite)) return(row_max + log(rowSums(exp(log_values - row_max))))
  out <- row_max
  out[finite] <- row_max[finite] +
    log(rowSums(exp(log_values[finite, , drop = FALSE] - row_max[finite])))
  out
}

#' Log class probabilities from a full logit coefficient matrix
#'
#' @param design Numeric design matrix, one row per unit.
#' @param coefficients Matrix with one column per class; the first column is
#'   the reference and is zero.
#' @return Matrix of log probabilities, one row per unit and column per class.
#' @noRd
.mixture_log_softmax <- function(design, coefficients) {
  scores <- design %*% coefficients
  scores - .mixture_row_lse(scores)
}

#' Weight a mixture-regression expectation by sampling weights
#'
#' Every posterior becomes its unit's weight times the posterior, so the
#' M-steps and the scores read weighted counts, and the log likelihood is the
#' weighted pseudo log likelihood. Integer weights equal duplication.
#'
#' @param spec The specification, carrying `sampling_weights` (one per
#'   independent unit) and `row_weights` (each row's unit weight).
#' @param expectation An unweighted expectation.
#' @return The weighted expectation.
#' @noRd
.mixture_weigh <- function(spec, expectation) {
  unit_weights <- spec$sampling_weights
  row_weights <- spec$row_weights
  # Classes at the row level are weighted per row; the units the likelihood
  # sums over are rows, or clusters when `id` names them.
  expectation$log_likelihood <- if (identical(spec$nesting, "observation")) {
    sum(row_weights * expectation$unit_log_likelihood)
  } else sum(unit_weights * expectation$unit_log_likelihood)
  expectation$tau <- expectation$tau * row_weights
  if (!is.null(expectation$group_tau)) {
    expectation$group_tau <- expectation$group_tau * unit_weights
  }
  if (!is.null(expectation$rho)) {
    expectation$rho <- expectation$rho * unit_weights
    expectation$tau_by_group_class <- lapply(expectation$tau_by_group_class,
                                             `*`, row_weights)
  }
  expectation
}

#' The E-step of a structure, without sampling weights
#'
#' @param spec The structure (a regression specification is one).
#' @param params The parameter list, carrying the structure's membership
#'   coefficients: `gamma`, or `delta` and `class_logits` for two levels.
#' @param log_density The block's log densities: units by classes for
#'   `"observation"` and `"two-level"`, rows by classes for `"group"`, whose
#'   rows the structure sums to units.
#' @return A list with `log_likelihood`, the per-unit `unit_log_likelihood`,
#'   the posteriors `tau` (one row per row of `log_density`), `log_prior`,
#'   `log_density` and, by nesting, `group_tau`, or `rho`,
#'   `tau_by_group_class`, `log_prior_by_group_class` and `log_eta`.
#' @noRd
.mixture_structure_expectation <- function(spec, params, log_density) {
  switch(spec$nesting,
    observation = {
      log_prior <- .mixture_log_softmax(spec$w, params$gamma)
      joint <- log_prior + log_density
      unit <- .mixture_row_lse(joint)
      list(log_likelihood = sum(unit), unit_log_likelihood = unit,
           tau = exp(joint - unit), log_prior = log_prior,
           log_density = log_density)
    },
    group = {
      log_prior <- .mixture_log_softmax(spec$w, params$gamma)
      joint <- log_prior + rowsum(log_density, spec$group_index,
                                  reorder = TRUE)
      unit <- .mixture_row_lse(joint)
      group_tau <- exp(joint - unit)
      list(log_likelihood = sum(unit), unit_log_likelihood = unit,
           tau = group_tau[spec$group_index, , drop = FALSE],
           group_tau = group_tau, log_prior = log_prior,
           log_density = log_density)
    },
    "two-level" = {
      log_eta <- .mixture_log_softmax(spec$v, params$delta)
      n_group_classes <- ncol(log_eta)
      # Per group class h: the row's marginal density given h, and the joint
      # class posterior given h.
      within <- lapply(seq_len(n_group_classes), function(h) {
        log_prior <- .mixture_log_softmax(
          .mixture_two_level_design(spec, h), params$class_logits)
        joint <- log_prior + log_density
        marginal <- .mixture_row_lse(joint)
        list(log_prior = log_prior, joint = joint, marginal = marginal)
      })
      group_log <- log_eta + vapply(within, function(piece) {
        drop(rowsum(piece$marginal, spec$group_index, reorder = TRUE))
      }, numeric(spec$n_groups))
      unit <- .mixture_row_lse(group_log)
      rho <- exp(group_log - unit)
      tau_by_group_class <- lapply(seq_len(n_group_classes), function(h) {
        rho[spec$group_index, h] *
          exp(within[[h]]$joint - within[[h]]$marginal)
      })
      list(log_likelihood = sum(unit), unit_log_likelihood = unit,
           tau = Reduce(`+`, tau_by_group_class), rho = rho,
           tau_by_group_class = tau_by_group_class,
           log_prior_by_group_class = lapply(within, `[[`, "log_prior"),
           log_eta = log_eta, log_density = log_density)
    })
}

#' Class-logit design of the two-level model for one group class
#'
#' Group class `h` contributes an indicator column for its own intercepts,
#' followed by the shared membership covariates.
#'
#' @param spec The model specification.
#' @param h Group class index.
#' @return An `n x (H + r)` matrix.
#' @noRd
.mixture_two_level_design <- function(spec, h) {
  indicator <- matrix(0, spec$n, spec$n_group_classes)
  indicator[, h] <- 1
  cbind(indicator, spec$w)
}

#' Solve a symmetric positive definite system, falling back to QR
#'
#' @param a Numeric matrix.
#' @param b Numeric vector.
#' @return The solution vector.
#' @noRd
.mixture_solve <- function(a, b) {
  factor <- tryCatch(chol(a), error = function(e) NULL)
  if (!is.null(factor)) return(backsolve(factor, forwardsolve(t(factor), b)))
  qr.solve(a, b, tol = 1e-12)
}

#' Weighted multinomial logit on soft counts
#'
#' Maximizes `sum_i sum_k counts_ik log p_ik(design_i)` over a coefficient
#' matrix whose first column is fixed at zero, by Newton-Raphson with step
#' halving. The counts need not sum to one per row: the two-level mixing step
#' passes joint posteriors whose row sums are the group-class posteriors.
#'
#' @param design Numeric design matrix.
#' @param counts Non-negative matrix, one column per outcome category.
#' @param start Full coefficient matrix (first column zero).
#' @param max_newton Maximum Newton iterations.
#' @return A list with the full `coefficients` matrix and `separation`.
#' @noRd
.mixture_multinomial <- function(design, counts, start, max_newton = 100L) {
  n_categories <- ncol(counts)
  r <- ncol(design)
  if (n_categories == 1L || r == 0L) {
    return(list(coefficients = start, separation = FALSE))
  }
  totals <- rowSums(counts)
  free <- seq.int(2L, n_categories)
  objective <- function(coefficients) {
    sum(counts * .mixture_log_softmax(design, coefficients))
  }
  newton_direction <- function(coefficients) {
    p <- exp(.mixture_log_softmax(design, coefficients))
    gradient <- as.vector(crossprod(design, counts[, free, drop = FALSE] -
                                      totals * p[, free, drop = FALSE]))
    index <- expand.grid(a = seq_along(free), b = seq_along(free))
    blocks <- lapply(seq_len(nrow(index)), function(cell) {
      a <- free[index$a[cell]]
      b <- free[index$b[cell]]
      w <- totals * p[, a] * ((a == b) - p[, b])
      crossprod(design, design * w)
    })
    # Assemble the (K-1) x (K-1) grid of r x r blocks.
    grid_rows <- lapply(seq_along(free), function(a) {
      do.call(cbind, blocks[index$a == a])
    })
    information <- do.call(rbind, grid_rows)
    ridge <- 1e-10 * (1 + max(abs(diag(information))))
    step <- .mixture_solve(information + diag(ridge, nrow(information)),
                          gradient)
    matrix(step, r, length(free))
  }
  coefficients <- start
  current <- objective(coefficients)
  iteration <- 0L
  # Sequential by nature: every Newton step starts from the previous point.
  repeat {
    iteration <- iteration + 1L
    direction <- newton_direction(coefficients)
    step <- 1
    accepted <- FALSE
    while (step > 1e-8 && !accepted) {
      candidate <- coefficients
      candidate[, free] <- coefficients[, free] + step * direction
      value <- objective(candidate)
      accepted <- is.finite(value) &&
        value >= current - 1e-12 * (1 + abs(current))
      if (!accepted) step <- step / 2
    }
    if (!accepted) break
    gain <- value - current
    coefficients <- candidate
    current <- value
    if (gain <= 1e-12 * (1 + abs(current)) || iteration >= max_newton) break
  }
  list(coefficients = coefficients, separation = any(abs(coefficients) > 30))
}

#' M-step for the mixing probabilities
#'
#' Without membership covariates the maximizer is closed form; the Newton
#' routine reaches it too, but the closed form is exact and cheaper.
#'
#' @param spec The model specification.
#' @param params The parameter list.
#' @param expectation The E-step result.
#' @return The parameter list with the mixing coefficients updated.
#' @noRd
.mixture_update_mixing <- function(spec, params, expectation) {
  closed_form <- function(counts) {
    shares <- pmax(colSums(counts), .Machine$double.xmin)
    matrix(log(shares) - log(shares[1L]), nrow = 1L)
  }
  switch(spec$nesting,
    observation = ,
    group = {
      counts <- if (identical(spec$nesting, "group")) expectation$group_tau else
        expectation$tau
      if (spec$intercept_only) {
        params$gamma <- closed_form(counts)
      } else {
        fit <- .mixture_multinomial(spec$w, counts, params$gamma)
        params$gamma <- fit$coefficients
        params$membership_separation <- fit$separation
      }
    },
    "two-level" = {
      if (spec$group_intercept_only) {
        params$delta <- closed_form(expectation$rho)
      } else {
        fit <- .mixture_multinomial(spec$v, expectation$rho, params$delta)
        params$delta <- fit$coefficients
        params$group_separation <- fit$separation
      }
      if (ncol(spec$w) == 0L) {
        # Intercepts only: each group class's class shares are closed form.
        params$class_logits <- do.call(rbind, lapply(
          expectation$tau_by_group_class, closed_form))
      } else {
        design <- do.call(rbind, lapply(seq_len(spec$n_group_classes),
                                        function(h) {
          .mixture_two_level_design(spec, h)
        }))
        counts <- do.call(rbind, expectation$tau_by_group_class)
        fit <- .mixture_multinomial(design, counts, params$class_logits)
        params$class_logits <- fit$coefficients
        params$membership_separation <- fit$separation
      }
    })
  params
}

#' Scores of a structure's membership coefficients
#'
#' By Fisher's identity a unit's score for a multinomial logit is its
#' design row times the posterior minus the prior class probabilities
#' (posterior counts minus their expectation). The two-level structure has
#' cluster-level scores for the group-class logits and unit-level scores for
#' the class logits, summed to clusters.
#'
#' @param spec The structure.
#' @param expectation Its weighted E-step.
#' @return A matrix with one row per independent unit (row, or cluster) and
#'   one column per free membership coefficient, in packing order.
#' @noRd
.mixture_structure_scores <- function(spec, expectation) {
  tau <- expectation$tau
  to_units <- function(scores) {
    if (identical(spec$nesting, "observation")) scores else
      rowsum(scores, spec$group_index, reorder = TRUE)
  }
  mixing <- switch(spec$nesting,
    observation = {
      residual <- tau - rowSums(tau) * exp(expectation$log_prior)
      do.call(cbind, lapply(seq.int(2L, length.out = spec$n_classes - 1L),
                            function(k) spec$w * residual[, k]))
    },
    group = {
      residual <- expectation$group_tau -
        rowSums(expectation$group_tau) * exp(expectation$log_prior)
      do.call(cbind, lapply(seq.int(2L, length.out = spec$n_classes - 1L),
                            function(k) spec$w * residual[, k]))
    },
    "two-level" = {
      rho <- expectation$rho
      group_residual <- rho - rowSums(rho) * exp(expectation$log_eta)
      delta_scores <- do.call(cbind, lapply(
        seq.int(2L, length.out = spec$n_group_classes - 1L),
        function(h) spec$v * group_residual[, h]))
      # Row-level class-logit scores, summed to groups below.
      residual_by_group_class <- lapply(
        seq_len(spec$n_group_classes), function(h) {
          expectation$tau_by_group_class[[h]] -
            rho[spec$group_index, h] *
            exp(expectation$log_prior_by_group_class[[h]])
        })
      logit_scores <- do.call(cbind, lapply(
        seq.int(2L, length.out = spec$n_classes - 1L), function(k) {
          intercepts <- vapply(residual_by_group_class, function(r) r[, k],
                               numeric(spec$n))
          intercepts <- matrix(intercepts, spec$n)
          slopes <- spec$w * Reduce(`+`, lapply(residual_by_group_class,
                                                function(r) r[, k]))
          cbind(intercepts, slopes)
        }))
      cbind(delta_scores %||% matrix(0, spec$n_groups, 0L),
            to_units(logit_scores %||% matrix(0, spec$n, 0L)))
    })
  # One class (and one group class) has no free membership coefficients.
  mixing %||% matrix(0, if (identical(spec$nesting, "observation"))
    spec$n else spec$n_groups, 0L)
}


#' Forward-backward pass of the Markov structure over occasions
#'
#' The structure of the transition models: a profile chain per group, with
#' group- and occasion-specific transition log probabilities (homogeneous
#' transitions repeat one matrix over groups and occasions). Occasions a
#' group was not observed at hold its state, so sequences of unequal length
#' share one layout.
#'
#' @param emission Offset-removed emission list (one groups x K matrix per
#'   occasion).
#' @param layout Sequence layout.
#' @param log_initial Groups x K log initial probabilities.
#' @param log_transition List over occasions 2..T of groups x K x K arrays.
#' @return A list with `alpha`, `beta` and `log_scaled`.
#' @noRd
.latents_forward_backward <- function(emission, layout, log_initial, log_transition) {
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


#' Log-sum-exp over the last dimension that tolerates rows of -Inf
#' @noRd
.latents_log_sum_exp <- function(values) {
  top <- .multilpa_row_max(values)
  finite <- is.finite(top)
  out <- rep(-Inf, nrow(values))
  if (any(finite)) {
    shifted <- exp(values[finite, , drop = FALSE] - top[finite])
    out[finite] <- top[finite] + log(rowSums(shifted))
  }
  out
}


#' The nested structure of the profile models
#'
#' Each group belongs to one group class, and each of its rows to one
#' profile, with profile probabilities that depend on the group class:
#' L = prod_j sum_h omega_h prod_{i in j} sum_k pi_k|h f_k(y_i) (Vermunt
#' 2003). One group class is the ordinary single-level mixture. The full
#' cross-level family adds a group-level density per group class.
#'
#' @param log_density Rows-by-profiles log densities, offset removed.
#' @param density_offset The per-row offset that was removed.
#' @param group_index Group of each row.
#' @param parameters A list with `profile_probabilities` (group classes by
#'   profiles), `group_probabilities` and optionally `group_log_density`
#'   (groups by group classes).
#' @return A list with `log_likelihood`, `group_log_likelihood`,
#'   `group_posteriors`, `subject_posteriors` and `joint` (the row posteriors
#'   joint with each group class).
#' @noRd
.multilpa_nested_structure <- function(log_density, density_offset, group_index,
                                       parameters) {
  n_types <- length(parameters$group_probabilities)
  # Each group type supplies a different prior over the same Gaussian profiles.
  conditional <- lapply(seq_len(n_types), function(group_type) {
    log_joint <- sweep(log_density, 2L,
                       log(parameters$profile_probabilities[group_type, ]), "+")
    log_marginal <- .multilpa_log_sum_exp(log_joint)
    posterior <- exp(sweep(log_joint, 1L, .multilpa_row_max(log_joint), "-"))
    list(log_marginal = log_marginal,
         posterior = posterior / rowSums(posterior))
  })
  group_scores <- matrix(vapply(seq_len(n_types), function(group_type) {
    as.numeric(rowsum(conditional[[group_type]]$log_marginal,
                      group_index, reorder = FALSE))
  }, numeric(max(group_index))), nrow = max(group_index), ncol = n_types)
  # The full cross-level family adds a group-level measurement density per
  # group class (manifest group means as between indicators); absent, nothing
  # changes.
  if (!is.null(parameters$group_log_density)) {
    stopifnot("group_log_density must be groups x group classes" =
                identical(dim(parameters$group_log_density), dim(group_scores)))
    group_scores <- group_scores + parameters$group_log_density
  }
  group_offset <- .multilpa_row_max(group_scores)
  group_scores <- sweep(sweep(group_scores, 1L, group_offset, "-"), 2L,
                        log(parameters$group_probabilities), "+")
  group_log_likelihood <- .multilpa_log_sum_exp(group_scores) + group_offset +
    as.numeric(rowsum(density_offset, group_index, reorder = FALSE))
  group_posteriors <- exp(sweep(group_scores, 1L, .multilpa_row_max(group_scores), "-"))
  group_posteriors <- group_posteriors / rowSums(group_posteriors)
  joint <- lapply(seq_len(n_types), function(group_type) {
    conditional[[group_type]]$posterior * group_posteriors[group_index, group_type]
  })
  subject_posteriors <- Reduce(`+`, joint)
  log_likelihood <- sum(group_log_likelihood)
  if (!is.finite(log_likelihood) || any(!is.finite(subject_posteriors))) {
    stop("Non-finite likelihood or posterior probabilities.")
  }
  list(log_likelihood = log_likelihood, group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors, subject_posteriors = subject_posteriors,
       joint = joint)
}


#' The nested structure with membership covariates
#'
#' As `.multilpa_nested_structure()`, with the profile and group-class
#' probabilities multinomial logits of covariates: the profile logits of
#' group class `h` read the design `profile_design[[h]]` (its own intercepts,
#' shared or class-specific slopes), and the group-class logits the group
#' design. The last class is the reference.
#'
#' @param log_density Rows-by-profiles log densities, offset removed.
#' @param density_offset The per-row offset that was removed.
#' @param group_index Group of each row.
#' @param profile_design List of profile design matrices, one per group class.
#' @param group_design Group-level design matrix.
#' @param beta Profile logit coefficients.
#' @param gamma Group logit coefficients.
#' @return As `.multilpa_nested_structure()`, plus `group_priors` and
#'   `profile_priors` (one rows-by-profiles matrix per group class).
#' @noRd
.multilpa_logit_structure <- function(log_density, density_offset, group_index,
                                      profile_design, group_design, beta, gamma) {
  conditional <- lapply(profile_design, function(design) {
    log_prior <- .multilpa_log_softmax(design, beta)
    prior <- exp(log_prior)
    scores <- log_density + log_prior
    offset <- .multilpa_row_max(scores)
    weights <- exp(sweep(scores, 1L, offset, "-"))
    total <- rowSums(weights)
    marginal <- offset + log(total)
    list(prior = prior, marginal = marginal,
         posterior = weights / total)
  })
  group_log_prior <- .multilpa_log_softmax(group_design, gamma)
  group_prior <- exp(group_log_prior)
  evidence <- matrix(vapply(seq_along(conditional), function(h) {
    as.vector(rowsum(conditional[[h]]$marginal, group_index, reorder = FALSE))
  }, numeric(nrow(group_design))), nrow(group_design))
  evidence_offset <- .multilpa_row_max(evidence)
  scores <- sweep(evidence, 1L, evidence_offset, "-") + group_log_prior
  score_offset <- .multilpa_row_max(scores)
  weights <- exp(sweep(scores, 1L, score_offset, "-"))
  total <- rowSums(weights)
  group_posteriors <- weights / total
  group_log_likelihood <- evidence_offset + score_offset + log(total) +
    as.vector(rowsum(density_offset, group_index, reorder = FALSE))
  joint <- lapply(seq_along(conditional), function(h) {
    conditional[[h]]$posterior * group_posteriors[group_index, h]
  })
  list(log_likelihood = sum(group_log_likelihood),
       group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors,
       subject_posteriors = Reduce(`+`, joint),
       joint = joint, group_priors = group_prior,
       profile_priors = lapply(conditional, `[[`, "prior"))
}


#' Weight an expectation step by sampling weights
#'
#' Pseudo maximum likelihood maximizes the weighted sum of the independent
#' units' log likelihoods. Its EM step is the ordinary one with every
#' posterior multiplied by its unit's weight, so integer weights reproduce the
#' fit to the data with each unit repeated that many times.
#'
#' @param expectation An unweighted expectation step.
#' @param weights Sampling weight per group.
#' @param group_index Group of each row.
#' @return The expectation with weighted posteriors, the weighted log
#'   likelihood, and `n_total`, the weighted row count.
#' @noRd
.multilpa_weigh <- function(expectation, weights, group_index) {
  row_weights <- weights[group_index]
  expectation$log_likelihood <- sum(weights * expectation$group_log_likelihood)
  expectation$group_posteriors <- expectation$group_posteriors * weights
  expectation$subject_posteriors <- expectation$subject_posteriors * row_weights
  expectation$joint <- lapply(expectation$joint, `*`, row_weights)
  expectation$n_total <- sum(row_weights)
  expectation
}


#' Pack a structure's free membership coefficients
#'
#' The reference class (and reference group class) has its logits fixed at
#' zero; the rest are named `membership.<class>.<term>` and, for two levels,
#' `group_membership.<group class>.<term>`, with the class logits of a
#' two-level structure carrying one intercept per group class
#' (`(Intercept):group_class_h`) and shared slopes.
#'
#' @param spec The structure.
#' @param params The parameter list.
#' @return A named numeric vector (empty for one class).
#' @noRd
.mixture_structure_pack <- function(spec, params) {
  # paste0() recycles a zero-length argument to one string, so an empty block
  # (one class, no group-class covariates) must be named explicitly as empty.
  labelled <- function(values, labels) {
    if (length(values) == 0L) numeric() else stats::setNames(values, labels)
  }
  class_names <- paste0("class_", seq_len(spec$n_classes))
  free_classes <- class_names[-1L]
  mixing <- if (identical(spec$nesting, "two-level")) {
    group_names <- paste0("group_class_", seq_len(spec$n_group_classes))
    delta <- params$delta[, -1L, drop = FALSE]
    logits <- params$class_logits[, -1L, drop = FALSE]
    logit_rows <- c(paste0("(Intercept):", group_names), colnames(spec$w))
    c(labelled(as.vector(delta), paste0(
      "group_membership.", rep(group_names[-1L], each = ncol(spec$v)), ".",
      rep(colnames(spec$v), spec$n_group_classes - 1L))),
      labelled(as.vector(logits), paste0(
        "membership.", rep(free_classes, each = length(logit_rows)), ".",
        rep(logit_rows, length(free_classes)))))
  } else {
    gamma <- params$gamma[, -1L, drop = FALSE]
    labelled(as.vector(gamma), paste0(
      "membership.", rep(free_classes, each = ncol(spec$w)), ".",
      rep(colnames(spec$w), length(free_classes))))
  }
  mixing
}

#' Unpack a structure's free membership coefficients
#'
#' @param spec The structure.
#' @param values The packed coefficients, in `.mixture_structure_pack()` order.
#' @param template A parameter list with the right shapes.
#' @return `template` with `gamma`, or `delta` and `class_logits`, replaced.
#' @noRd
.mixture_structure_unpack <- function(spec, values, template) {
  k <- spec$n_classes
  position <- 0L
  take <- function(count) {
    taken <- values[position + seq_len(count)]
    position <<- position + count
    unname(taken)
  }
  # Filled in place, so the template's shapes and names are kept.
  params <- template
  if (identical(spec$nesting, "two-level")) {
    params$delta[] <- cbind(0, matrix(take(ncol(spec$v) * (spec$n_group_classes - 1L)),
                                      ncol(spec$v)))
    params$class_logits[] <- cbind(0, matrix(take(nrow(template$class_logits) * (k - 1L)),
                                             nrow(template$class_logits)))
  } else {
    params$gamma[] <- cbind(0, matrix(take(ncol(spec$w) * (k - 1L)), ncol(spec$w)))
  }
  stopifnot("the membership coefficients must be consumed exactly" =
              position == length(values))
  params
}
