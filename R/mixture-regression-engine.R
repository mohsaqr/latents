# Estimation engine for finite mixtures of regressions.
#
# One likelihood covers three nestings. Write f_k(y_i) for the density of row i
# under the regression of class k.
#
#   "observation"  every row has its own class:
#                  L = prod_i sum_k pi_ik f_k(y_i)
#   "group"        every row of a group shares one class:
#                  L = prod_j sum_k pi_jk prod_{i in j} f_k(y_i)
#   "two-level"    rows have their own class, and a group class h shifts how
#                  probable each class is for the rows of its groups:
#                  L = prod_j sum_h eta_jh prod_{i in j} sum_k pi_ik|h f_k(y_i)
#
# The class probabilities are multinomial logits of the membership covariates,
# with the first class as the reference. In the two-level model each group
# class has its own intercepts and the covariate slopes are shared.
#
# Everything is estimated by EM. The regression M-step maximizes the
# posterior-weighted complete-data likelihood by Newton-Raphson, which for the
# Gaussian family is a single weighted least-squares solve; the mixing M-step
# is a weighted multinomial logit fitted to the posterior counts.
#
# References: DeSarbo, W. S., & Cron, W. L. (1988). A maximum likelihood
# methodology for clusterwise linear regression. Journal of Classification, 5,
# 249-282. Wedel, M., & DeSarbo, W. S. (1995). A mixture likelihood approach
# for generalized linear models. Journal of Classification, 12, 21-55.
# Vermunt, J. K. (2003). Multilevel latent class models. Sociological
# Methodology, 33, 213-239.

#' Stable log(1 + exp(x))
#' @param x Numeric vector.
#' @return Numeric vector of the same length.
#' @noRd
.mixture_log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

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

#' Linear predictors of every class
#' @param spec The model specification from `.mixture_spec()`.
#' @param params The parameter list.
#' @return An `n x K` matrix.
#' @noRd
.mixture_linear_predictors <- function(spec, params) {
  eta <- spec$x %*% params$beta
  if (ncol(spec$z) > 0L) eta <- eta + drop(spec$z %*% params$common)
  eta + spec$offset
}

#' Log density of every row under every class
#' @param spec The model specification.
#' @param params The parameter list.
#' @return An `n x K` matrix of log densities.
#' @noRd
.mixture_log_density <- function(spec, params) {
  eta <- .mixture_linear_predictors(spec, params)
  y <- spec$y
  switch(spec$family,
    gaussian = {
      sigma2 <- matrix(params$sigma2, nrow(eta), ncol(eta), byrow = TRUE)
      -0.5 * (log(2 * pi * sigma2) + (y - eta)^2 / sigma2)
    },
    binomial = y * eta - spec$trials * .mixture_log1pexp(eta) +
      spec$log_normalizer,
    poisson = y * eta - exp(eta) + spec$log_normalizer)
}

#' E-step: likelihood and posterior class probabilities
#'
#' @param spec The model specification.
#' @param params The parameter list.
#' @return A list with `log_likelihood`, the per-unit `unit_log_likelihood`,
#'   the row posteriors `tau` (`n x K`) and, by nesting, `group_tau`
#'   (`G x K`), or `rho` (`G x H`) and `tau_by_group_class` (a list of `H`
#'   `n x K` joint posteriors). `log_density` is kept for the scores.
#' @noRd
.mixture_expectation <- function(spec, params, weighted = TRUE) {
  expectation <- .mixture_unweighted_expectation(spec, params)
  if (!isTRUE(weighted) || is.null(spec$sampling_weights)) return(expectation)
  .mixture_weigh(spec, expectation)
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

#' The E-step without sampling weights
#' @noRd
.mixture_unweighted_expectation <- function(spec, params) {
  log_density <- .mixture_log_density(spec, params)
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

#' Family pieces for Newton-Raphson on the regression block
#'
#' @param family Family name.
#' @param y Response (successes for the binomial).
#' @param trials Binomial trials.
#' @param eta Linear predictor matrix, `n x K`.
#' @param sigma2 Gaussian variances, one per class.
#' @return A list of `gradient` (d loglik / d eta) and `weight`
#'   (-d2 loglik / d eta2), both `n x K`.
#' @noRd
.mixture_family_derivatives <- function(family, y, trials, eta, sigma2) {
  switch(family,
    gaussian = {
      precision <- matrix(1 / sigma2, nrow(eta), ncol(eta), byrow = TRUE)
      list(gradient = (y - eta) * precision, weight = precision)
    },
    binomial = {
      p <- stats::plogis(eta)
      list(gradient = y - trials * p, weight = trials * p * (1 - p))
    },
    poisson = {
      mu <- exp(eta)
      list(gradient = y - mu, weight = mu)
    })
}

#' Block-diagonal matrix from a list of square blocks
#' @param blocks List of numeric matrices.
#' @return A numeric matrix.
#' @noRd
.mixture_block_diagonal <- function(blocks) {
  sizes <- vapply(blocks, nrow, integer(1))
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  out <- matrix(0, sum(sizes), sum(sizes))
  placed <- Map(function(block, from, to) {
    list(index = seq.int(from, to), block = block)
  }, blocks, starts, ends)
  # Each block writes to its own diagonal window, so the order is immaterial.
  Reduce(function(matrix_so_far, piece) {
    matrix_so_far[piece$index, piece$index] <- piece$block
    matrix_so_far
  }, placed, out)
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

#' Posterior-weighted complete-data log likelihood of the regression block
#' @noRd
.mixture_regression_objective <- function(spec, params, tau) {
  sum(tau * .mixture_log_density(spec, params))
}

#' M-step for the regression coefficients
#'
#' Newton-Raphson on the posterior-weighted complete-data log likelihood,
#' `sum_ik tau_ik log f_k(y_i)`, over the class-specific coefficients and the
#' coefficients shared by every class, holding the Gaussian variances fixed.
#' The objective is concave in the coefficients for all three canonical-link
#' families, and a step is halved until it does not decrease the objective, so
#' the update can only increase the EM likelihood.
#'
#' @param spec The model specification.
#' @param params The current parameter list.
#' @param tau Row posteriors, `n x K`.
#' @param max_newton Maximum Newton iterations.
#' @return The parameter list with `beta` and `common` updated, and
#'   `separation` set when a class's coefficients diverged.
#' @noRd
.mixture_update_regression <- function(spec, params, tau, max_newton = 50L) {
  n_classes <- spec$n_classes
  p <- ncol(spec$x)
  q <- ncol(spec$z)
  pack <- function(beta, common) c(as.vector(beta), common)
  unpack <- function(theta, pars) {
    pars$beta <- matrix(theta[seq_len(p * n_classes)], p, n_classes,
                        dimnames = dimnames(pars$beta))
    if (q > 0L) pars$common <- theta[p * n_classes + seq_len(q)]
    pars
  }
  objective <- function(pars) .mixture_regression_objective(spec, pars, tau)
  newton_step <- function(pars) {
    eta <- .mixture_linear_predictors(spec, pars)
    pieces <- .mixture_family_derivatives(spec$family, spec$y, spec$trials,
                                         eta, pars$sigma2)
    score_weight <- tau * pieces$gradient
    info_weight <- tau * pieces$weight
    gradient <- c(as.vector(crossprod(spec$x, score_weight)),
                  if (q > 0L) rowSums(crossprod(spec$z, score_weight)) else
                    numeric())
    class_blocks <- lapply(seq_len(n_classes), function(k) {
      crossprod(spec$x, spec$x * info_weight[, k])
    })
    information <- .mixture_block_diagonal(class_blocks)
    if (q > 0L) {
      cross <- do.call(rbind, lapply(seq_len(n_classes), function(k) {
        crossprod(spec$x, spec$z * info_weight[, k])
      }))
      common_block <- crossprod(spec$z, spec$z * rowSums(info_weight))
      information <- rbind(cbind(information, cross),
                           cbind(t(cross), common_block))
    }
    # A ridge far below the data scale keeps an empty class solvable without
    # moving a well-determined solution.
    ridge <- 1e-10 * (1 + max(abs(diag(information))))
    .mixture_solve(information + diag(ridge, nrow(information)), gradient)
  }
  current <- objective(params)
  iteration <- 0L
  # Newton iterations are inherently sequential; each needs the previous point.
  repeat {
    iteration <- iteration + 1L
    direction <- newton_step(params)
    theta <- pack(params$beta, params$common)
    step <- 1
    accepted <- FALSE
    while (step > 1e-8 && !accepted) {
      candidate <- unpack(theta + step * direction, params)
      value <- objective(candidate)
      accepted <- is.finite(value) && value >= current - 1e-12 * (1 + abs(current))
      if (!accepted) step <- step / 2
    }
    if (!accepted) break
    gain <- value - current
    params <- candidate
    current <- value
    # The Gaussian objective is quadratic: one full step is the exact optimum.
    if (identical(spec$family, "gaussian") && isTRUE(all.equal(step, 1))) break
    if (gain <= 1e-12 * (1 + abs(current)) || iteration >= max_newton) break
  }
  # Divergence is a property of the logit and log links; a Gaussian
  # coefficient of any size is on the outcome's own scale.
  params$separation <- !identical(spec$family, "gaussian") &&
    (any(abs(params$beta) > 30) || (q > 0L && any(abs(params$common) > 30)))
  params
}

#' M-step for the Gaussian variances
#' @param spec The model specification.
#' @param params The parameter list.
#' @param tau Row posteriors.
#' @return The parameter list with `sigma2` updated.
#' @noRd
.mixture_update_variance <- function(spec, params, tau) {
  if (!identical(spec$family, "gaussian")) return(params)
  residual2 <- (spec$y - .mixture_linear_predictors(spec, params))^2
  weighted <- colSums(tau * residual2)
  sizes <- colSums(tau)
  sigma2 <- if (identical(spec$variance, "equal")) {
    rep(sum(weighted) / sum(sizes), spec$n_classes)
  } else weighted / pmax(sizes, .Machine$double.eps)
  params$sigma2 <- pmax(sigma2, spec$min_variance)
  params
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

#' One full M-step
#' @noRd
.mixture_maximization <- function(spec, params, expectation) {
  tau <- expectation$tau
  params <- .mixture_update_regression(spec, params, tau)
  params <- .mixture_update_variance(spec, params, tau)
  .mixture_update_mixing(spec, params, expectation)
}

#' Run EM from a starting point
#'
#' @param spec The model specification.
#' @param params Starting parameters.
#' @param max_iter Maximum EM iterations.
#' @param tol Relative convergence tolerance on the log likelihood.
#' @return A list with `params`, `expectation`, `converged`, `iterations`,
#'   `history` and `degenerate` (a class emptied or a variance hit the floor).
#' @noRd
.mixture_em <- function(spec, params, max_iter, tol) {
  expectation <- .mixture_expectation(spec, params)
  history <- expectation$log_likelihood
  converged <- FALSE
  iteration <- 0L
  # EM is a fixed-point iteration: each step needs the previous posteriors.
  while (iteration < max_iter && !converged) {
    params <- .mixture_maximization(spec, params, expectation)
    updated <- .mixture_expectation(spec, params)
    improvement <- updated$log_likelihood - expectation$log_likelihood
    if (!is.finite(updated$log_likelihood)) break
    if (improvement < -1e-7 * (1 + abs(expectation$log_likelihood))) {
      stop(errorCondition(sprintf(paste(
        "EM likelihood decreased by %.3g, beyond numerical roundoff. This is a",
        "defect; please report it with the call that produced it."),
        -improvement), class = "latents_em_decrease", call = NULL))
    }
    converged <- abs(improvement) <= tol * (1 + abs(expectation$log_likelihood))
    iteration <- iteration + 1L
    history <- c(history, updated$log_likelihood)
    expectation <- updated
  }
  sizes <- colSums(expectation$tau)
  degenerate <- any(sizes < ncol(spec$x) + 1) ||
    (identical(spec$family, "gaussian") &&
       any(params$sigma2 <= spec$min_variance * (1 + 1e-8)))
  list(params = params, expectation = expectation, converged = converged,
       iterations = iteration, history = history, degenerate = degenerate)
}
