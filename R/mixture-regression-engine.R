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
    poisson = matrix(stats::dpois(y, exp(eta), log = TRUE), nrow(eta), ncol(eta)),
    negative_binomial = {
      # NB2 with size r = 1 / alpha: log Gamma(y + r) - log Gamma(r)
      # + r log(r / (r + mu)) + y log(mu / (r + mu)) - log y!.
      size <- matrix(1 / params$sigma2, nrow(eta), ncol(eta), byrow = TRUE)
      matrix(stats::dnbinom(y, size = size, mu = exp(eta), log = TRUE),
             nrow(eta), ncol(eta))
    },
    ordinal = .ordinal_log_density(spec, params, eta))
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

#' The E-step without sampling weights: the regression block's log densities
#' under the specification's structure
#' @noRd
.mixture_unweighted_expectation <- function(spec, params) {
  .mixture_structure_expectation(spec, params, .mixture_log_density(spec, params))
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
    },
    negative_binomial = {
      # The observed information, r mu (r + y) / (r + mu)^2, is positive for
      # every count, so the log-link objective is concave.
      mu <- exp(eta)
      size <- matrix(1 / sigma2, nrow(eta), ncol(eta), byrow = TRUE)
      list(gradient = size * (y - mu) / (size + mu),
           weight = size * mu * (size + y) / (size + mu)^2)
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
  # The ordinal family's thresholds enter the density outside the linear
  # predictor, so it has its own Newton step over thresholds and slopes.
  if (identical(spec$family, "ordinal")) {
    return(.ordinal_update(spec, params, tau, max_newton))
  }
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

#' M-step for the Gaussian variances or the negative-binomial dispersions
#'
#' `sigma2` holds each class's dispersion parameter: the residual variance of
#' the Gaussian family, the dispersion alpha of the negative binomial
#' (variance mu + alpha mu^2); binomial and Poisson classes have none.
#'
#' @param spec The model specification.
#' @param params The parameter list.
#' @param tau Row posteriors.
#' @return The parameter list with `sigma2` updated.
#' @noRd
.mixture_update_variance <- function(spec, params, tau) {
  if (identical(spec$family, "negative_binomial")) {
    return(.mixture_update_dispersion(spec, params, tau))
  }
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

#' M-step for the negative-binomial dispersions
#'
#' Modified Newton on the log dispersion (one per class, or one shared), the
#' means held: the curvature's sign is reflected so every step ascends, and a
#' step is halved until the posterior-weighted log likelihood does not fall.
#' Where the overdispersion score at alpha = 0, sum tau ((y - mu)^2 - y) / 2,
#' is not positive the maximum is the Poisson limit (Karush-Kuhn-Tucker), and
#' the dispersion is set to its floor.
#'
#' @param spec The model specification.
#' @param params The parameter list (`sigma2` holds the dispersions).
#' @param tau Row posteriors.
#' @param max_newton Most Newton steps.
#' @return The parameter list with `sigma2` updated.
#' @noRd
.mixture_update_dispersion <- function(spec, params, tau, max_newton = 50L) {
  mu <- exp(.mixture_linear_predictors(spec, params))
  shared <- identical(spec$variance, "equal")
  floor <- log(.latents_min_dispersion)
  objective <- function(alpha) {
    point <- params
    point$sigma2 <- alpha
    .mixture_regression_objective(spec, point, tau)
  }
  log_alpha <- log(pmax(params$sigma2, .latents_min_dispersion))
  if (shared) log_alpha <- log_alpha[1L]
  expand <- function(values) rep_len(exp(values), spec$n_classes)
  current <- objective(expand(log_alpha))
  iteration <- 0L
  # Newton steps refine one estimate in sequence; nothing to vectorize.
  while (iteration < max_newton) {
    derivatives <- .latents_negative_binomial_dispersion_derivatives(
      spec$y, mu, expand(log_alpha))
    gradient <- colSums(tau * derivatives$score)
    curvature <- colSums(tau * derivatives$curvature)
    if (shared) {
      gradient <- sum(gradient)
      curvature <- sum(curvature)
    }
    held <- log_alpha <= floor + 1e-12 & gradient < 0
    if (all(held | abs(gradient) <= 1e-10 * (1 + sum(tau)))) break
    direction <- ifelse(held, 0, gradient / pmax(abs(curvature),
                                                 1e-10 * (1 + abs(gradient))))
    direction <- direction * min(1, 5 / max(abs(direction)))
    step <- 1
    accepted <- FALSE
    while (step > 1e-12 && !accepted) {
      candidate <- pmax(log_alpha + step * direction, floor)
      value <- objective(expand(candidate))
      accepted <- is.finite(value) && value >= current - 1e-12 * (1 + abs(current))
      if (!accepted) step <- step / 2
    }
    if (!accepted) break
    log_alpha <- candidate
    current <- value
    iteration <- iteration + 1L
  }
  # The Poisson limit is a boundary maximum EM approaches only slowly.
  score_at_zero <- colSums(tau * ((spec$y - mu)^2 - spec$y)) / 2
  boundary <- if (shared) sum(score_at_zero) <= 0 else score_at_zero <= 0
  log_alpha[boundary] <- floor
  params$sigma2 <- expand(log_alpha)
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

#' EM settings of the regression and growth mixtures
#'
#' Their M-steps solve inner Newton (or conditional) maximizations, whose
#' rounding allows a looser decrease guard than a closed-form step; a step to
#' a non-finite likelihood (a separating logit) ends the run, and a decrease
#' is a classed defect.
#' @noRd
.mixture_em_settings <- function() {
  .latents_em_settings(decrease_tolerance = 1e-7, decrease_condition = "classed",
                       stop_on_nonfinite = TRUE)
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
  fit <- .latents_em(params, function(point) .mixture_expectation(spec, point),
                     function(point, expectation) {
                       .mixture_maximization(spec, point, expectation)
                     },
                     max_iter, tol, .mixture_em_settings())
  params <- fit$state
  expectation <- fit$expectation
  converged <- fit$converged
  iteration <- fit$iterations
  history <- fit$history
  sizes <- colSums(expectation$tau)
  # A class needs more rows than its own parameters: the slopes and an
  # intercept, or an ordinal class's C - 1 thresholds.
  class_parameters <- ncol(spec$x) +
    if (identical(spec$family, "ordinal")) spec$n_categories - 1L else 1L
  degenerate <- any(sizes < class_parameters) ||
    (identical(spec$family, "gaussian") &&
       any(params$sigma2 <= spec$min_variance * (1 + 1e-8)))
  list(params = params, expectation = expectation, converged = converged,
       iterations = iteration, history = history, degenerate = degenerate)
}
