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
  design <- spec$x
  common_design <- spec$z
  coefficients <- params$beta
  intercept <- match("(Intercept)", colnames(design), nomatch = 0L)
  if (intercept > 0L && nrow(design) > 0L &&
      (any(abs(design) > 1e6) || any(abs(common_design) > 1e6))) {
    center <- function(matrix) {
      if (ncol(matrix) == 0L) return(numeric())
      means <- colMeans(matrix)
      spread <- apply(abs(sweep(matrix, 2L, means, "-")), 2L, max)
      means[abs(means) <= 1e6 | abs(means) <= 1e4 * spread] <- 0
      means
    }
    x_center <- center(design)
    z_center <- center(common_design)
    if (any(c(x_center, z_center) != 0)) {
      # Avoid a different large-intercept cancellation on every row. Any
      # unavoidable rounding of the raw intercept remains one common shift.
      coefficients[intercept, ] <- coefficients[intercept, ] +
        as.vector(crossprod(x_center, coefficients)) + sum(z_center * params$common)
      design <- sweep(design, 2L, x_center, "-")
      if (ncol(common_design) > 0L) common_design <- sweep(common_design, 2L, z_center, "-")
    }
  }
  eta <- design %*% coefficients
  if (ncol(common_design) > 0L) eta <- eta + drop(common_design %*% params$common)
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
    negative_binomial = .mixture_negative_binomial_log_density(
      y, exp(eta), matrix(params$sigma2, nrow(eta), ncol(eta), byrow = TRUE)),
    ordinal = .ordinal_log_density(spec, params, eta))
}

#' Negative-binomial (NB2) log density, accurate near the Poisson limit
#'
#' With size r = 1 / alpha, log Gamma(y + r) - log Gamma(r) + r log(r / (r +
#' mu)) + y log(mu / (r + mu)) - log y!. Where alpha max(y, mu) < 1e-3 that
#' difference of large gamma terms carries rounding noise of about 1e-9 per
#' row (1e-7 over a few hundred rows at alpha = 1e-8), larger than the change
#' the dispersion step must resolve there; the Poisson log density plus the
#' fourth-order expansion of the NB-to-Poisson log ratio in alpha (whose
#' derivatives `.latents_negative_binomial_dispersion_derivatives()` uses) is
#' exact to about 4e-16 instead (review 2026-10-26, round 2, N4).
#'
#' @param y Counts, one per row.
#' @param mu Means, rows by classes.
#' @param alpha Dispersions, rows by classes.
#' @return A rows-by-classes matrix of log densities.
#' @noRd
.mixture_negative_binomial_log_density <- function(y, mu, alpha) {
  count <- matrix(y, nrow(mu), ncol(mu))
  out <- matrix(stats::dnbinom(count, size = 1 / alpha, mu = mu, log = TRUE),
                nrow(mu), ncol(mu))
  small <- alpha * pmax(count, mu) < 1e-3
  if (any(small)) {
    k <- count[small]
    m <- mu[small]
    a <- alpha[small]
    sums <- list(k * (k - 1) / 2, k * (k - 1) * (2 * k - 1) / 6,
                 (k * (k - 1) / 2)^2,
                 k * (k - 1) * (2 * k - 1) * (3 * k^2 - 3 * k - 1) / 30)
    ratio <- Reduce(`+`, lapply(seq_len(4L), function(order) {
      (-1)^(order + 1) * a^order *
        ((sums[[order]] - k * m^order) / order + m^(order + 1) / (order + 1))
    }))
    out[small] <- stats::dpois(k, m, log = TRUE) + ratio
  }
  out
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
  intercept <- match("(Intercept)", colnames(spec$x), nomatch = 0L)
  x_center <- numeric(p)
  z_center <- numeric(q)
  if (intercept > 0L) {
    supported <- rowSums(tau) > 0
    center <- function(design) {
      if (ncol(design) == 0L) return(numeric())
      rows <- design[supported, , drop = FALSE]
      means <- colMeans(rows)
      spread <- apply(abs(sweep(rows, 2L, means, "-")), 2L, max)
      # Constancy is read on the raw values: the rounded mean of an exactly
      # constant column (three 0.1s) leaves a spread of one ulp, and centring
      # would turn an aliased constant into a column of rounding noise.
      constant <- apply(rows, 2L, function(column) {
        length(column) == 0L || max(column) == min(column)
      })
      means[constant | abs(means) <= 1e4 * spread] <- 0
      means
    }
    x_center <- center(spec$x)
    z_center <- center(spec$z)
    # An intercept absorbs these shifts exactly. Construct cross-products in
    # this chart before they can lose the small variance around a large origin.
    params$beta[intercept, ] <- params$beta[intercept, ] +
      as.vector(crossprod(x_center, params$beta)) + sum(z_center * params$common)
    spec$x <- sweep(spec$x, 2L, x_center, "-")
    if (q > 0L) spec$z <- sweep(spec$z, 2L, z_center, "-")
  }
  if (identical(spec$family, "gaussian") && q == 0L && any(x_center != 0)) {
    # Absolute weighted QR avoids repeatedly rounding a very large raw
    # intercept in successive Newton increments. The class variance is a
    # common scalar in each independent least-squares problem and cancels.
    # An empty or rank-deficient class has no unique least-squares solution;
    # it is left to the ridge-stabilized Newton step below.
    systems <- lapply(seq_len(n_classes), function(class) {
      largest <- max(tau[, class])
      if (!(largest > 0)) return(NULL)
      weight <- sqrt(tau[, class] / largest)
      design <- spec$x * weight
      scale <- sqrt(colSums(design^2))
      if (!all(scale > 0)) return(NULL)
      decomposition <- qr(sweep(design, 2L, scale, "/"), tol = 1e-12)
      if (decomposition$rank < p) return(NULL)
      list(decomposition = decomposition, weight = weight, scale = scale)
    })
    if (!any(vapply(systems, is.null, logical(1)))) {
      coefficients <- lapply(systems, function(system) {
        coefficient <- qr.coef(system$decomposition,
                               (spec$y - spec$offset) * system$weight) / system$scale
        coefficient[intercept] <- coefficient[intercept] - sum(x_center * coefficient)
        coefficient
      })
      params$beta[,] <- do.call(cbind, coefficients)
      params$separation <- FALSE
      return(params)
    }
  }
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
    # Solve in information-standardized coordinates when predictor units
    # differ greatly. A ridge based on the largest raw diagonal can otherwise
    # overwhelm the intercept (and turn an ordinary Gaussian MLE into a ridge
    # estimate). Keep the ordinary-scale path's arithmetic unchanged.
    diagonal <- diag(information)
    ridge <- 1e-10 * (1 + max(abs(diagonal)))
    if (any(ridge > 1e-7 * pmax(diagonal, .Machine$double.xmin))) {
      scale <- sqrt(pmax(diagonal, .Machine$double.xmin))
      standardized <- information / tcrossprod(scale)
      direction <- tryCatch(.mixture_solve(standardized, gradient / scale),
                            error = function(condition) {
        .mixture_solve(standardized + diag(1e-10, nrow(information)),
                       gradient / scale)
      })
      return(as.vector(direction / scale))
    }
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
  # A raw slope changes with predictor units. Detect extreme fitted linear
  # predictors, which are invariant to that change of coordinates.
  supported <- (spec$row_weights %||% rep(1, spec$n)) > 0
  params$separation <- !identical(spec$family, "gaussian") &&
    any(abs(.mixture_linear_predictors(spec, params)[supported, , drop = FALSE]) > 30)
  if (intercept > 0L) {
    params$beta[intercept, ] <- params$beta[intercept, ] -
      as.vector(crossprod(x_center, params$beta)) - sum(z_center * params$common)
  }
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
#' is not positive the Poisson limit is a local maximum (Karush-Kuhn-Tucker),
#' and the dispersion is set to its floor unless the Newton point is higher.
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
  # The Poisson limit is a boundary maximum EM approaches only slowly. A
  # non-positive score at zero makes it a local maximum only: the profile in
  # alpha can rise again to a higher interior maximum, so a class moves to the
  # floor only when that does not lower its own objective.
  score_at_zero <- colSums(tau * ((spec$y - mu)^2 - spec$y)) / 2
  boundary <- if (shared) sum(score_at_zero) <= 0 else score_at_zero <= 0
  if (any(boundary)) {
    class_objective <- function(values) {
      point <- params
      point$sigma2 <- expand(values)
      by_class <- colSums(tau * .mixture_log_density(spec, point))
      if (shared) sum(by_class) else by_class
    }
    at_floor <- log_alpha
    at_floor[boundary] <- floor
    reached <- class_objective(log_alpha)
    gain <- class_objective(at_floor) - reached
    boundary <- boundary & gain >= -1e-12 * (1 + abs(reached))
  }
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
.mixture_em_settings <- function(stop_scale = NULL) {
  .latents_em_settings(decrease_tolerance = 1e-7, decrease_condition = "classed",
                       stop_on_nonfinite = TRUE, stop_scale = stop_scale)
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
                     max_iter, tol,
                     # A unit-free stop (review 2026-10-26 round 2, K2-6): one
                     # outcome per row of positive weight.
                     .mixture_em_settings(stop_scale = if (is.null(spec$row_weights))
                       length(spec$y) else sum(spec$row_weights > 0)))
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
