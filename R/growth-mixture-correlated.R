# Growth mixture models with correlated residuals (phase G5 of todo.MD):
# mixture_regression(random =, residual_correlation = "ar" | "brownian").
# Within class k, person i's outcomes have covariance
#   V_ik = Z_i G_k Z_i' + W_i + sigma_k^2 I,
# where W_i is a continuous-time autoregressive process,
#   W[j, l] = s^2 exp(-rho |t_j - t_l|),
# or Brownian motion, W[j, l] = s^2 min(t_j, t_l), shared by the classes:
# lcmm::hlme(cor = AR(time)) and hlme(cor = BM(time)) (Proust-Lima, Philipps
# and Liquet 2017). W is not diagonal, so the q x q identities of the
# uncorrelated engine do not apply; each class's V is factorized per pattern
# of times and random design, which the persons sharing it reuse.

#' Names of the residual-process coordinates
#' @noRd
.growth_correlation_names <- function(spec) {
  switch(spec$residual_correlation %||% "none",
         none = character(),
         ar = c("residual_correlation.log_rate", "residual_correlation.log_sd"),
         brownian = "residual_correlation.log_sd")
}

#' Whether a growth specification has correlated residuals
#' @noRd
.growth_is_correlated <- function(spec) {
  !identical(spec$residual_correlation %||% "none", "none")
}

#' The residual-process covariance at a vector of times
#' @param spec The growth specification.
#' @param correlation The process parameters: `c(rate, sd)` or `c(sd)`.
#' @param times One person's times.
#' @return A square matrix.
#' @noRd
.growth_process_covariance <- function(spec, correlation, times) {
  switch(spec$residual_correlation,
         ar = correlation[["sd"]]^2 * exp(-correlation[["rate"]] * abs(outer(times, times, "-"))),
         brownian = correlation[["sd"]]^2 * outer(times, times, pmin))
}

#' Persons grouped by their times and random design
#'
#' Two persons with the same times and the same rows of the random design have
#' the same marginal covariance in every class.
#'
#' @param spec The growth specification, with `residual_times`.
#' @param rows Each person's rows (a list).
#' @return A list of patterns, each with `members`, `random` (the pattern's
#'   random design) and `times`.
#' @noRd
.growth_dense_patterns <- function(spec, rows) {
  keys <- vapply(rows, function(r) {
    paste(sprintf("%a", c(spec$residual_times[r], spec$random_design[r, ])),
          collapse = "\r")
  }, character(1))
  groups <- split(seq_along(rows), factor(keys, levels = unique(keys)))
  lapply(unname(groups), function(members) {
    r <- rows[[members[1L]]]
    list(members = members,
         random = spec$random_design[r, , drop = FALSE],
         times = spec$residual_times[r])
  })
}

#' Every class's marginal covariance on every pattern
#' @return A list over classes of lists over patterns of matrices.
#' @noRd
.growth_dense_covariances <- function(spec, stats, params) {
  covariances <- .growth_covariances(spec, params)
  lapply(seq_len(spec$n_classes), function(k) {
    lapply(stats$dense, function(pattern) {
      z <- pattern$random
      z %*% covariances[[k]] %*% t(z) +
        .growth_process_covariance(spec, params$correlation, pattern$times) +
        diag(params$sigma2[k], length(pattern$times))
    })
  })
}

#' Class log densities and random-effect posteriors with correlated residuals
#'
#' @return A list over classes with `log_density` (per person), `mean`
#'   (persons by random effects: posterior means of b), `u` (each row's
#'   V^-1 r, laid out by row), `inverse` (per pattern) and `covariance`
#'   (random by random by patterns: posterior covariances).
#' @noRd
.growth_dense_terms <- function(spec, stats, params, marginal = NULL) {
  marginal <- marginal %||% .growth_dense_covariances(spec, stats, params)
  covariances <- .growth_covariances(spec, params)
  design <- stats$raw$design
  y <- stats$raw$y
  rows <- stats$raw$rows
  n_persons <- length(rows)
  n_random <- ncol(spec$random_design)
  lapply(seq_len(spec$n_classes), function(k) {
    residual <- y - as.vector(design %*% .growth_coefficients(params, k))
    log_density <- numeric(n_persons)
    means <- matrix(0, n_persons, n_random)
    u <- numeric(length(y))
    inverses <- vector("list", length(stats$dense))
    posterior_covariance <- array(0, c(n_random, n_random, length(stats$dense)))
    invisible(lapply(seq_along(stats$dense), function(g) {
      pattern <- stats$dense[[g]]
      factor <- tryCatch(chol(marginal[[k]][[g]]), error = function(condition) {
        stop(errorCondition(paste(
          "A marginal covariance is numerically singular (a process variance",
          "far beyond the residual variance):", conditionMessage(condition)),
          class = "latents_singular_covariance", call = NULL))
      })
      inverse <- chol2inv(factor)
      members <- pattern$members
      at <- do.call(rbind, rows[members])
      r <- matrix(residual[at], length(members))
      solved <- r %*% inverse
      log_density[members] <<- -0.5 * (ncol(r) * log(2 * pi) +
                                         2 * sum(log(diag(factor))) + rowSums(r * solved))
      u[at] <<- solved
      gz <- covariances[[k]] %*% t(pattern$random)
      means[members, ] <<- solved %*% t(gz)
      inverses[[g]] <<- inverse
      posterior_covariance[, , g] <<- covariances[[k]] - gz %*% inverse %*% t(gz)
    }))
    list(log_density = log_density, mean = means, u = u, inverse = inverses,
         covariance = posterior_covariance)
  })
}

#' Per-person scores of the growth mixture with correlated residuals
#'
#' Within class k, d log f / d beta = X'V^-1 r and, for any coordinate of the
#' covariance, d log f / d theta = (u' dV u - tr(V^-1 dV)) / 2 with
#' u = V^-1 r. The derivative dV of every class's marginal covariance on
#' every pattern is taken by central differences of the data-free map from
#' the coordinates to V; the class terms are weighted by the person's
#' posterior.
#'
#' @return A units-by-parameters matrix, as `.growth_scores()`.
#' @noRd
.growth_dense_scores <- function(spec, stats, params, expectation) {
  theta <- .growth_pack(spec, params)
  names_theta <- names(theta)
  structure <- .growth_structure(spec)
  membership <- names_theta %in% .growth_structure_names(structure)
  coefficient <- startsWith(names_theta, "coefficient.")
  covariance_index <- which(!membership & !coefficient)
  posterior <- expectation$posterior
  terms <- expectation$dense
  n_persons <- nrow(posterior)
  person <- spec$group_index
  design <- stats$raw$design
  n_specific <- ncol(spec$x)
  n_common <- ncol(spec$z)
  # X'V^-1 r per person and class.
  xu <- lapply(seq_len(spec$n_classes), function(k) {
    rowsum(design * terms[[k]]$u, person, reorder = TRUE)
  })
  beta_scores <- do.call(cbind, lapply(seq_len(spec$n_classes), function(k) {
    xu[[k]][, seq_len(n_specific), drop = FALSE] * posterior[, k]
  }))
  common_scores <- if (n_common > 0L) {
    Reduce(`+`, lapply(seq_len(spec$n_classes), function(k) {
      xu[[k]][, n_specific + seq_len(n_common), drop = FALSE] * posterior[, k]
    }))
  }
  rows <- stats$raw$rows
  covariance_scores <- vapply(covariance_index, function(j) {
    h <- 1e-6 * max(1, abs(theta[[j]]))
    up <- .growth_dense_covariances(spec, stats, .growth_unpack(spec, replace(theta, j, theta[j] + h), params))
    down <- .growth_dense_covariances(spec, stats, .growth_unpack(spec, replace(theta, j, theta[j] - h), params))
    Reduce(`+`, lapply(seq_len(spec$n_classes), function(k) {
      score <- numeric(n_persons)
      invisible(lapply(seq_along(stats$dense), function(g) {
        derivative <- (up[[k]][[g]] - down[[k]][[g]]) / (2 * h)
        members <- stats$dense[[g]]$members
        at <- do.call(rbind, rows[members])
        u <- matrix(terms[[k]]$u[at], length(members))
        score[members] <<- 0.5 * (rowSums((u %*% derivative) * u) -
                                    sum(terms[[k]]$inverse[[g]] * derivative))
      }))
      score * posterior[, k]
    }))
  }, numeric(n_persons)) |> matrix(n_persons)
  # The covariance coordinates sit between the coefficients and the
  # membership logits, in .growth_names() order.
  block <- cbind(beta_scores, common_scores, covariance_scores)
  scores <- if (is.null(spec$cluster)) {
    prior <- exp(.mixture_log_softmax(spec$w, params$gamma))
    membership_scores <- if (spec$n_classes > 1L) {
      do.call(cbind, lapply(seq_len(spec$n_classes)[-1L], function(k) {
        spec$w * (posterior[, k] - prior[, k])
      }))
    }
    cbind(block, membership_scores)
  } else {
    cbind(rowsum(block, spec$cluster_index, reorder = TRUE),
          .mixture_structure_scores(structure, expectation$structure))
  }
  dimnames(scores) <- list(NULL, names_theta)
  scores
}

#' Bounds of the growth coordinates
#'
#' The same floors as EM's M-step (GR2-8): the residual variance, the process
#' variance and each random effect's conditional variance (its squared
#' log-Cholesky diagonal) at `min_variance`, in standard coordinates. The process rate is kept
#' within a factor of 10^4 of the reciprocal median time gap, beyond which an
#' autoregressive process cannot be told from white noise (fast) or from a
#' random intercept (slow).
#'
#' @return A list with `lower` and `upper`, named like the coordinates.
#' @noRd
.growth_correlated_bounds <- function(spec, theta, gap) {
  names_theta <- names(theta)
  floor <- 0.5 * log(spec$min_variance)
  lower <- stats::setNames(rep(-Inf, length(theta)), names_theta)
  upper <- stats::setNames(rep(Inf, length(theta)), names_theta)
  variance_like <- startsWith(names_theta, "log_sigma.") |
    startsWith(names_theta, "random.log_cholesky.") |
    names_theta == "residual_correlation.log_sd"
  lower[variance_like] <- floor
  rate <- names_theta == "residual_correlation.log_rate"
  lower[rate] <- log(1e-4 / gap)
  upper[rate] <- log(1e4 / gap)
  list(lower = lower, upper = upper)
}

#' Fit the correlated-residual model from a fit without the correlation
#'
#' The model without correlated residuals, fitted from every start, places
#' the classes; the full likelihood is then maximized within the bounds of
#' `.growth_correlated_bounds()` by L-BFGS-B with the analytic scores, from
#' that fit with half of each residual variance moved to the process (at two
#' process rates), and finished by Newton steps on the free coordinates. A
#' coordinate left on a bound is a variance at its floor: the fit is a
#' boundary fit, warned about with `latents_random_boundary`, and judged
#' converged on its other coordinates.
#'
#' @param spec The correlated specification.
#' @param uncorrelated The fit without correlated residuals.
#' @param max_iter,tol Iteration bound and relative tolerance.
#' @return The fit, with `spec`, `stats`, `params` and the likelihood updated.
#' @noRd
.growth_correlated_refit <- function(spec, uncorrelated, max_iter, tol) {
  stats <- .growth_statistics(spec)
  base <- uncorrelated$params
  gaps <- unlist(lapply(stats$raw$rows, function(r) diff(sort(spec$residual_times[r]))))
  gap <- if (any(gaps > 0)) stats::median(gaps[gaps > 0]) else 1
  floor <- spec$min_variance
  starts <- lapply(c(0.5, 2) / gap, function(rate) {
    params <- base
    params$sigma2 <- pmax(base$sigma2 / 2, 2 * floor)
    sd <- sqrt(max(mean(base$sigma2) / 2, 2 * floor))
    params$correlation <- if (identical(spec$residual_correlation, "ar"))
      c(rate = rate, sd = sd) else c(sd = sd)
    params
  })
  if (identical(spec$residual_correlation, "brownian")) starts <- starts[1L]
  weights <- spec$sampling_weights %||% 1
  bounds <- .growth_correlated_bounds(spec, .growth_pack(spec, starts[[1L]]), gap)
  finishes <- lapply(starts, function(start) {
    theta <- pmin(pmax(.growth_pack(spec, start), bounds$lower), bounds$upper)
    # V = Z G Z' + W + sigma^2 I is positive definite in exact arithmetic, but
    # a line search can try a process variance so far beyond sigma^2 that it
    # is singular in floating point, or overflows: either is a rejected step.
    objective <- function(v) {
      value <- tryCatch(
        .growth_expectation(spec, stats, .growth_unpack(spec, v, start))$log_likelihood,
        latents_singular_covariance = function(condition) -Inf)
      if (is.finite(value)) -value else .Machine$double.xmax
    }
    gradient <- function(v) {
      current <- .growth_unpack(spec, v, start)
      tryCatch(-colSums(.growth_scores(spec, stats, current) * weights),
               latents_singular_covariance = function(condition) rep(0, length(v)))
    }
    fitted <- stats::optim(theta, objective, gradient, method = "L-BFGS-B",
                           lower = bounds$lower, upper = bounds$upper,
                           control = list(maxit = max(100L, max_iter),
                                          factr = 10, pgtol = 0))
    list(params = .growth_unpack(spec, fitted$par, start), value = -fitted$value,
         iterations = as.integer(fitted$counts[["gradient"]]))
  })
  best <- finishes[[which.max(vapply(finishes, `[[`, numeric(1), "value"))]]
  # L-BFGS-B stops on its relative tolerance a little short of the maximum,
  # as EM does: Newton steps on the analytic score finish it.
  newton <- .growth_newton_finish(spec, stats, best$params, weights, bounds)
  expectation <- .growth_expectation(spec, stats, newton$params)
  if (length(newton$at_bound) > 0L) {
    warning(warningCondition(sprintf(paste(
      "With correlated residuals, %s %s at %s floor (`min_variance`) or bound,",
      "or not identified: the data cannot separate the residual variance, the",
      "process and the random effects there. Their standard errors are not",
      "reported."),
      paste(sprintf("`%s`", newton$at_bound), collapse = ", "),
      if (length(newton$at_bound) == 1L) "is" else "are",
      if (length(newton$at_bound) == 1L) "its" else "their"),
      class = "latents_random_boundary", call = NULL))
  }
  if (!newton$converged) {
    warning(warningCondition(sprintf(paste(
      "The correlated-residual model did not converge (largest score %.3g",
      "standard errors from its root); the estimates are not a maximum."),
      newton$displacement), class = "latents_unconverged", call = NULL))
  }
  fit <- uncorrelated
  fit$spec <- spec
  fit$stats <- stats
  fit$params <- newton$params
  fit$expectation <- expectation
  fit$log_likelihood <- expectation$log_likelihood
  fit$n_parameters <- .growth_count_parameters(spec)
  degenerate <- .growth_degenerate_covariances(spec, newton$params)
  fit$degenerate_classes <- degenerate
  fit$residual_boundary <- newton$at_bound
  fit$random_boundary <- length(degenerate) > 0L || length(newton$at_bound) > 0L
  fit$converged <- newton$converged
  fit$iterations <- uncorrelated$iterations + best$iterations + newton$iterations
  fit$history <- c(uncorrelated$history, fit$log_likelihood)
  fit
}

#' Newton steps on the growth likelihood
#'
#' The last stage of every growth fit, with or without correlated residuals.
#' The information is the central-difference Jacobian of the analytic total
#' score, restricted to the coordinates not on a bound; each step is halved
#' until the likelihood does not decrease, and clipped to the bounds.
#' Converged when every free score times its standard error (the estimate's
#' distance from the root in standard-error units) is below 1e-4, a
#' hundredth of the displacement at which parameter_inference() warns.
#'
#' @return A list with `params`, `iterations`, `converged`, `displacement`
#'   and `at_bound` (the coordinates left on a bound).
#' @noRd
.growth_newton_finish <- function(spec, stats, params, weights, bounds,
                                 max_steps = 25L) {
  log_likelihood <- function(v) {
    tryCatch(.growth_expectation(spec, stats, .growth_unpack(spec, v, params))$log_likelihood,
             latents_singular_covariance = function(condition) -Inf)
  }
  total <- function(v) {
    colSums(.growth_scores(spec, stats, .growth_unpack(spec, v, params)) * weights)
  }
  on_bound <- function(v, score) {
    width <- 1e-8 * (1 + abs(v))
    (v <= bounds$lower + width & score <= 0) | (v >= bounds$upper - width & score >= 0)
  }
  # The Newton direction and displacement on the free coordinates.
  newton_step <- function(v) {
    score <- total(v)
    free <- !on_bound(v, score)
    information <- .inference_score_information(total, v, 1e-5, scheme = "central",
                                                h_rule = "plus1")[free, free, drop = FALSE]
    inverse <- tryCatch(solve(information), error = function(condition) NULL)
    if (is.null(inverse) || any(diag(inverse) <= 0)) {
      return(list(direction = NULL, displacement = Inf, free = free))
    }
    direction <- numeric(length(v))
    direction[free] <- as.vector(inverse %*% score[free])
    list(direction = direction, free = free,
         displacement = max(abs(score[free]) * sqrt(diag(inverse))))
  }
  theta <- .growth_pack(spec, params)
  value <- log_likelihood(theta)
  iterations <- 0L
  step <- newton_step(theta)
  # Each step needs the last one's result, so they run in sequence.
  while (iterations < max_steps && !is.null(step$direction) &&
         step$displacement >= 1e-8) {  # nolint: undesirable_operator.
    accepted <- FALSE
    multiplier <- 1
    while (!accepted && multiplier > 1e-9) {  # nolint: undesirable_operator.
      candidate <- pmin(pmax(theta + multiplier * step$direction, bounds$lower),
                        bounds$upper)
      candidate_value <- log_likelihood(candidate)
      accepted <- is.finite(candidate_value) &&
        candidate_value >= value - 1e-10 * (1 + abs(value))
      multiplier <- multiplier / 2
    }
    if (!accepted) break
    theta <- candidate
    value <- candidate_value
    iterations <- iterations + 1L
    step <- newton_step(theta)
  }
  score <- total(theta)
  at_bound <- names(theta)[on_bound(theta, score)]
  displacement <- step$displacement
  if (displacement >= 1e-4) {
    # A process the data cannot tell from white noise (a fast rate) or from
    # nothing leaves a flat ridge in its two coordinates: the likelihood is
    # at its maximum in every other direction. That is reported as a
    # boundary of the process, not as non-convergence, when the other
    # coordinates meet the criterion with the process held.
    process <- startsWith(names(theta), "residual_correlation.")
    rest <- !process & !(names(theta) %in% at_bound)
    information <- .inference_score_information(total, theta, 1e-5, scheme = "central",
                                                h_rule = "plus1")[rest, rest, drop = FALSE]
    inverse <- tryCatch(solve(information), error = function(condition) NULL)
    held <- if (is.null(inverse) || any(diag(inverse) <= 0)) Inf else
      max(abs(score[rest]) * sqrt(diag(inverse)))
    if (held < 1e-4) {
      at_bound <- union(at_bound, names(theta)[process])
      displacement <- held
    }
  }
  list(params = .growth_unpack(spec, theta, params), iterations = iterations,
       displacement = displacement, converged = displacement < 1e-4,
       at_bound = at_bound)
}

#' Resolve and check the residual time column of a correlated growth model
#'
#' By default the first numeric variable of `random`, or else of the
#' formula's predictors.
#' @return The column name.
#' @noRd
.growth_residual_time <- function(residual_correlation, residual_time, random, formula,
                                  data) {
  bad <- function(message, class = "latents_bad_argument") {
    stop(errorCondition(message, class = class, call = NULL))
  }
  if (is.null(random)) {
    bad("`residual_correlation` applies to a growth mixture: give `random` as well.")
  }
  if (is.null(residual_time)) {
    # As for trajectory plots: the variables of `random`, then the predictors.
    candidates <- unique(c(all.vars(random), all.vars(formula[[3L]])))
    numeric_ones <- candidates[vapply(candidates, function(v) is.numeric(data[[v]]),
                                      logical(1))]
    if (length(numeric_ones) == 0L) {
      bad(paste("`residual_correlation` needs the time of each observation: name it",
                "with `residual_time`."))
    }
    residual_time <- numeric_ones[1L]
  }
  if (!is.character(residual_time) || length(residual_time) != 1L ||
      !residual_time %in% names(data) || !is.numeric(data[[residual_time]])) {
    bad("`residual_time` must name a numeric column of `data`.")
  }
  if (identical(residual_correlation, "brownian") &&
      any(data[[residual_time]] < 0, na.rm = TRUE)) {
    bad(paste("Brownian motion starts at time zero: `residual_time` must not be",
              "negative."), "latents_bad_data")
  }
  residual_time
}
