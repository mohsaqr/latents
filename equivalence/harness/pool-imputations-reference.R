# From-scratch references for pool_imputations(), in base R.
#
# Two things are written here without calling the package.
#
# 1. One imputation's fit. A single-level Gaussian latent profile model with
#    diagonal, class-varying variances:
#
#      log L = sum_i log sum_k pi_k(x_i) prod_j N(y_ij; mu_kj, s2_kj),
#
#    with pi_k constant (free pi_1, ..., pi_{K-1}; pi_K = 1 - sum) or a
#    multinomial logit against the last profile, pi_k(x) proportional to
#    exp(a_k + x'b_k), a_K = b_K = 0. Its gradient is written analytically;
#    the maximum is reached by Newton steps from a supplied start, with the
#    Hessian a Richardson-extrapolated Jacobian of that gradient (numDeriv),
#    and the observed-information covariance is the inverse of minus that
#    Hessian. Parameters are on the natural scale the package reports (means,
#    variances, probabilities, logit coefficients), so at the maximum the
#    covariance is the one any reparameterisation's delta method gives.
#
# 2. The relabelling of one imputation onto another. Profiles are matched by
#    the permutation that minimises the squared distance between profile means
#    (searched exhaustively over all K! permutations), and every reported
#    parameter is carried through it by one linear map: a probability, mean or
#    variance moves with its profile; a logit against the last profile is
#    re-expressed against the new last profile, b_k' = b_perm(k) - b_perm(K).
#
# Rubin's rules themselves are not written here: they are mice's
# pool.scalar(), called by the suite.

#' The parameters of the model, free and reported, in one fixed order
#'
#' Keys follow the package's documented reporting (`?pool_imputations`):
#' `level`, `outcome`, `term`, `parameter`.
#'
#' @param vars Indicator names.
#' @param n_profiles Number of profiles K.
#' @param covariates Names of the membership covariates; none for constant
#'   class probabilities.
#' @return A list: `reported`, a `data.frame` of keys with `block`
#'   (`"mean"`, `"variance"`, `"probability"`, `"coefficient"`, `"group"`),
#'   `profile` (integer) and `free` (logical); and `terms`, the coefficient
#'   terms.
pool_reference_layout <- function(vars, n_profiles, covariates = character()) {
  stopifnot("`vars` must name at least one indicator" = is.character(vars) && length(vars) > 0L,
            "`n_profiles` must be at least two" = length(n_profiles) == 1L && n_profiles >= 2L)
  profiles <- seq_len(n_profiles)
  measurement <- function(parameter) {
    data.frame(level = "measurement", outcome = rep(sprintf("profile_%d", profiles),
                                                    each = length(vars)),
               term = rep(vars, times = n_profiles), parameter = parameter,
               block = parameter, profile = rep(profiles, each = length(vars)),
               free = TRUE, stringsAsFactors = FALSE)
  }
  terms <- c("(Intercept)", covariates)
  membership <- if (length(covariates) == 0L) {
    rbind(
      data.frame(level = "profile", outcome = sprintf("profile_%d", profiles),
                 term = "group_class_1", parameter = "probability", block = "probability",
                 profile = profiles, free = profiles < n_profiles, stringsAsFactors = FALSE),
      # One group class: its probability is the constant one.
      data.frame(level = "group", outcome = "group_class_1", term = NA_character_,
                 parameter = "probability", block = "group", profile = NA_integer_,
                 free = FALSE, stringsAsFactors = FALSE))
  } else {
    data.frame(level = "profile",
               outcome = rep(sprintf("profile_%d", profiles[-n_profiles]), each = length(terms)),
               term = rep(terms, times = n_profiles - 1L), parameter = "coefficient",
               block = "coefficient", profile = rep(profiles[-n_profiles], each = length(terms)),
               free = TRUE, stringsAsFactors = FALSE)
  }
  reported <- rbind(measurement("mean"), measurement("variance"), membership)
  row.names(reported) <- NULL
  list(reported = reported, terms = terms, vars = vars, n_profiles = n_profiles,
       covariates = covariates)
}

#' Unpack a free-parameter vector
#' @param theta Free parameters in the layout's order.
#' @param layout From [pool_reference_layout()].
#' @return A list with `mean` and `variance` (K x p) and `coefficients`
#'   ((1 + q) x (K - 1)) or `probability` (length K).
.pool_reference_unpack <- function(theta, layout) {
  free <- layout$reported[layout$reported$free, ]
  stopifnot("`theta` must have one value per free parameter" = length(theta) == nrow(free))
  k <- layout$n_profiles
  p <- length(layout$vars)
  block <- function(name) theta[free$block == name]
  out <- list(mean = matrix(block("mean"), k, p, byrow = TRUE),
              variance = matrix(block("variance"), k, p, byrow = TRUE))
  if (length(layout$covariates) == 0L) {
    out$probability <- c(block("probability"), 1 - sum(block("probability")))
  } else {
    out$coefficients <- matrix(block("coefficient"), length(layout$terms), k - 1L)
  }
  out
}

#' Log prior class probabilities and per-class log densities
#' @param parts From [.pool_reference_unpack()].
#' @param y Rows-by-indicators numeric matrix.
#' @param design Rows-by-terms design matrix (intercept first), or `NULL`.
#' @return A list with `log_prior`, `log_density` (rows x K) and `prior`.
.pool_reference_terms <- function(parts, y, design) {
  n <- nrow(y)
  k <- nrow(parts$mean)
  log_density <- vapply(seq_len(k), function(profile) {
    rowSums(stats::dnorm(y, mean = rep(parts$mean[profile, ], each = n),
                         sd = rep(sqrt(parts$variance[profile, ]), each = n), log = TRUE))
  }, numeric(n))
  log_prior <- if (is.null(parts$coefficients)) {
    matrix(log(parts$probability), n, k, byrow = TRUE)
  } else {
    eta <- cbind(design %*% parts$coefficients, 0)
    eta - .pool_log_sum_exp(eta)
  }
  list(log_prior = log_prior, log_density = log_density, prior = exp(log_prior))
}

#' Row-wise log-sum-exp
#' @param a A numeric matrix.
#' @return A numeric vector, one value per row.
.pool_log_sum_exp <- function(a) {
  top <- apply(a, 1L, max)
  top + log(rowSums(exp(a - top)))
}

#' The observed-data log-likelihood
#' @param theta Free parameters in the layout's order.
#' @param layout From [pool_reference_layout()].
#' @param y,design As for [.pool_reference_terms()].
#' @return One number.
pool_reference_log_likelihood <- function(theta, layout, y, design = NULL) {
  pieces <- .pool_reference_terms(.pool_reference_unpack(theta, layout), y, design)
  sum(.pool_log_sum_exp(pieces$log_prior + pieces$log_density))
}

#' The analytic gradient of [pool_reference_log_likelihood()]
#'
#' With posteriors w_ik: d/d mu_kj = sum_i w_ik (y_ij - mu_kj) / s2_kj;
#' d/d s2_kj = sum_i w_ik ((y_ij - mu_kj)^2 - s2_kj) / (2 s2_kj^2);
#' d/d pi_k = sum_i (w_ik / pi_k - w_iK / pi_K); d/d b_k = X'(w_k - pi_k(x)).
#' @inheritParams pool_reference_log_likelihood
#' @return A numeric vector in the layout's free order.
pool_reference_gradient <- function(theta, layout, y, design = NULL) {
  parts <- .pool_reference_unpack(theta, layout)
  pieces <- .pool_reference_terms(parts, y, design)
  joint <- pieces$log_prior + pieces$log_density
  w <- exp(joint - .pool_log_sum_exp(joint))
  k <- layout$n_profiles
  per_profile <- lapply(seq_len(k), function(profile) {
    residual <- sweep(y, 2L, parts$mean[profile, ])
    v <- parts$variance[profile, ]
    list(mean = colSums(w[, profile] * residual) / v,
         variance = colSums(w[, profile] * sweep(residual^2, 2L, v)) / (2 * v^2))
  })
  membership <- if (is.null(parts$coefficients)) {
    colSums(w[, -k, drop = FALSE]) / parts$probability[-k] - sum(w[, k]) / parts$probability[k]
  } else {
    as.vector(crossprod(design, w[, -k, drop = FALSE] - pieces$prior[, -k, drop = FALSE]))
  }
  c(unlist(lapply(per_profile, `[[`, "mean")), unlist(lapply(per_profile, `[[`, "variance")),
    membership)
}

#' Carry a start to the maximum by Newton steps, and invert the information
#'
#' The analytic gradient is first checked against numDeriv's numerical one
#' near the start, so the maximum it defines is the likelihood's. Steps stop when
#' the largest relative step falls below `1e-11`; failing to get there in
#' `max_steps` is an error, never a silently unconverged reference.
#'
#' @param start Free parameters to start from.
#' @param layout From [pool_reference_layout()].
#' @param y,design As for [.pool_reference_terms()].
#' @param max_steps Bound on the Newton steps.
#' @return A list: `theta`, `log_likelihood`, `covariance` (free x free),
#'   `gradient` (at the maximum) and `steps`.
pool_reference_maximise <- function(start, layout, y, design = NULL, max_steps = 30L) {
  gradient <- function(theta) pool_reference_gradient(theta, layout, y, design)
  # Checked 2% away from the start, where the gradient is large: at a maximum
  # both are zero and their gap is numDeriv's rounding (about 1e-6 here).
  away <- 0.98 * start
  numerical <- numDeriv::grad(pool_reference_log_likelihood, away, layout = layout, y = y,
                              design = design)
  stopifnot("the analytic gradient must agree with the numerical one" =
              max(abs(numerical - gradient(away))) / max(abs(numerical), 1) < 1e-7)
  hessian <- function(theta) {
    h <- numDeriv::jacobian(gradient, theta)
    (h + t(h)) / 2
  }
  # Recursion rather than a loop: Newton converges quadratically from a
  # converged start, so the depth is a handful of steps.
  step_from <- function(theta, steps) {
    step <- solve(hessian(theta), gradient(theta))
    theta <- theta - step
    if (max(abs(step) / pmax(abs(theta), 1)) < 1e-11) return(list(theta = theta, steps = steps))
    if (steps >= max_steps) {
      stop(errorCondition(sprintf("Newton refinement did not converge in %d steps", max_steps),
                          class = "multilpa_bad_reference", call = NULL))
    }
    step_from(theta, steps + 1L)
  }
  reached <- step_from(start, 1L)
  information <- -hessian(reached$theta)
  list(theta = reached$theta,
       log_likelihood = pool_reference_log_likelihood(reached$theta, layout, y, design),
       covariance = solve(information), gradient = gradient(reached$theta),
       steps = reached$steps)
}

#' Every reported parameter, its covariance and standard error
#'
#' Derived rows: pi_K = 1 - sum of the free probabilities (Jacobian -1 on
#' each), and the single group class's probability, the constant one.
#'
#' @param theta,covariance A maximum and its covariance, free order.
#' @param layout From [pool_reference_layout()].
#' @return A list: `estimate` and `covariance` over `layout$reported`.
pool_reference_reported <- function(theta, covariance, layout) {
  reported <- layout$reported
  free_rows <- which(reported$free)
  jacobian <- matrix(0, nrow(reported), length(free_rows))
  jacobian[cbind(free_rows, seq_along(free_rows))] <- 1
  estimate <- numeric(nrow(reported))
  estimate[free_rows] <- theta
  last <- which(reported$block == "probability" & !reported$free)
  if (length(last) == 1L) {
    probability_columns <- which(reported$block[free_rows] == "probability")
    jacobian[last, probability_columns] <- -1
    estimate[last] <- 1 - sum(theta[probability_columns])
  }
  estimate[reported$block == "group"] <- 1
  list(estimate = estimate, covariance = jacobian %*% covariance %*% t(jacobian))
}

#' The profile permutation matching one fit's means to another's
#' @param means K x p means to relabel.
#' @param target K x p means to match.
#' @return An integer vector `perm`: new profile k is old profile `perm[k]`.
pool_reference_permutation <- function(means, target) {
  k <- nrow(means)
  permutations <- function(v) {
    if (length(v) <= 1L) return(list(v))
    do.call(c, lapply(seq_along(v), function(i) lapply(permutations(v[-i]), function(rest) c(v[i], rest))))
  }
  candidates <- permutations(seq_len(k))
  distance <- vapply(candidates, function(perm) sum((means[perm, , drop = FALSE] - target)^2),
                     numeric(1))
  candidates[[which.min(distance)]]
}

#' The linear map carrying reported parameters through a relabelling
#'
#' Row r of the result gives relabelled parameter r as a combination of the
#' original reported parameters.
#'
#' @param layout From [pool_reference_layout()].
#' @param perm From [pool_reference_permutation()].
#' @return A reported x reported matrix.
pool_reference_relabel <- function(layout, perm) {
  reported <- layout$reported
  k <- layout$n_profiles
  key <- function(profile, row) {
    paste(reported$block[row], sprintf("profile_%d", profile), reported$term[row], sep = "\r")
  }
  own <- paste(reported$block, reported$outcome, reported$term, sep = "\r")
  n <- nrow(reported)
  unit <- function(column) replace(numeric(n), column, 1)
  rows <- vapply(seq_len(n), function(row) {
    profile <- reported$profile[row]
    if (is.na(profile)) return(unit(row))
    if (!identical(reported$block[row], "coefficient")) {
      return(unit(match(key(perm[profile], row), own)))
    }
    # b_k' = b_perm(k) - b_perm(K), with b_K = 0 not a parameter.
    from <- if (perm[profile] != k) unit(match(key(perm[profile], row), own)) else numeric(n)
    base <- if (perm[k] != k) unit(match(key(perm[k], row), own)) else numeric(n)
    from - base
  }, numeric(n))
  map <- t(rows)
  stopifnot("every relabelled parameter must come from reported ones" =
              !anyNA(map) && all(rowSums(abs(map)) > 0))
  map
}

#' The name vcov() gives a parameter
#'
#' `level.parameter.outcome.term`, with a missing term left out, as
#' `?pool_imputations`'s vcov() method names its rows.
#' @param keys A `data.frame` with `level`, `outcome`, `term`, `parameter`.
#' @return A character vector.
pool_reference_names <- function(keys) {
  ifelse(is.na(keys$term),
         paste(keys$level, keys$parameter, keys$outcome, sep = "."),
         paste(keys$level, keys$parameter, keys$outcome, keys$term, sep = "."))
}
