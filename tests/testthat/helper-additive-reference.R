# Independent references for the additive engine. None of these use the
# sufficient-statistic formulas: they build the full covariance of a group's
# stacked raw observations, condition a joint Gaussian, or integrate.

#' Simulate raw ratings from the additive model
#' @return A data frame with `group`, `y1`, ..., one row per observation.
additive_draw <- function(parameters, sizes, seed) {
  # Scope the seed: restore the caller's RNG state on exit.
  old_seed <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  on.exit(if (is.null(old_seed)) {
    rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(seed)
  d <- ncol(parameters$means)
  classes <- sample.int(nrow(parameters$means), length(sizes), replace = TRUE,
                        prob = parameters$weights)
  blocks <- lapply(seq_along(sizes), function(j) {
    intercept <- stats::rnorm(d, parameters$means[classes[j], ],
                              sqrt(parameters$between[classes[j], ]))
    noise <- matrix(stats::rnorm(sizes[j] * d), sizes[j], d) %*%
      diag(sqrt(parameters$within), d)
    ratings <- sweep(noise, 2L, intercept, "+")
    colnames(ratings) <- paste0("y", seq_len(d))
    data.frame(group = sprintf("g%02d", j), ratings)
  })
  do.call(rbind, blocks)
}

#' Log density of one group's stacked observations under one class, by
#' Cholesky factorisation of the full (n d) x (n d) covariance.
additive_dense_density <- function(x, mu, between, within) {
  n <- nrow(x)
  covariance <- kronecker(diag(within, length(within)), diag(n)) +
    kronecker(diag(between, length(between)), matrix(1, n, n))
  residual <- as.vector(x) - rep(mu, each = n)
  factor <- chol(covariance)
  z <- backsolve(factor, residual, transpose = TRUE)
  -0.5 * (length(residual) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(z^2))
}

#' Posterior mean and variance of each indicator's group intercept by
#' conditioning the joint Gaussian of (B_r, y_1r, ..., y_nr).
additive_dense_moments <- function(x, mu, between, within) {
  n <- nrow(x)
  vapply(seq_len(ncol(x)), function(r) {
    covariance <- diag(within[r], n) + matrix(between[r], n, n)
    gain <- solve(covariance, rep(between[r], n))
    c(mean = mu[r] + sum(gain * (x[, r] - mu[r])),
      variance = between[r] - between[r] * sum(gain))
  }, numeric(2))
}

#' Split a rating frame into one matrix per group, in first-occurrence order.
additive_blocks <- function(data, vars, id) {
  groups <- unique(data[[id]])
  lapply(groups, function(g) as.matrix(subset(data, data[[id]] == g)[, vars]))
}
