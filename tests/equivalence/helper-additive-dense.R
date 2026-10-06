# Dense reference densities for the additive families (equivalence only).

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
