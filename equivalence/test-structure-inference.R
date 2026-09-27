# Wald standard errors of all fourteen covariance structures against the
# independent reference in tests/testthat/helper-structure-reference.R: its own
# likelihood, its own chart of each structure, and no analytic derivative
# anywhere. The Hessian is the reference's central second difference at
# h = 2e-4 (about eps^(1/4), the step that balances truncation against
# rounding for a log likelihood of order 1e3); the delta-method Jacobian is
# numDeriv's Richardson extrapolation. numDeriv's own `hessian()` is not used
# for the curvature: at its default steps it loses about 1e-4 relative
# precision on these likelihoods (tried and measured on the VVV fit, where
# the package's chart is the long-validated log-Cholesky one).
source(file.path("..", "tests", "testthat", "helper-structure-reference.R"),
       local = TRUE)

structure_reference_data <- function() {
  set.seed(5)
  n_groups <- 50L
  size <- 6L
  unit <- rep(seq_len(n_groups), each = size)
  kind <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size)
  high <- stats::rbinom(length(unit), 1L, ifelse(kind == 1L, 0.8, 0.2))
  a <- stats::rnorm(length(unit), 2.5 * high, 1 + 0.5 * high)
  b <- stats::rnorm(length(unit), -1.5 * high, 1.2 - 0.4 * high) +
    0.4 * a * (1 - high)
  c <- stats::rnorm(length(unit), high, 0.8)
  data.frame(unit = unit, a = a, b = b, c = c)
}

structure_reference_check <- function(n_group_classes, id) {
  data <- structure_reference_data()
  numderiv_jacobian <- function(f, x) numDeriv::jacobian(f, x)
  worst <- vapply(latents:::.multilpa_structures(), function(code) {
    fit <- suppressWarnings(do.call(multilpa, c(
      list(data = data, vars = c("a", "b", "c"), id = id, n_profiles = 2L,
           n_group_classes = n_group_classes, n_starts = 4L, seed = 2L,
           tol = 1e-14, max_iter = 20000L),
      latents:::.multilpa_structure_arguments(code))))
    expect_true(fit$converged, label = code)
    ours <- parameter_inference(fit)
    reference <- reference_structure_se(fit, data, jacobian = numderiv_jacobian)
    expect_equal(reference$log_likelihood, fit$log_likelihood, tolerance = 1e-10,
                 info = code)
    informative <- reference$standard_error > 1e-8
    relative <- abs(ours$standard_error - reference$standard_error)[informative] /
      reference$standard_error[informative]
    expect_lt(max(relative), 1e-5, label = code)
    max(relative)
  }, numeric(1))
  worst
}

test_that("single-level Wald SEs of all fourteen structures match numDeriv", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  worst <- structure_reference_check(1L, NULL)
  message("single-level max relative SE difference by structure:\n",
          paste(sprintf("%s %.2e", names(worst), worst), collapse = "\n"))
})

test_that("two-level Wald SEs of all fourteen structures match numDeriv", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  worst <- structure_reference_check(2L, "unit")
  message("two-level max relative SE difference by structure:\n",
          paste(sprintf("%s %.2e", names(worst), worst), collapse = "\n"))
})
