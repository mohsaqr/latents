# Cluster-robust variances carry the CR1 small-sample factor G / (G - 1), the
# default of sandwich::vcovCL(), at every robust path.

# A hand-built CR1 reference: G / (G - 1) (X'X)^-1 M (X'X)^-1 with M the
# cross-product of per-cluster score sums. The comparison with
# sandwich::vcovCL() itself is in tests/equivalence/.
cr1_reference <- function(x, residuals, groups) {
  scores <- rowsum(x * residuals, groups, reorder = FALSE)
  n_groups <- nrow(scores)
  bread <- solve(crossprod(x))
  n_groups / (n_groups - 1) * bread %*% crossprod(scores) %*% bread
}

test_that("the sandwich is the hand-built CR1 variance", {
  set.seed(5)
  groups <- rep(seq_len(12), each = 6)
  x <- cbind(1, stats::rnorm(72))
  y <- drop(x %*% c(1, 0.4)) + rep(stats::rnorm(12), each = 6) + stats::rnorm(72)
  residuals <- stats::residuals(stats::lm(y ~ x[, 2L]))
  scores <- rowsum(x * residuals, groups, reorder = FALSE)
  ours <- .inference_sandwich(solve(crossprod(x)), crossprod(scores), nrow(scores))
  expect_equal(ours, cr1_reference(x, residuals, groups), tolerance = 1e-12)
})

test_that("a one-class robust regression gives the hand-built CR1 standard errors", {
  skip_on_cran()
  set.seed(11)
  data <- data.frame(id = rep(seq_len(25), each = 8), x = stats::rnorm(200))
  data$y <- 1 + 0.5 * data$x + rep(stats::rnorm(25), each = 8) + stats::rnorm(200)
  fit <- mixture_regression(y ~ x, data = data, n_classes = 1, id = "id",
                            vcov_type = "robust", n_starts = 1, seed = 1)
  ours <- as.data.frame(get_results(fit, "coefficients"))$std_error
  x <- cbind(1, data$x)
  reference <- sqrt(diag(cr1_reference(x, stats::residuals(stats::lm(y ~ x, data)), data$id)))
  expect_equal(ours, reference, tolerance = 1e-4)
})

test_that("the factor scales the variance by G / (G - 1) and needs two units", {
  expect_equal(.inference_cluster_factor(20), 20 / 19)
  bread <- diag(2)
  meat <- crossprod(matrix(c(1, -1, 2, -2, 0.5, -0.5), 3L, 2L))
  expect_equal(.inference_sandwich(bread, meat, 3) / meat, matrix(1.5, 2L, 2L))
  expect_error(.inference_cluster_factor(1))
  expect_error(.inference_sandwich(bread, meat, 1))
})

test_that("the correction vanishes as the number of units grows", {
  factors <- vapply(c(5, 50, 500, 5000), .inference_cluster_factor, numeric(1))
  expect_true(all(diff(factors) < 0))
  expect_lt(factors[[4L]] - 1, 1e-3)
})
