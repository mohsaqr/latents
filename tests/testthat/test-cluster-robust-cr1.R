# Cluster-robust variances carry the CR1 small-sample factor G / (G - 1), the
# default of sandwich::vcovCL(), at every robust path.

test_that("the sandwich equals sandwich::vcovCL() with its default adjustment", {
  skip_if_not_installed("sandwich")
  set.seed(5)
  groups <- rep(seq_len(12), each = 6)
  x <- cbind(1, stats::rnorm(72))
  y <- drop(x %*% c(1, 0.4)) + rep(stats::rnorm(12), each = 6) + stats::rnorm(72)
  model <- stats::lm(y ~ x[, 2L])
  bread <- solve(crossprod(x))
  scores <- rowsum(x * stats::residuals(model), groups, reorder = FALSE)
  ours <- .inference_sandwich(bread, crossprod(scores), nrow(scores))
  reference <- sandwich::vcovCL(model, cluster = groups, type = "HC0")
  expect_equal(unname(ours), unname(reference), tolerance = 1e-10)
})

test_that("a one-class robust regression reproduces sandwich::vcovCL()", {
  skip_if_not_installed("sandwich")
  set.seed(11)
  data <- data.frame(id = rep(seq_len(25), each = 8), x = stats::rnorm(200))
  data$y <- 1 + 0.5 * data$x + rep(stats::rnorm(25), each = 8) + stats::rnorm(200)
  fit <- mixture_regression(y ~ x, data = data, n_classes = 1, id = "id",
                            vcov_type = "robust", n_starts = 1, seed = 1)
  ours <- as.data.frame(get_results(fit, "coefficients"))$std_error
  reference <- sqrt(diag(sandwich::vcovCL(stats::lm(y ~ x, data), cluster = ~id,
                                          type = "HC0")))
  expect_equal(ours, unname(reference), tolerance = 1e-4)
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
