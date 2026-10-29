# Cluster-robust variances against sandwich::vcovCL(), whose default
# adjustment is the CR1 factor G / (G - 1).

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


# Moved from tests/testthat/test-cluster-robust-cr1.R: comparisons with independent references.

# A hand-built CR1 reference: G / (G - 1) (X'X)^-1 M (X'X)^-1 with M the
# cross-product of per-cluster score sums. The comparison with
# sandwich::vcovCL() itself is at the top of this file.
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
