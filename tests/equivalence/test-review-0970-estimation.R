# Comparisons with independent references moved from tests/testthat/test-review-0970-estimation.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("direct effects fit correlated full-rank covariates by weighted QR", {
  set.seed(10)
  data <- data.frame(z = stats::rnorm(200), w = stats::rnorm(200))
  data$w <- data$z + 1e-4 * data$w
  data$y <- 1 + 2 * data$z - 3 * data$w + stats::rnorm(200)
  fit <- lpa(data, "y", 1, direct_effects = "y ~ z + w", n_starts = 1,
             search = "kmeans", tol = 1e-10)
  reference <- stats::coef(stats::lm(y ~ z + w, data))
  estimates <- c(as.vector(fit$means), unname(fit$direct_effects))
  expect_lt(max(abs(estimates - unname(reference))), 1e-6)
})
