# Comparisons with independent references moved from tests/testthat/test-audit-extension-repairs.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("Gaussian regression MLE is stable under extreme predictor units", {
  set.seed(423)
  data <- data.frame(x = rnorm(100))
  data$y <- 1 + 2 * data$x + rnorm(100)
  data$x <- data$x * 1e10
  fit <- mixture_regression(y ~ x, data, 1, n_starts = 0, vcov_type = "none")
  reference <- stats::glm(y ~ x, data, family = stats::gaussian())
  expect_true(fit$converged)
  expect_equal(fit$log_likelihood, as.numeric(logLik(reference)), tolerance = 1e-10)
  expect_equal(as.vector(fit$params$beta), unname(coef(reference)), tolerance = 1e-9)
  expect_equal(predict(fit)$fitted, unname(fitted(reference)), tolerance = 1e-8)
})

test_that("binomial separation warning does not depend on predictor units", {
  set.seed(423)
  data <- data.frame(x = rnorm(200))
  data$y <- rbinom(200, 1, plogis(.3 + .8 * data$x))
  ordinary <- mixture_regression(y ~ x, data, 1, family = "binomial",
                                 n_starts = 0, vcov_type = "none")
  data$x <- data$x * .01
  expect_no_warning(scaled <- mixture_regression(y ~ x, data, 1,
    family = "binomial", n_starts = 0, vcov_type = "none"))
  reference <- stats::glm(y ~ x, data, family = stats::binomial())
  expect_gt(abs(scaled$params$beta[2L, 1L]), 30)
  expect_equal(scaled$log_likelihood, as.numeric(logLik(reference)), tolerance = 1e-9)
  expect_equal(predict(scaled)$fitted, predict(ordinary)$fitted, tolerance = 1e-7)
  expect_false(scaled$params$separation)
  data$weight <- 1
  padded <- rbind(data, data.frame(x = 1000, y = 0, weight = 0))
  expect_no_warning(excluded <- mixture_regression(y ~ x, padded, 1,
    family = "binomial", weights = "weight", n_starts = 0, vcov_type = "none"))
  expect_false(excluded$params$separation)
})

test_that("Gaussian regression information is stable around a large predictor origin", {
  set.seed(9001)
  data <- data.frame(x = rnorm(200))
  data$y <- 2 + .4 * data$x + rnorm(200, sd = .7)
  ordinary <- mixture_regression(y ~ x, data, 1, n_starts = 0, vcov_type = "none")
  invisible(lapply(c(1e8, 1e10, 1e12), function(origin) {
    moved_data <- transform(data, x = x + origin)
    moved <- mixture_regression(y ~ x, moved_data, 1, n_starts = 0, vcov_type = "none")
    # Evaluate the independent GLM in a centered chart too; raw QR can itself
    # mistake a large-origin full-rank design for a constant predictor.
    moved_data$centered <- moved_data$x - mean(moved_data$x)
    reference <- glm(y ~ centered, moved_data, family = gaussian())
    expect_true(moved$converged)
    expect_equal(moved$log_likelihood, as.numeric(logLik(reference)), tolerance = 1e-8)
    expect_equal(moved$params$beta[2L, 1L], coef(reference)[[2L]], tolerance = 1e-10)
    # The original intercept is stored as a double on its original scale.
    expect_lt(max(abs(predict(moved)$fitted - fitted(reference))),
              max(1e-8, .Machine$double.eps * origin))
  }))
  expect_true(ordinary$converged)
})
