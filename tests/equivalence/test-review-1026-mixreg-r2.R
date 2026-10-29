# Comparisons with independent references moved from tests/testthat/test-review-1026-mixreg-r2.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("standard errors do not depend on a predictor's origin (N6)", {
  set.seed(3)
  n <- 300
  d <- data.frame(t = stats::rnorm(n))
  d$y <- 1 + 2 * d$t + stats::rnorm(n)
  d$count <- stats::rpois(n, exp(0.5 + 0.3 * d$t))
  for (shift in c(1e4, 5e4, 1e6)) {
    shifted <- d
    shifted$x <- d$t + shift
    gaussian <- mixture_regression(y ~ x, shifted, 1, n_starts = 0)
    reference <- stats::glm(y ~ x, data = shifted)
    expect_equal(as.data.frame(gaussian)$std_error[2],
                 sqrt(stats::vcov(reference)[2, 2] * (n - 2) / n), tolerance = 1e-7)
    poisson <- mixture_regression(count ~ x, shifted, 1, family = "poisson", n_starts = 0)
    reference <- stats::glm(count ~ x, data = shifted, family = stats::poisson)
    expect_equal(as.data.frame(poisson)$std_error[2],
                 sqrt(stats::vcov(reference)[2, 2]), tolerance = 1e-6)
  }
  # Two classes: the slope's errors are the same at origin 0 and 1e5.
  set.seed(2)
  cls <- stats::rbinom(n, 1, 0.4) + 1
  d$z <- ifelse(cls == 1, 1 + 2 * d$t, -1 - d$t) + stats::rnorm(n, sd = 0.5)
  errors_at <- function(shift) {
    shifted <- d
    shifted$x <- d$t + shift
    fit <- mixture_regression(z ~ x, shifted, 2, n_starts = 4, seed = 1, tol = 1e-12)
    table <- as.data.frame(fit)
    table$std_error[table$term == "x"]
  }
  expect_equal(errors_at(1e5), errors_at(0), tolerance = 1e-8)
})
