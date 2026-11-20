# Comparisons with independent references moved from tests/testthat/test-review-1026-mixreg-r2.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("standard errors do not depend on a predictor's origin (N6)", {
  set.seed(3)
  n <- 300
  d <- data.frame(t = stats::rnorm(n))
  d$y <- 1 + 2 * d$t + stats::rnorm(n)
  d$count <- stats::rpois(n, exp(0.5 + 0.3 * d$t))
  slope_error <- function(fit) {
    table <- as.data.frame(fit)
    table$std_error[table$term == "x"]
  }
  shifts <- c(1e4, 5e4, 1e6)
  shifted_at <- function(shift) transform(d, x = t + shift)
  origins <- sprintf("origin %g", shifts)
  # One gaussian class at each origin against glm(); the maximum-likelihood
  # variance divides by n, glm()'s dispersion by n - 2.
  expect_reference(
    stats::setNames(vapply(shifts, function(shift) {
      slope_error(mixture_regression(y ~ x, shifted_at(shift), 1, n_starts = 0))
    }, numeric(1)), origins),
    stats::setNames(vapply(shifts, function(shift) {
      reference <- stats::glm(y ~ x, data = shifted_at(shift))
      sqrt(stats::vcov(reference)[2, 2] * (n - 2) / n)
    }, numeric(1)), origins),
    reference = "stats",
    covers = c("mixture_regression()", "as.data.frame.latents_mixture_regression()"),
    quantity = "slope standard error, one gaussian class", tolerance = 1e-7,
    measure = "relative")
  # One Poisson class at each origin against glm(); 1e-6, the original tolerance.
  expect_reference(
    stats::setNames(vapply(shifts, function(shift) {
      slope_error(mixture_regression(count ~ x, shifted_at(shift), 1, family = "poisson",
                                     n_starts = 0))
    }, numeric(1)), origins),
    stats::setNames(vapply(shifts, function(shift) {
      reference <- stats::glm(count ~ x, data = shifted_at(shift), family = stats::poisson)
      sqrt(stats::vcov(reference)[2, 2])
    }, numeric(1)), origins),
    reference = "stats",
    covers = c("mixture_regression()", "as.data.frame.latents_mixture_regression()"),
    quantity = "slope standard error, one poisson class", tolerance = 1e-6,
    measure = "relative")
  # Two classes: the slope's errors are the same at origin 0 and 1e5 (a unit
  # check of invariance).
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
