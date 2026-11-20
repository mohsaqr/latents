# Comparisons with independent references moved from tests/testthat/test-review-followup.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("near-Poisson negative-binomial dispersion reaches the profile maximum", {
  set.seed(2)
  x1 <- stats::rnorm(2000L)
  y <- stats::rnbinom(2000L, size = 1e9, mu = 500)
  data <- data.frame(x1 = x1, c1 = y)
  fit <- lpa(data, c("x1", "c1"), 1, count = "c1",
             count_model = "negative_binomial", n_starts = 1)
  # One profile: the NB2 mean MLE is the sample mean, so the dispersion MLE is
  # the maximum of the profile likelihood in alpha alone.
  profile <- function(log_alpha) {
    sum(stats::dnbinom(y, size = exp(-log_alpha), mu = mean(y), log = TRUE))
  }
  reference <- stats::optimize(profile, c(log(1e-9), 0), maximum = TRUE, tol = 1e-12)
  expect_true(isTRUE(fit$converged))
  # Before the fix the Newton step stalled at 4.2e-7 against 1.4e-5.
  fitted <- log(as.vector(fit$count_dispersion))
  # optimize() is an independent maximiser of the profile likelihood written
  # above; 1e-4 relative on log alpha (about -11), the original tolerance.
  expect_reference(fitted, reference$maximum,
                   reference = "independent: NB2 profile likelihood in the dispersion, maximised by stats::optimize()",
                   covers = "lpa()", quantity = "log dispersion, near-Poisson counts, one profile",
                   tolerance = 1e-4, measure = "relative")
  expect_gte(profile(fitted), reference$objective - 1e-8)
  expect_true(all(is.finite(diag(vcov(fit)))))
})

test_that("the negative-binomial M-step ascends from either side of an indefinite region", {
  set.seed(2)
  stats::rnorm(2000L)
  y <- stats::rnbinom(2000L, size = 1e9, mu = 500)
  posteriors <- matrix(1, length(y), 1L)
  reference <- stats::optimize(function(log_alpha) {
    sum(stats::dnbinom(y, size = exp(-log_alpha), mu = mean(y), log = TRUE))
  }, c(log(1e-9), 0), maximum = TRUE, tol = 1e-12)$maximum
  # Starts below the maximum (where the objective is convex in log alpha) and
  # above it (the moment start, whose first full step overshot) both arrive.
  starts <- c(1e-8, 1e-7, 0.05, 1)
  reached <- vapply(starts, function(start) {
    point <- latents:::.latents_negative_binomial_maximize(y, posteriors, mean(y) * 1.01, start)
    log(point$dispersion)
  }, numeric(1))
  # The same maximum from every start; 1e-4 relative, the original tolerance.
  expect_reference(reached, rep(reference, length(starts)),
                   reference = "independent: NB2 profile likelihood in the dispersion, maximised by stats::optimize()",
                   covers = "lpa()", quantity = "log dispersion of the M-step from each start",
                   tolerance = 1e-4, measure = "relative", key = sprintf("start %g", starts),
                   note = "Package side: the internal negative-binomial M-step (.latents_negative_binomial_maximize) behind lpa(count_model = \"negative_binomial\").")
})
