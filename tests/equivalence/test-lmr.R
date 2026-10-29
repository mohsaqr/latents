# Comparisons with independent references moved from tests/testthat/test-lmr.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("the Imhof tail matches exact weighted chi-square probabilities", {
  # One weight, or equal weights, give a scaled chi-square.
  invisible(lapply(c(0.05, 0.5, 3, 9, 30), \(q) {
    expect_equal(.latents_imhof_upper(q, 1), stats::pchisq(q, 1, lower.tail = FALSE),
                 tolerance = 1e-8)
    expect_equal(.latents_imhof_upper(q, rep(1, 3)),
                 stats::pchisq(q, 3, lower.tail = FALSE), tolerance = 1e-8)
    expect_equal(.latents_imhof_upper(q, 2.5), stats::pchisq(q / 2.5, 1, lower.tail = FALSE),
                 tolerance = 1e-8)
  }))
  # Two weights, against the convolution of their distributions: P(a X + b Y >
  # q) with y = t^2, so the chi-square(1) density times dy is 2 dnorm(t) dt,
  # smooth at zero. One weight small, so the tail decays slowly; one
  # negative, as an eigenvalue of Vuong's matrix can be.
  two <- function(q, a, b) {
    f <- function(t) stats::pchisq(pmax((q - b * t^2) / a, 0), 1, lower.tail = FALSE) *
      2 * stats::dnorm(t)
    if (b > 0) {
      kink <- sqrt(q / b)
      stats::integrate(f, 0, kink, rel.tol = 1e-12)$value +
        2 * stats::pnorm(kink, lower.tail = FALSE)
    } else stats::integrate(f, 0, Inf, rel.tol = 1e-12)$value
  }
  grid <- expand.grid(q = c(1e-6, 0.5, 5.5, 30, 80), a = c(0.9, 3, 0.01),
                      b = c(0.05, -0.4, -3, 50))
  expect_equal(mapply(function(q, a, b) .latents_imhof_upper(q, c(a, b)),
                      grid$q, grid$a, grid$b),
               mapply(two, grid$q, grid$a, grid$b), tolerance = 1e-10)
  # A statistic of zero with weights of both signs has no oscillation to
  # extrapolate: P(a X_m - b Y_k > 0) = P(F(m, k) > b k / (a m)), exactly.
  # Extrapolating partial sums returned 5.6e-7 for 4.2e-7 here.
  cases <- expand.grid(m = c(1, 3), b = c(0.1, 10, 1e4), k = c(1, 3, 6))
  expect_equal(mapply(function(m, b, k) .latents_imhof_upper(0, c(rep(1, m), rep(-b, k))),
                      cases$m, cases$b, cases$k),
               with(cases, stats::pf(b * k / m, m, k, lower.tail = FALSE)),
               tolerance = 1e-10)
  # The accuracy is absolute: about 1e-12, here on a probability of 4e-7.
  expect_lt(abs(.latents_imhof_upper(0, c(1, -1e4, -1e4, -1e4)) -
                  stats::pf(3e4, 1, 3, lower.tail = FALSE)), 1e-12)
  # The probability does not depend on the weights' units.
  invisible(lapply(c(1e-12, 1e-6, 1e6), \(size) {
    expect_equal(.latents_imhof_upper(0, size * c(1, -2)),
                 stats::pf(2, 1, 1, lower.tail = FALSE), tolerance = 1e-10)
    expect_equal(.latents_imhof_upper(5 * size, size * c(2, 1, 0.5, -0.4)),
                 .latents_imhof_upper(5, c(2, 1, 0.5, -0.4)), tolerance = 1e-10)
  }))
})

test_that("VLMR and adjusted LMR agree with Latent GOLD", {
  skip_on_cran()
  # Latent GOLD 6.1, `latent Cluster nominal 2:3` with one variance per
  # indicator shared by the classes, on these simulated data (2 true classes):
  # LR 8.3027, adjusted 8.0439, reference mean 4.1216 and standard deviation
  # 2.8763, p-values 0.0848 (VLMR) and 0.0941 (adjusted).
  set.seed(2024)
  n <- 500
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  data <- data.frame(y1 = stats::rnorm(n, c(0, 1.6)[class]),
                     y2 = stats::rnorm(n, c(0, 1.4)[class]),
                     y3 = stats::rnorm(n, c(0, 1.2)[class]),
                     y4 = stats::rnorm(n, c(0, 1)[class]))
  models <- quietly(enumerate_lpa(data, paste0("y", 1:4), n_profiles = 2:3,
                                  model = "EEI", n_starts = 50, seed = 1, tol = 1e-12))
  test <- get_results(models, "lmr")
  within <- function(actual, printed) {
    expect_true(abs(actual - printed) <= 5e-5, info = format(actual, digits = 10))
  }
  within(test$lr, 8.3027)
  within(test$adjusted_lr, 8.0439)
  within(test$reference_mean, 4.1216)
  within(test$reference_sd, 2.8763)
  within(test$vlmr_p, 0.0848)
  within(test$almr_p, 0.0941)
  # The adjustment is LR / (1 + 1 / (d log n)), d the added parameters.
  expect_equal(test$adjusted_lr, test$lr / (1 + 1 / (5 * log(500))), tolerance = 1e-12)
})
