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
  # stats::glm(family = gaussian()): least squares by QR.
  expect_reference(as.numeric(logLik(fit)), as.numeric(stats::logLik(reference)), reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, predictor in units of 1e10",
                   tolerance = 1e-10, measure = "relative")
  expect_reference(as.vector(fit$params$beta), stats::coef(reference), reference = "stats",
                   covers = "mixture_regression()", quantity = "coefficients, predictor in units of 1e10",
                   tolerance = 1e-9, measure = "relative")
  # Scaled, not relative: fitted values 1 + 2x cross zero, where a relative
  # measure would judge rounding noise; scaled is relative for |value| >= 1.
  expect_reference(predict(fit)$fitted, unname(stats::fitted(reference)), reference = "stats",
                   covers = c("mixture_regression()", "predict.latents_mixture_regression()"),
                   quantity = "fitted values, predictor in units of 1e10",
                   tolerance = 1e-8, measure = "scaled")
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
  expect_reference(as.numeric(logLik(scaled)), as.numeric(stats::logLik(reference)), reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, logistic regression with predictor in units of 0.01",
                   tolerance = 1e-9, measure = "relative")
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
  origins <- c(1e8, 1e10, 1e12)
  moved_fits <- lapply(origins, function(origin) {
    moved_data <- transform(data, x = x + origin)
    moved <- mixture_regression(y ~ x, moved_data, 1, n_starts = 0, vcov_type = "none")
    # Evaluate the independent GLM in a centered chart too; raw QR can itself
    # mistake a large-origin full-rank design for a constant predictor.
    moved_data$centered <- moved_data$x - mean(moved_data$x)
    reference <- stats::glm(y ~ centered, moved_data, family = stats::gaussian())
    expect_true(moved$converged)
    list(moved = moved, reference = reference,
         actual = data.frame(origin = origin, log_likelihood = as.numeric(logLik(moved)),
                             slope = moved$params$beta[2L, 1L]),
         expected = data.frame(origin = origin,
                               log_likelihood = as.numeric(stats::logLik(reference)),
                               slope = stats::coef(reference)[[2L]]))
  })
  stack <- function(part) do.call(rbind, lapply(moved_fits, `[[`, part))
  # stats::glm(family = gaussian()) on the centered predictor.
  expect_reference(stack("actual")[c("origin", "log_likelihood")],
                   stack("expected")[c("origin", "log_likelihood")], reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, predictor at a large origin",
                   tolerance = 1e-8, measure = "relative", key = "origin")
  expect_reference(stack("actual")[c("origin", "slope")],
                   stack("expected")[c("origin", "slope")], reference = "stats",
                   covers = "mixture_regression()", quantity = "slope, predictor at a large origin",
                   tolerance = 1e-10, measure = "relative", key = "origin")
  # The original intercept is stored as a double on its original scale, so the
  # fitted values carry the rounding of the origin: max(1e-8, eps * origin)
  # absolute, the original bound, one call per origin.
  fitted_at <- function(i) {
    expect_reference(predict(moved_fits[[i]]$moved)$fitted,
                     unname(stats::fitted(moved_fits[[i]]$reference)), reference = "stats",
                     covers = c("mixture_regression()", "predict.latents_mixture_regression()"),
                     quantity = sprintf("fitted values, predictor at origin %g", origins[i]),
                     tolerance = max(1e-8, .Machine$double.eps * origins[i]), measure = "absolute")
  }
  fitted_at(1L)
  fitted_at(2L)
  fitted_at(3L)
  expect_true(ordinary$converged)
})
