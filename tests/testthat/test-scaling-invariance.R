# Unit invariance of a covariate fit, and of the standard errors built from it.
#
# Multiplying a covariate by a constant, or an indicator by a constant, is a
# reparameterization of the same model. The fitted log likelihood therefore
# cannot move (an indicator carries the density's Jacobian, which is a known
# constant), and an estimate carried back to the original units cannot move
# either. These are property tests of that invariance, not fixed-value
# regression tests: they would hold for any correct implementation.
#
# Tolerances. The optimization coordinates are the design columns divided by
# their own root-mean-square, so rescaling a column leaves the scaled design
# unchanged up to one rounding of a product and a quotient. The EM path is then
# identical and agreement is at machine precision, which is why the covariate
# assertions below are made at 1e-10 -- tighter than the testthat default of
# sqrt(.Machine$double.eps), not looser. The standard errors are the only place
# a looser tolerance is used: they come from a finite-difference Hessian, whose
# step is chosen on the scaled coordinates, so they agree to about seven digits
# rather than to machine precision. That is measured, not assumed.

# One covariate frame, and the same frame with its covariate in other units.
.scaling_covariate_frame <- function(seed = 22L, n = 600L) {
  set.seed(seed)
  z <- stats::rnorm(n)
  state <- stats::rbinom(n, 1L, stats::plogis(0.3 + 1.5 * z))
  data.frame(g = rep(seq_len(n / 20L), each = 20L), z = z,
             y = 5 * state + stats::rnorm(n, sd = 0.4))
}

# Warnings are collected rather than suppressed, so a test can say which
# qualification the fit raised instead of hiding all of them.
.scaling_fit_recording <- function(expression) {
  raised <- character()
  value <- withCallingHandlers(expression, warning = function(condition) {
    raised <<- c(raised, class(condition)[[1L]])
    invokeRestart("muffleWarning")
  })
  list(fit = value, warnings = raised)
}

test_that("the multinomial M-step is invariant to the units of a covariate", {
  set.seed(101)
  design <- cbind(1, stats::rnorm(200))
  truth <- matrix(c(0.2, 1.4), 2L, 1L)
  counts <- .multilpa_softmax(design, truth) * 3
  solve_at <- function(unit) {
    rescaled <- design
    rescaled[, 2L] <- design[, 2L] * unit
    step <- .multilpa_weighted_logits(rescaled, counts, matrix(0, 2L, 1L))
    c(intercept = unname(step$coefficients[1L, 1L]),
      slope_original_units = unname(step$coefficients[2L, 1L]) * unit,
      converged = as.numeric(step$converged),
      scaled_score = step$scaled_score)
  }
  units <- c(1, 1e-7, 1e3, 1e6)
  solved <- vapply(units, solve_at, numeric(4))

  expect_true(all(solved["converged", ] == 1))
  # The step really is solved, not merely returned: the unit-free score is at
  # most the tolerance the step is judged against.
  expect_true(all(solved["scaled_score", ] <= 1e-8 * (1 + sum(counts))))
  expect_equal(unname(solved["intercept", ]),
               rep(unname(solved["intercept", 1L]), length(units)),
               tolerance = 1e-10)
  expect_equal(unname(solved["slope_original_units", ]),
               rep(unname(solved["slope_original_units", 1L]), length(units)),
               tolerance = 1e-10)
  # It recovers the coefficients the counts were generated from.
  expect_equal(unname(solved[c("intercept", "slope_original_units"), 1L]),
               as.vector(truth), tolerance = 1e-5)
})

test_that("a covariate fit is invariant to the units its covariate arrives in", {
  frame <- .scaling_covariate_frame()
  units <- c(1, 1e-7, 1e3, 1e6)
  fitted <- lapply(units, function(unit) {
    rescaled <- frame
    rescaled$z <- frame$z * unit
    .scaling_fit_recording(
      fit_covariates(rescaled, "y", "g", n_profiles = 2, n_group_classes = 1,
                     profile_covariates = "z", n_starts = 1, seed = 3))
  })
  summarise <- function(recorded, unit) {
    fit <- recorded$fit
    c(log_likelihood = unname(fit$log_likelihood),
      slope_original_units = unname(fit$profile_coefficients[2L, 1L]) * unit,
      mean_1 = unname(fit$means[1L, 1L]),
      variance_1 = unname(fit$variances[1L, 1L]),
      converged = as.numeric(fit$converged))
  }
  table <- vapply(seq_along(units), function(i) summarise(fitted[[i]], units[[i]]),
                  numeric(5))

  expect_true(all(table["converged", ] == 1))
  expect_equal(unname(table["log_likelihood", ]),
               rep(unname(table["log_likelihood", 1L]), length(units)), tolerance = 1e-10)
  expect_equal(unname(table["slope_original_units", ]),
               rep(unname(table["slope_original_units", 1L]), length(units)),
               tolerance = 1e-10)
  expect_equal(unname(table["mean_1", ]),
               rep(unname(table["mean_1", 1L]), length(units)), tolerance = 1e-10)
  expect_equal(unname(table["variance_1", ]),
               rep(unname(table["variance_1", 1L]), length(units)), tolerance = 1e-10)
  # The effect is real, so a fit that quietly lost it would be caught here even
  # if every scale agreed on some other value.
  expect_lt(unname(table["slope_original_units", 1L]), -1)

  # Nothing is suppressed: at 1e-7 units the coefficient genuinely is above the
  # extreme-logit threshold, and the fit says so, at every other unit it is not.
  expect_equal(fitted[[1L]]$warnings, character())
  expect_equal(fitted[[3L]]$warnings, character())
  expect_equal(fitted[[4L]]$warnings, character())
  expect_true(length(fitted[[2L]]$warnings) >= 1L)
})

test_that("full-covariance inference survives an indicator held in large units", {
  set.seed(394)
  frame <- data.frame(g = rep(seq_len(60L), each = 5L), a = stats::rnorm(300),
                      b = stats::rnorm(300))
  units <- c(1, 1e3, 1e6)
  fitted <- lapply(units, function(unit) {
    rescaled <- frame
    rescaled$b <- frame$b * unit
    .scaling_fit_recording(
      fit_covariates(rescaled, c("a", "b"), "g", n_profiles = 1,
                     n_group_classes = 1, covariance_model = "full",
                     n_starts = 1, tol = 1e-10, seed = 1))
  })
  expect_true(all(vapply(fitted, function(recorded) recorded$fit$converged,
                         logical(1))))
  expect_equal(unlist(lapply(fitted, `[[`, "warnings")), character())

  # Rescaling one indicator moves the log likelihood by the density's Jacobian
  # and by nothing else.
  likelihoods <- vapply(fitted, function(recorded) recorded$fit$log_likelihood,
                        numeric(1))
  expect_equal(likelihoods - likelihoods[[1L]], -nrow(frame) * log(units),
               tolerance = 1e-10)

  # Inference is available at every unit, not only at the convenient one.
  inference <- lapply(fitted, function(recorded) parameter_inference(recorded$fit))
  expect_true(all(vapply(inference, function(table) {
    all(is.finite(table$standard_error))
  }, logical(1))))

  # Every estimate and every standard error carries the units of the parameter
  # it belongs to: a mean in `b` scales once, a covariance of `b` with itself
  # twice, and anything in `a` alone not at all.
  unit_power <- function(term) {
    parts <- unlist(strsplit(term, ":", fixed = TRUE))
    sum(parts == "b")
  }
  compare_at <- function(index) {
    joined <- merge(inference[[1L]], inference[[index]],
                    by = c("level", "outcome", "term", "parameter"),
                    suffixes = c("_base", "_scaled"))
    expect_gt(nrow(joined), 0L)
    factors <- unname(units[[index]]^vapply(joined$term, unit_power, numeric(1)))
    expect_equal(joined$estimate_scaled, joined$estimate_base * factors,
                 tolerance = 1e-8)
    expect_equal(joined$standard_error_scaled,
                 joined$standard_error_base * factors, tolerance = 1e-5)
  }
  compare_at(2L)
  compare_at(3L)
})

test_that("covariance coordinate scales name the units of each coordinate", {
  variances <- rbind(c(4, 9), c(1, 100))
  both <- .multilpa_covariance_coordinate_scale(variances, 2L, seq_len(2L))
  # Lower triangle of a 2x2, column-major: (1,1), (2,1), (2,2) per profile.
  expect_equal(both, c(1, 3, 1, 1, 10, 1))
  shared <- .multilpa_covariance_coordinate_scale(variances, 2L, 1L)
  expect_equal(shared, c(1, 3, 1))
  expect_equal(.multilpa_covariance_coordinate_scale(variances, 0L, 1L), numeric(0))
  expect_error(.multilpa_covariance_coordinate_scale(variances, 2L, 5L))
})

test_that("design scales are proportional to a column's units and never zero", {
  design <- cbind(1, c(-2, 0, 2, 4), 0)
  scale <- .multilpa_design_scale(design)
  expect_equal(scale[[1L]], 1)
  expect_equal(scale[[2L]], sqrt(mean(c(-2, 0, 2, 4)^2)))
  # A column with no variation at all still gets a usable scale rather than a
  # division by zero.
  expect_equal(scale[[3L]], 1)
  rescaled <- design
  rescaled[, 2L] <- design[, 2L] * 1e8
  expect_equal(.multilpa_design_scale(rescaled)[[2L]], scale[[2L]] * 1e8,
               tolerance = 1e-12)
  expect_error(.multilpa_design_scale(cbind(1, c(1, Inf))))
})

test_that("a non-finite membership design is refused by class, not by message", {
  design <- cbind(1, c(0.5, NA_real_, -0.5))
  counts <- matrix(c(0.6, 0.4, 0.5, 0.5, 0.2, 0.8), 3L, 2L, byrow = TRUE)
  expect_error(.multilpa_weighted_logits(design, counts, matrix(0, 2L, 1L)),
               class = "multilpa_bad_data")
  expect_error(
    .multilpa_weighted_logits(cbind(1, c(0, 1, 2)),
                              matrix(c(0.5, Inf, 0.5, 0.5, 0.5, 0.5), 3L, 2L),
                              matrix(0, 2L, 1L)),
    class = "multilpa_bad_data")
})
