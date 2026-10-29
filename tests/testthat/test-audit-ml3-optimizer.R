# Fixed-posterior regression for a regular equal-variance local maximum.
# The small fixture retains simulated data/posteriors and an independently
# optimized reference point, never read by the production optimizer.
# The comparisons with independent references are in tests/equivalence/test-audit-ml3-optimizer.R.

test_that("equal-variance step three preserves the observed-class Gaussian solution", {
  pieces <- list(posteriors = rbind(cbind(rep(1, 30), 0), cbind(rep(0, 30), 1)),
                 modal = rep(1:2, each = 30), n_classes = 2L)
  values <- c(rep(c(-2, -1, 0), 10), rep(c(3, 4, 5), 10))
  fit <- .multilpa_step_ml_fit(pieces, values, "equal")
  expect_equal(fit$means, c(-1, 4), tolerance = 1e-12)
  expect_equal(fit$variances, rep(2 / 3, 2), tolerance = 1e-12)
  expected <- sum(stats::dnorm(values, rep(c(-1, 4), each = 30), sqrt(2 / 3),
                                log = TRUE)) - 60 * log(2)
  expect_equal(fit$log_likelihood, expected, tolerance = 1e-12)
  expect_true(fit$converged)
  expect_false(fit$at_floor)
})

test_that("equal-variance optimizer basins respect outcome units and origin", {
  skip_on_cran()
  case <- readRDS(test_path("fixtures", "ml3-equal-rare-class.rds"))
  baseline <- .multilpa_step_ml_fit(case$pieces, case$values, "equal")
  invisible(lapply(list(c(scale = 1e-4, offset = 0), c(scale = 1e4, offset = 0),
                       c(scale = 1, offset = 1e8)), function(change) {
    values <- change[["offset"]] + change[["scale"]] * case$values
    fit <- .multilpa_step_ml_fit(case$pieces, values, "equal")
    expect_lt(max(abs((fit$means - change[["offset"]]) / change[["scale"]] -
                       baseline$means) / case$known_standard_error), 1e-6)
    expect_equal(fit$variances / change[["scale"]]^2, baseline$variances,
                 tolerance = 1e-7)
    expect_lt(abs(fit$log_likelihood + length(values) * log(change[["scale"]]) -
                    baseline$log_likelihood), 1e-6)
    expect_true(fit$converged)
    expect_false(fit$at_floor)
  }))
})
