# Comparisons with independent references moved from tests/testthat/test-review-1015-regression.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("a factor predictor can be shared by every class through `common`", {
  shared <- mixture_regression(score ~ hours + strategy, study_hours, n_classes = 1,
                               common = "strategy", n_starts = 0, seed = 1,
                               vcov_type = "none")
  reference <- stats::lm(score ~ hours + strategy, study_hours)
  n <- nrow(study_hours)
  maximum <- -n / 2 * (log(2 * pi * mean(stats::residuals(reference)^2)) + 1)
  expect_equal(shared$log_likelihood, maximum, tolerance = 1e-8)
  expect_true("coefficient.common.strategysurface" %in% names(coef(shared)))
  # A term is matched by its variables, in whatever order they are written.
  crossed <- mixture_regression(score ~ hours * sleep, study_hours, n_classes = 2,
                                common = ~ sleep:hours, n_starts = 1, seed = 1,
                                vcov_type = "none")
  expect_true("coefficient.common.hours:sleep" %in% names(coef(crossed)))
  expect_error(mixture_regression(score ~ hours, study_hours, n_classes = 2,
                                  common = "strategy", n_starts = 1),
               class = "latents_bad_argument")
})
