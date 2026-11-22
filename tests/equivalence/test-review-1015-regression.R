# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-review-1015-regression.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("a factor predictor can be shared by every class through `common`", {
  shared <- mixture_regression(score ~ hours + strategy, study_hours, n_classes = 1,
                               common = "strategy", n_starts = 0, seed = 1,
                               vcov_type = "none")
  reference <- stats::lm(score ~ hours + strategy, study_hours)
  n <- nrow(study_hours)
  # The Gaussian maximum likelihood at lm()'s least-squares fit, written out.
  maximum <- -n / 2 * (log(2 * pi * mean(stats::residuals(reference)^2)) + 1)
  expect_reference(as.numeric(logLik(shared)), maximum, reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one class with a shared factor", tolerance = 1e-8,
                   measure = "relative",
                   note = "Reference: -n/2 (log(2 pi RSS/n) + 1) from the residuals of stats::lm(score ~ hours + strategy).")
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
