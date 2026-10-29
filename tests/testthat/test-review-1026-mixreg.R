# ---- Review 2026-10-26: mixture-regression mathematics ----
#
# Calibration: the one-class limit is the ordinary GLM (glm(), MASS::glm.nb(),
# MASS::polr()), and the observed-information standard errors agree with a
# numerical Hessian of the log likelihood. Invariants: class relabelling,
# translating or rescaling a predictor. Regressions: the large-origin Gaussian
# least-squares shortcut with an empty or rank-deficient class, and an aliased
# constant column whose rounded mean is not exact.
# The comparisons with independent references are in tests/equivalence/test-review-1026-mixreg.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

.review_mixreg_data <- function(n = 240L, seed = 2026L) {
  .latents_local_seed(seed)
  frame <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  class_two <- stats::rbinom(n, 1L, stats::plogis(0.3 + 0.8 * frame$z)) == 1L
  frame$y <- ifelse(class_two, -1 - frame$x, 1 + 2 * frame$x) +
    stats::rnorm(n, sd = ifelse(class_two, 1, 0.5))
  frame$count <- stats::rpois(n, exp(ifelse(class_two, 1.5 - 0.3 * frame$x,
                                             0.2 + 0.8 * frame$x)))
  frame$exposure <- stats::runif(n, 1, 3)
  frame$successes <- stats::rbinom(n, 5L, stats::plogis(-0.2 + 0.7 * frame$x))
  frame$failures <- 5L - frame$successes
  frame$ordinal <- cut(ifelse(class_two, -1.5 * frame$x, 2 * frame$x) +
                         stats::rlogis(n), c(-Inf, -1, 0.5, 1.5, Inf),
                       ordered_result = TRUE)
  frame
}

test_that("relabelling the classes leaves the likelihood and the fitted means unchanged", {
  frame <- .review_mixreg_data()
  fit <- mixture_regression(y ~ x, frame, 2, membership = ~z, n_starts = 2, seed = 1,
                            vcov_type = "none")
  swapped <- .mixture_permute(fit$spec, fit$params, c(2L, 1L))
  original <- .mixture_expectation(fit$spec, fit$params)
  relabelled <- .mixture_expectation(fit$spec, swapped)
  expect_equal(relabelled$log_likelihood, original$log_likelihood, tolerance = 1e-12)
  expect_equal(relabelled$tau, original$tau[, 2:1], tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(.mixture_class_means(fit$spec, swapped),
               .mixture_class_means(fit$spec, fit$params)[, 2:1], ignore_attr = TRUE)
})

test_that("translating a predictor by 1e7 changes only the intercepts", {
  frame <- .review_mixreg_data()
  frame$shifted <- frame$x + 1e7
  compare <- function(formula, shifted_formula, family) {
    near <- mixture_regression(formula, frame, 2, family = family, n_starts = 2,
                               seed = 1, vcov_type = "none")
    far <- mixture_regression(shifted_formula, frame, 2, family = family,
                              n_starts = 2, seed = 1, vcov_type = "none")
    expect_equal(far$log_likelihood, near$log_likelihood, tolerance = 1e-9)
    expect_equal(unname(far$params$beta[2L, ]), unname(near$params$beta[2L, ]),
                 tolerance = 1e-8)
    expect_equal(.mixture_linear_predictors(far$spec, far$params),
                 .mixture_linear_predictors(near$spec, near$params), tolerance = 1e-8)
  }
  compare(y ~ x, y ~ shifted, "gaussian")
  compare(count ~ x, count ~ shifted, "poisson")
})

test_that("rescaling a predictor rescales its slopes and nothing else", {
  frame <- .review_mixreg_data()
  frame$scaled <- frame$x * 1e6
  compare <- function(response, family) {
    unit <- mixture_regression(stats::reformulate("x", response), frame, 2,
                               family = family, n_starts = 2, seed = 1,
                               vcov_type = "none")
    large <- mixture_regression(stats::reformulate("scaled", response), frame, 2,
                                family = family, n_starts = 2, seed = 1,
                                vcov_type = "none")
    expect_equal(large$log_likelihood, unit$log_likelihood, tolerance = 1e-9)
    slope <- function(fit, term) unname(fit$params$beta[term, ])
    expect_equal(slope(large, "scaled") * 1e6, slope(unit, "x"), tolerance = 1e-6)
  }
  compare("y", "gaussian")
  compare("count", "poisson")
  # The ordinal M-step solved a raw-ridge system whose ridge, scaled by the
  # slope's curvature, swamped the thresholds' and stalled short of the MLE.
  compare("ordinal", "ordinal")
})

test_that("an overfitted large-origin gaussian mixture fails only with a classed condition", {
  .latents_local_seed(17L)
  frame <- data.frame(x = round(stats::rnorm(20), 1) + 1e7)
  frame$y <- round(1 + 0.5 * (frame$x - 1e7) + stats::rnorm(20))
  outcome <- tryCatch(
    suppressWarnings(mixture_regression(y ~ x, frame, 4, n_starts = 10, seed = 17,
                                        vcov_type = "none")),
    latents_no_valid_start = function(condition) "no valid start")
  expect_true(identical(outcome, "no valid start") ||
                inherits(outcome, "latents_mixture_regression"))
})

test_that("an aliased constant column with an inexact mean adds nothing", {
  frame <- .review_mixreg_data(n = 99L)
  frame$constant <- 0.1
  # Where colMeans() accumulates in double precision (arm64 macOS) the mean of
  # 0.1 rounds, so a spread-is-zero test would miss the column; where it uses
  # extended precision (x86 Linux, Windows) the mean is exact. The fit must
  # ignore the column either way.
  for_family <- function(response, family) {
    plain <- mixture_regression(stats::reformulate("x", response), frame, 1,
                                family = family, n_starts = 0, vcov_type = "none")
    aliased <- mixture_regression(stats::reformulate(c("x", "constant"), response),
                                  frame, 1, family = family, n_starts = 0,
                                  vcov_type = "none")
    expect_equal(aliased$log_likelihood, plain$log_likelihood, tolerance = 1e-9)
    expect_equal(unname(aliased$params$beta["x", ]), unname(plain$params$beta["x", ]),
                 tolerance = 1e-8)
  }
  for_family("y", "gaussian")
  for_family("count", "poisson")
})
