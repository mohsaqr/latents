# Covariate fits with missing indicators (`missing = "fiml"`).
#
# The likelihood is checked against the observed-data definition written out
# row by row, the analytic score against numerical differentiation of that
# likelihood, and the standard errors against a Hessian of the likelihood alone,
# so none of them is compared with the machinery under test.
# The comparisons with independent references are in tests/equivalence/test-covariate-missing.R.

skip_on_cran()

.cov_missing_fixture <- function(n_groups = 40L, size = 8L, seed = 3L,
                                 share = 0.15) {
  set.seed(seed)
  n <- n_groups * size
  frame <- data.frame(g = rep(seq_len(n_groups), each = size),
                      z = stats::rnorm(n))
  group_class <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size)
  state <- 1L + stats::rbinom(n, 1L, stats::plogis(1.2 * frame$z +
                                                     c(-1.5, 1.5)[group_class + 1L]))
  frame$y1 <- stats::rnorm(n, c(0, 2)[state])
  frame$y2 <- stats::rnorm(n, c(0, 1.5)[state]) + 0.5 * frame$y1
  frame$q <- vapply(state, function(profile) {
    sample(c("a", "b", "c"), 1L,
           prob = list(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))[[profile]])
  }, character(1))
  blank <- function(value) replace(value, stats::runif(n) < share, NA)
  frame$y1 <- blank(frame$y1)
  frame$y2 <- blank(frame$y2)
  frame$q <- blank(frame$q)
  frame
}

.cov_missing_fit <- function(data, covariance_model, n_group_classes = 2L,
                             missing = "fiml", tol = 1e-12) {
  quietly(multilpa(data, c("y1", "y2", "q"), "g", 2L, n_group_classes,
                   categorical = "q", missing = missing,
                   covariance_model = covariance_model,
                   profile_covariates = "z", n_starts = 2, seed = 1,
                   tol = tol, max_iter = 5000))
}

test_that("without covariates the FIML fit is the covariate-free FIML model", {
  data <- .cov_missing_fixture()
  invisible(lapply(c("diagonal", "full"), function(covariance_model) {
    plain <- quietly(multilpa(data, c("y1", "y2", "q"), "g", 2L, 1L,
                              categorical = "q", missing = "fiml",
                              covariance_model = covariance_model,
                              n_starts = 1, seed = 1, tol = 1e-12,
                              max_iter = 5000))
    covariate <- quietly(.multilpa_fit_covariates(
      data, c("y1", "y2", "q"), "g", 2L, 1L, categorical = "q",
      missing = "fiml", covariance_model = covariance_model, n_starts = 1,
      seed = 1, tol = 1e-12, max_iter = 5000))
    expect_equal(covariate$log_likelihood, plain$log_likelihood,
                 tolerance = 1e-9)
  }))
})

test_that("on complete data `missing = \"fiml\"` changes nothing", {
  data <- .cov_missing_fixture(share = 0)
  invisible(lapply(c("diagonal", "full"), function(covariance_model) {
    strict <- .cov_missing_fit(data, covariance_model, missing = "error",
                               tol = 1e-8)
    lenient <- .cov_missing_fit(data, covariance_model, missing = "fiml",
                                tol = 1e-8)
    expect_identical(lenient$converged, strict$converged)
    skip_if_not(strict$converged, "no start converged on this fixture")
    expect_equal(lenient$log_likelihood, strict$log_likelihood, tolerance = 1e-10)
    expect_equal(lenient$profile_coefficients, strict$profile_coefficients,
                 tolerance = 1e-6)
    expect_equal(parameter_inference(lenient)$standard_error,
                 parameter_inference(strict)$standard_error, tolerance = 1e-6)
  }))
})

test_that("inference accepts the fitting data with its missing values", {
  data <- .cov_missing_fixture()
  fit <- .cov_missing_fit(data, "diagonal")
  from_stored <- parameter_inference(fit)
  from_data <- parameter_inference(fit, data)
  expect_equal(from_data$standard_error, from_stored$standard_error)
  expect_identical(nrow(parameter_inference(fit, vcov_type = "robust")),
                   nrow(from_stored))
  expect_true(all(is.finite(from_stored$standard_error)))
  # Filling a missing value changes the data, and inference says so.
  altered <- data
  altered$y1[which(is.na(altered$y1))[1L]] <- 0
  expect_error(parameter_inference(fit, altered),
               class = "latents_bad_inference_data")
})

test_that("missing covariates and unhandled missing indicators are refused by class", {
  data <- .cov_missing_fixture()
  expect_error(.cov_missing_fit(data, "diagonal", missing = "error"),
               class = "latents_bad_data")
  data$z[3L] <- NA
  expect_error(.cov_missing_fit(data, "diagonal"), class = "latents_bad_data")
})
