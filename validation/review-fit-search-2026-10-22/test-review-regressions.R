# Standalone review regressions, deliberately outside tests/testthat.
# Run from the repository root with testthat::test_file().
# These assert corrected behavior and are expected to fail on this snapshot.
pkgload::load_all(quiet = TRUE)

capture_review_fit <- function(expr) {
  warnings <- character()
  fit <- withCallingHandlers(expr, warning = function(condition) {
    warnings <<- c(warnings, class(condition)[1L])
    invokeRestart("muffleWarning")
  }, message = function(condition) invokeRestart("muffleMessage"))
  list(fit = fit, warnings = warnings)
}

testthat::test_that("full-covariance LTA identifies an active eigenvalue bound", {
  set.seed(12)
  n <- 180L
  location <- sample(c(-5, 5), n, replace = TRUE)
  residual <- rnorm(n)
  data <- data.frame(id = rep(1:60, each = 3), time = rep(1:3, 60),
                     a = location + residual,
                     b = location + residual + rnorm(n, sd = .03))
  fitted <- capture_review_fit(lta(data, c("a", "b"), "id", 2, time = "time",
    model = "VVV", transitions = "occasion", n_starts = 1, seed = 1,
    max_iter = 300, min_variance = .01))
  minimum <- apply(fitted$fit$measurement[[1L]]$covariances, 3L, function(value) {
    min(eigen(value, symmetric = TRUE, only.values = TRUE)$values)
  })
  testthat::expect_true(fitted$fit$converged)
  testthat::expect_equal(unname(minimum), rep(.01, 2L), tolerance = 1e-9)
  testthat::expect_true(fitted$fit$boundary)
  testthat::expect_error(parameter_inference(fitted$fit),
                         class = "latents_boundary_fit")
})

testthat::test_that("an excluded outlier does not derail all default starts", {
  set.seed(803)
  data <- data.frame(y = c(rnorm(50, -3, .5), rnorm(50, 3, .5)), w = 1)
  fit <- function(value) capture_review_fit(lpa(
    rbind(data, data.frame(y = value, w = 0)), "y", 2,
    weights = "w", n_starts = 10, seed = 1))$fit
  ordinary <- fit(0)
  outlier <- fit(10000)
  testthat::expect_true(ordinary$converged && outlier$converged)
  testthat::expect_equal(outlier$log_likelihood, ordinary$log_likelihood,
                         tolerance = 1e-7)
})

testthat::test_that("non-Gaussian rows cannot hide a thin Gaussian profile", {
  set.seed(73)
  data <- as.data.frame(matrix(rnorm(80 * 5), 80, 5))
  names(data) <- letters[1:5]
  data[1:3, ] <- data[1:3, ] + 10
  data$categ <- c(rep("rare", 3), rep("common", 77))
  empty <- data[rep(1L, 300L), ]
  empty[, letters[1:5]] <- NA_real_
  start <- list(means = rbind(rep(10, 5), rep(0, 5)),
    variances = matrix(1, 2, 5),
    profile_probabilities = matrix(c(3, 77) / 80, 1), group_probabilities = 1,
    response_probabilities = list(categ =
      rbind(c(1e-10, 1 - 1e-10), c(1 - 1e-10, 1e-10))))
  fitted <- capture_review_fit(lpa(rbind(data, empty), c(letters[1:5], "categ"),
    2, categorical = "categ", start = start, n_starts = 1, seed = 1,
    missing = "fiml"))
  testthat::expect_equal(unname(colSums(fitted$fit$subject_posteriors[1:80, ])),
                         c(3, 77), tolerance = 1e-6)
  testthat::expect_true("latents_thin_profiles" %in% fitted$warnings)
})

testthat::test_that("cross-level thinness uses weighted Gaussian support", {
  set.seed(7)
  data <- data.frame(g = rep(1:40, each = 2), a = rnorm(80), b = rnorm(80))
  data[1:6, c("a", "b")] <- data[1:6, c("a", "b")] + 20
  data$weight <- ifelse(data$g <= 3, .01, 1)
  fitted <- capture_review_fit(multilpa(data, c("a", "b"), "g", 2, 1,
    family = "full_cross_level", weights = "weight", n_starts = 1,
    seed = 1, max_iter = 300))
  support <- colSums(fitted$fit$subject_posteriors *
    fitted$fit$sampling_weights[fitted$fit$sufficient_statistics$index])
  testthat::expect_true(min(support) < 1)
  testthat::expect_true(fitted$fit$converged)
  testthat::expect_true("latents_thin_profiles" %in% fitted$warnings)
})

testthat::test_that("cross-level thinness also checks the group Gaussian mixture", {
  set.seed(19)
  data <- data.frame(g = rep(1:40, each = 3), a = rnorm(120), b = rnorm(120))
  data[data$g <= 3, c("a", "b")] <- data[data$g <= 3, c("a", "b")] + 20
  fitted <- capture_review_fit(multilpa(data, c("a", "b"), "g",
    n_profiles = 1, n_group_classes = 2, family = "full_cross_level",
    between_variance = "varying", n_starts = 1, seed = 1, max_iter = 300))
  testthat::expect_equal(unname(colSums(fitted$fit$group_posteriors)), c(3, 37),
                         tolerance = 1e-6)
  testthat::expect_true("latents_thin_profiles" %in% fitted$warnings)
})
