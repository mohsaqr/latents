#' Build a known-variance Gaussian diagnostic fixture
#' @return Data, posteriors and two full-rank Gaussian parameter blocks.
gaussian_diagnostics_fixture <- function() {
  data <- rbind(c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  parameters <- list(
    means = matrix(0, 2L, 2L),
    variances = matrix(1, 2L, 2L,
      dimnames = list(c("flat", "regular"), c("a", "b"))),
    covariances = array(c(matrix(c(1, .999, .999, 1), 2L), diag(2L)),
                        c(2L, 2L, 2L)))
  list(data = data, posteriors = matrix(rep(c(.3, .7), each = 4L), 4L),
       parameters = parameters)
}

test_that("Gaussian diagnostics describe a flat direction with moderate support", {
  fixture <- gaussian_diagnostics_fixture()
  parameters_before <- fixture$parameters
  rows <- rep(seq_len(4L), 20L)
  result <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data[rows, ], fixture$posteriors[rows, ],
    context = "occasion 2")
  expect_equal(result$profile, c("flat", "regular"))
  expect_equal(result$context, rep("occasion 2", 2L))
  expect_equal(result$gaussian_mass, c(24, 56))
  expect_equal(result$positive_rows, rep(80L, 2L))
  expect_equal(result$weight_effective_n, rep(80, 2L))
  expect_equal(result$min_variance_ratio, c(1, 1))
  expect_equal(result$min_standardized_eigenvalue, c(.001, 1), tolerance = 1e-12)
  expect_false(any(result$boundary))
  expect_identical(fixture$parameters, parameters_before)
})

test_that("standardized Gaussian ratios tolerate extreme changes of units", {
  fixture <- gaussian_diagnostics_fixture()
  reference <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data, fixture$posteriors)
  scales <- list(c(10, .1), c(-2, 3), rep(1e100, 2L), rep(1e-100, 2L),
                 c(1e100, 1e-100))
  invisible(lapply(scales, function(scale) {
    parameters <- fixture$parameters
    parameters$means <- sweep(parameters$means, 2L, scale, "*")
    parameters$variances <- sweep(parameters$variances, 2L, scale^2, "*")
    parameters$covariances <- array(vapply(seq_len(2L), function(profile) {
      parameters$covariances[, , profile] * outer(scale, scale)
    }, matrix(0, 2L, 2L)), c(2L, 2L, 2L))
    result <- .latents_gaussian_diagnostics(
      parameters, sweep(fixture$data, 2L, scale, "*"), fixture$posteriors)
    # The absolute lower-bound flag deliberately depends on the chosen units.
    columns <- c("gaussian_mass", "positive_rows", "weight_effective_n",
                 "min_variance_ratio", "min_standardized_eigenvalue")
    expect_equal(result[columns], reference[columns], tolerance = 1e-12)
  }))
})

test_that("tiny positive posterior mass does not imply a singular covariance", {
  fixture <- gaussian_diagnostics_fixture()
  fixture$parameters$covariances[, , 1L] <- diag(2L)
  posteriors <- matrix(rep(c(1e-12, 1 - 1e-12), each = 4L), 4L)
  result <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data, posteriors)
  expect_equal(result$gaussian_mass[1L], 4e-12, tolerance = 1e-15)
  expect_equal(result$weight_effective_n, rep(4, 2L))
  expect_equal(result$positive_rows, rep(4L, 2L))
  expect_equal(result$min_standardized_eigenvalue, c(1, 1))
  expect_false(any(result$boundary))
})

test_that("normalized Gaussian mass and weight effective n remain distinct", {
  fixture <- gaussian_diagnostics_fixture()
  reference <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data, fixture$posteriors, weights = rep(1, 4L),
    context = "individual profiles")
  # Normalizing to all eight units doubles active weights. Excluded extreme
  # rows cannot supply Gaussian support, change ratios, or change Kish n.
  appended <- .latents_gaussian_diagnostics(
    fixture$parameters,
    rbind(fixture$data, matrix(rep(c(1e100, -1e100), each = 4L), 4L)),
    rbind(fixture$posteriors, fixture$posteriors),
    weights = c(rep(2, 4L), rep(0, 4L)), context = "individual profiles")
  expect_equal(appended$gaussian_mass, reference$gaussian_mass * 2)
  expect_equal(appended$positive_rows, reference$positive_rows)
  expect_equal(appended$weight_effective_n, reference$weight_effective_n)
  expect_equal(appended$min_variance_ratio, reference$min_variance_ratio)
  expect_equal(appended$min_standardized_eigenvalue,
               reference$min_standardized_eigenvalue)
  expect_equal(appended$context, reference$context)
})

test_that("Gaussian diagnostics mask missing rows but retain partial observations", {
  fixture <- gaussian_diagnostics_fixture()
  data <- fixture$data
  data[1L, ] <- NA_real_
  data[2L, 2L] <- NA_real_
  data[3L, 1L] <- NA_real_
  result <- .latents_gaussian_diagnostics(
    fixture$parameters, data, fixture$posteriors, weights = c(5, 1, 2, 3))
  expect_equal(result$gaussian_mass, c(1.8, 4.2))
  expect_equal(result$positive_rows, rep(3L, 2L))
  expect_equal(result$weight_effective_n, rep(36 / 14, 2L))
  expect_equal(result$min_variance_ratio, rep(1 / .96, 2L))
  no_gaussian <- .latents_gaussian_diagnostics(
    list(), matrix(numeric(), 4L, 0L), fixture$posteriors,
    context = "categorical only")
  expect_s3_class(no_gaussian, "data.frame")
  expect_equal(nrow(no_gaussian), 0L)
  expect_identical(names(no_gaussian), names(result))
})

test_that("Gaussian diagnostics report unavailable reference variation explicitly", {
  fixture <- gaussian_diagnostics_fixture()
  fixture$data[, 2L] <- NA_real_
  result <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data, fixture$posteriors)
  expect_true(all(is.na(result$min_variance_ratio)))
  expect_true(all(is.na(result$min_standardized_eigenvalue)))
  expect_equal(result$gaussian_mass, c(1.2, 2.8))
  expect_output(.latents_print_gaussian_diagnostics(result),
                "insufficient observed reference variation")
  expect_output(.latents_print_gaussian_diagnostics(NULL), NA)
})

test_that("Gaussian diagnostic labels identify the component and density context", {
  fixture <- gaussian_diagnostics_fixture()
  individual <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data, fixture$posteriors,
    context = "individual profiles")
  group <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data, fixture$posteriors,
    context = "group classes")
  expect_equal(individual$profile, group$profile)
  expect_equal(group$context, rep("group classes", 2L))
  expect_output(.latents_print_gaussian_diagnostics(group),
                "0.001 \\(group classes, flat\\)")
  rownames(fixture$parameters$variances) <- NULL
  unnamed <- .latents_gaussian_diagnostics(
    fixture$parameters, fixture$data, fixture$posteriors)
  expect_equal(unnamed$profile, c("profile_1", "profile_2"))
})
