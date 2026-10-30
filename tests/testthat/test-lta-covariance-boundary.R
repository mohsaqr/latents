#' Generate Gaussian transition data with a controlled covariance geometry
#' @param geometry A flat covariance, narrow marginal, or interior covariance.
#' @return A three-occasion data frame with two well-separated profiles.
skip_on_cran()

lta_covariance_boundary_data <- function(geometry) {
  stopifnot(geometry %in% c("flat", "narrow", "interior"))
  set.seed(12)
  n <- 180L
  location <- sample(c(-5, 5), n, replace = TRUE)
  z <- rnorm(n)
  noise <- rnorm(n)
  data.frame(
    id = rep(seq_len(n / 3L), each = 3L), time = rep(seq_len(3L), n / 3L),
    a = location + if (geometry == "narrow") .03 * z else z,
    b = location + if (geometry == "flat") z + .03 * noise else noise)
}

test_that("general LTA detects a full covariance eigenvalue at its lower bound", {
  data <- lta_covariance_boundary_data("flat")
  invisible(lapply(c("invariant", "occasion"), function(measurement) {
    expect_warning(fit <- lta(
      data, c("a", "b"), "id", 2, time = "time", model = "VVV",
      transitions = "occasion", measurement = measurement,
      n_starts = 1, seed = 1, max_iter = 300, min_variance = .01),
      class = "latents_boundary")
    expect_true(fit$converged)
    expect_false(fit$probability_boundary)
    expect_true(fit$boundary)
    expect_length(fit$measurement, if (measurement == "occasion") 3L else 1L)
    invisible(lapply(fit$measurement, function(block) {
      # Marginal variances are interior; only the correlated direction is bound.
      expect_gt(min(block$variances), .1)
      eigenvalues <- apply(block$covariances, 3L, function(covariance) {
        min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values)
      })
      expect_equal(unname(eigenvalues), rep(.01, 2L), tolerance = 1e-10)
    }))
    expect_error(parameter_inference(fit), class = "latents_boundary_fit")
    expect_error(vcov(fit), class = "latents_boundary_fit")
  }))
})

test_that("general LTA retains diagonal variance boundary protection", {
  data <- lta_covariance_boundary_data("narrow")
  invisible(lapply(c("invariant", "occasion"), function(measurement) {
    expect_warning(fit <- lta(
      data, c("a", "b"), "id", 2, time = "time", model = "VVI",
      transitions = "occasion", measurement = measurement,
      n_starts = 1, seed = 1, max_iter = 300, min_variance = .01),
      class = "latents_boundary")
    expect_true(fit$converged)
    expect_false(fit$probability_boundary)
    expect_true(fit$boundary)
    invisible(lapply(fit$measurement, function(block) {
      expect_null(block$covariances)
      expect_equal(unname(block$variances[, 1L]), rep(.01, 2L), tolerance = 1e-10)
    }))
    expect_error(parameter_inference(fit), class = "latents_boundary_fit")
  }))
})

test_that("general LTA interior full covariances are not flagged as bound", {
  data <- lta_covariance_boundary_data("interior")
  invisible(lapply(c("invariant", "occasion"), function(measurement) {
    expect_no_warning(fit <- lta(
      data, c("a", "b"), "id", 2, time = "time", model = "VVV",
      transitions = "occasion", measurement = measurement,
      n_starts = 1, seed = 1, max_iter = 300, min_variance = .01))
    expect_true(fit$converged)
    expect_false(fit$probability_boundary)
    expect_false(fit$boundary)
    invisible(lapply(fit$measurement, function(block) {
      eigenvalues <- apply(block$covariances, 3L, function(covariance) {
        min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values)
      })
      expect_gt(min(eigenvalues), .1)
    }))
  }))
})
