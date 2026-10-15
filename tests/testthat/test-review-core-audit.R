# Regressions from the 2026-10-13 function audit: preprocessing contracts and
# covariance-family invariants must also hold where a variance floor is active.

core_audit_centering_data <- function() {
  set.seed(40)
  id <- rep(seq_len(50), each = 8)
  z <- stats::rnorm(length(id))
  class <- stats::rbinom(length(id), 1, stats::plogis(0.5 * z))
  data.frame(id = id, z = z,
    y = 100 + rep(stats::rnorm(50, sd = 3), each = 8) + 2 * class +
      stats::rnorm(length(id), sd = 0.7))
}

core_audit_centering_fit <- function(data, centering = "person", ...) {
  multilpa(data, "y", "id", 2, n_group_classes = 1,
    profile_covariates = "z", centering = centering,
    n_starts = 1, seed = 2, tol = 1e-10, ...)
}

test_that("covariate centering records the model and preserves raw-data inference", {
  data <- core_audit_centering_data()
  fit <- core_audit_centering_fit(data)
  expect_identical(fit$centering, "person")
  expect_lt(max(abs(rowsum(fit$indicator_data, fit$group_index))), 1e-10)
  expect_equal(get_results(fit, "data")$y, data$y, tolerance = 1e-12)
  supplied <- parameter_inference(fit, data)
  stored <- parameter_inference(fit)
  expect_equal(supplied, stored, tolerance = 1e-10)
  wrong <- data
  wrong$y[1] <- wrong$y[1] + 1
  expect_error(parameter_inference(fit, wrong), class = "latents_bad_inference_data")
  manual <- data
  manual$y <- data$y - stats::ave(data$y, data$id)
  reference <- core_audit_centering_fit(manual, "none")
  expect_equal(fit$means, reference$means, tolerance = 1e-10)
  expect_equal(fit$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  grand <- core_audit_centering_fit(data, "grand")
  plain <- core_audit_centering_fit(data, "none")
  expect_equal(grand$means + mean(data$y), plain$means, tolerance = 1e-7)
  expect_equal(grand$log_likelihood, plain$log_likelihood, tolerance = 1e-10)
})

test_that("covariate centering composes with FIML, full covariance and weights", {
  data <- core_audit_centering_data()
  data$y2 <- 20 + data$y / 10 + stats::rnorm(nrow(data), sd = 0.5)
  data$w <- rep(rep(1:3, length.out = 50), each = 8)
  data$y2[seq(1, nrow(data), by = 7)] <- NA_real_
  fit <- multilpa(data, c("y", "y2"), "id", 2, n_group_classes = 1,
    profile_covariates = "z", centering = "person", covariance_model = "full",
    missing = "fiml", weights = "w", n_starts = 1, seed = 2, tol = 1e-10)
  expect_identical(fit$centering, "person")
  returned <- get_results(fit, "data")
  expect_equal(returned$y, data$y, tolerance = 1e-12)
  expect_equal(returned$y2, data$y2, tolerance = 1e-12)
  expect_equal(parameter_inference(fit, data), parameter_inference(fit), tolerance = 1e-9)
})

test_that("all covariance structures retain their defining constraints at a floor", {
  angle <- 0.4
  rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
  scatter <- list(diag(c(0.001, 10)), rotation %*% diag(c(1, 0.01)) %*% t(rotation))
  weights <- c(10, 10)
  start <- array(rep(diag(0.2, 2), 2), c(2, 2, 2))
  cost <- function(blocks) sum(vapply(1:2, function(k) {
    weights[k] * as.numeric(determinant(blocks[[k]], logarithm = TRUE)$modulus) +
      sum(solve(blocks[[k]]) * scatter[[k]])
  }, numeric(1)))
  invisible(lapply(latents:::.multilpa_structures(), function(code) {
    result <- latents:::.multilpa_structure_covariances(scatter, weights, code,
      0.1, start = start)
    blocks <- lapply(1:2, function(k) matrix(result[, , k], 2, 2))
    values <- t(vapply(blocks, function(block) {
      eigen(block, symmetric = TRUE, only.values = TRUE)$values
    }, numeric(2)))
    volumes <- sqrt(apply(values, 1, prod))
    shapes <- values / volumes
    expect_gte(min(values), 0.1 * (1 - 1e-10))
    expect_lte(cost(blocks), cost(list(start[, , 1], start[, , 2])) + 1e-8)
    if (substr(code, 1, 1) == "E") expect_equal(volumes[1], volumes[2], tolerance = 1e-9)
    if (substr(code, 2, 2) == "E") expect_equal(shapes[1, ], shapes[2, ], tolerance = 1e-9)
    if (substr(code, 3, 3) == "E") {
      expect_equal(blocks[[1]] %*% blocks[[2]], blocks[[2]] %*% blocks[[1]], tolerance = 1e-9)
    }
    if (code %in% c("EII", "VII")) expect_equal(values[, 1], values[, 2])
  }))
})

test_that("structured floor-active profile models converge without an EM decrease", {
  set.seed(7)
  data <- data.frame(x = c(stats::rnorm(100, -10, 0.005), stats::rnorm(100, 10, 0.3)),
                     y = c(stats::rnorm(100, 0, 1), stats::rnorm(100, 0, 0.3)))
  invisible(lapply(c("EVI", "VEI", "EVV", "VEV", "VEE", "EVE", "VVE"), function(code) {
    expect_warning(fit <- lpa(data, c("x", "y"), 2, model = code,
      min_variance = 0.05, n_starts = 1, seed = 1, tol = 1e-10), class = "latents_boundary")
    expect_true(fit$converged)
    expect_true(all(diff(fit$log_likelihood_history) >= -1e-8))
  }))
})

test_that("constrained prior updates retain volume and shape at their own floors", {
  set.seed(20)
  data <- cbind(stats::rnorm(100, sd = 0.001), stats::rnorm(100))
  posterior <- cbind(rep(0.3, 100), rep(0.7, 100))
  invisible(lapply(c("EVI", "VEI", "VEV"), function(code) {
    prior <- prior_control(scale = if (code == "VEV") diag(1e-5, 2) else 1e-5)
    hyper <- latents:::.multilpa_prior_parameters(prior, data, 2, code)
    result <- latents:::.multilpa_prior_maximize(data, posterior, code, hyper, 0.1)
    values <- if (is.null(result$covariances)) result$variances else
      t(vapply(1:2, function(k) eigen(result$covariances[, , k], symmetric = TRUE,
        only.values = TRUE)$values, numeric(2)))
    volume <- sqrt(apply(values, 1, prod))
    expect_gte(min(values), 0.1 * (1 - 1e-10))
    if (code == "EVI") expect_equal(volume[1], volume[2], tolerance = 1e-9)
    else expect_equal(values[1, ] / volume[1], values[2, ] / volume[2], tolerance = 1e-9)
  }))
})

test_that("positive-definite starts outside the covariance family are not retained", {
  blocks <- list(diag(c(0.1, 10)), matrix(c(1, 0.1, 0.1, 0.3), 2))
  weights <- c(10, 20)
  start <- array(unlist(blocks), c(2, 2, 2))
  scatter <- Map(`*`, blocks, weights)
  # These are unrestricted optima, so accepting them as a cheaper fallback
  # would necessarily discard the requested equal-volume/shared-shape model.
  invisible(lapply(c("EEV", "VEV", "EVV", "VEE", "EVE", "VVE"), function(code) {
    result <- latents:::.multilpa_bounded_covariances(scatter, weights, code,
      0.1, start, 1000L, 1e-10)
    fitted <- lapply(1:2, function(k) matrix(result[, , k], 2))
    values <- t(vapply(fitted, function(block) eigen(block, symmetric = TRUE,
      only.values = TRUE)$values, numeric(2)))
    volume <- sqrt(apply(values, 1, prod))
    if (substr(code, 1, 1) == "E") expect_equal(volume[1], volume[2], tolerance = 1e-8)
    if (substr(code, 2, 2) == "E") expect_equal(values[1, ] / volume[1],
      values[2, ] / volume[2], tolerance = 1e-8)
    if (substr(code, 3, 3) == "E") expect_equal(fitted[[1]] %*% fitted[[2]],
      fitted[[2]] %*% fitted[[1]], tolerance = 1e-8)
  }))
})
