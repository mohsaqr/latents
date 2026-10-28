# Inference must transform with the units of the measured indicators. Exact
# rescaling of one fitted object isolates derivatives/inversion from fitting.

audit_lta_panel <- function() {
  withr::local_seed(8)
  n <- 80L
  path <- unlist(lapply(seq_len(n), function(i) {
    Reduce(function(previous, step) {
      if (stats::runif(1) < 0.85) previous else 3L - previous
    }, seq_len(4L), sample(1:2, 1L), accumulate = TRUE)
  }))
  data.frame(id = rep(seq_len(n), each = 5L), time = rep(1:5, n),
             z = stats::rnorm(n * 5L),
             y1 = stats::rnorm(n * 5L, c(0, 2.5)[path]),
             y2 = stats::rnorm(n * 5L, c(0, 1.5)[path], 0.8))
}

audit_lta_rescale <- function(fit, units) {
  fit$x <- sweep(fit$x, 2L, units, "*")
  fit$centers <- fit$centers * units
  fit$measurement <- lapply(fit$measurement, function(block) {
    block$means <- sweep(block$means, 2L, units, "*")
    block$variances <- sweep(block$variances, 2L, units^2, "*")
    block
  })
  # This fit has no active bounds; the rescaled flag stays false. Its stored
  # likelihood is not consumed by inference and need not be translated.
  fit
}

test_that("general LTA covariance is equivariant to measurement units", {
  fit <- lta(audit_lta_panel(), c("y1", "y2"), "id", 2, time = "time",
             transition_covariates = "z", n_starts = 1, seed = 2, tol = 1e-10)
  expect_true(fit$converged)
  reference <- lapply(c("observed", "robust", "opg"), function(type) vcov(fit, type = type))
  names(reference) <- c("observed", "robust", "opg")
  invisible(lapply(list(c(1e-4, 1e-4), c(1e-5, 1e-5), c(1e6, 1e6),
                        c(1e-5, 1e6)), function(units) {
    changed <- audit_lta_rescale(fit, units)
    map <- rep(1, length(coef(fit)))
    means <- startsWith(names(coef(fit)), "mean.")
    map[means] <- rep(units, length.out = sum(means))
    invisible(lapply(names(reference), function(type) {
      actual <- vcov(changed, type = type) / tcrossprod(map)
      # Standardize each entry by its marginal standard errors: a fixed
      # absolute tolerance could conceal the entire covariance at small units.
      scale <- sqrt(diag(reference[[type]]))
      expect_lt(max(abs(actual - reference[[type]]) / tcrossprod(scale)), 1e-6)
      actual_table <- parameter_inference(changed, vcov_type = type)
      reference_table <- parameter_inference(fit, vcov_type = type)
      expect_equal(actual_table$estimate, reference_table$estimate, tolerance = 1e-10)
      expect_equal(actual_table$standard_error, reference_table$standard_error,
                   tolerance = 1e-6)
    }))
  }))
})

test_that("weighted public LTA refits retain inference in tiny outcome units", {
  data <- audit_lta_panel()
  data$weight <- 1 + data$id %% 3L
  fit <- function(data) lta(data, c("y1", "y2"), "id", 2, time = "time",
                            transition_covariates = "z", weights = "weight",
                            n_starts = 1, seed = 2, tol = 1e-12)
  original <- fit(data)
  data[c("y1", "y2")] <- data[c("y1", "y2")] * 1e-5
  changed <- fit(data)
  expect_true(changed$converged)
  jacobian <- sum(original$sampling_weights) * 5 * 2 * log(1e-5)
  expect_equal(changed$log_likelihood + jacobian, original$log_likelihood,
               tolerance = 1e-9)
  means <- startsWith(names(coef(original)), "mean.")
  map <- ifelse(means, 1e-5, 1)
  reference <- vcov(original)
  expected_scale <- sqrt(diag(reference))
  expect_lt(max(abs(vcov(changed) / tcrossprod(map) - reference) /
                  tcrossprod(expected_scale)), 1e-5)
  expect_equal(parameter_inference(changed)$standard_error,
               parameter_inference(original)$standard_error, tolerance = 1e-5)
})
