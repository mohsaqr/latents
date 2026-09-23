test_that("mixed full covariance inference matches Gaussian and binomial information", {
  set.seed(7601)
  dat <- data.frame(g = rep(seq_len(30), each = 8),
                    a = rnorm(240), b = rnorm(240), u = sample(1:2, 240, TRUE))
  dat$b <- dat$b + 0.4 * dat$a
  fit <- multilpa(dat, c("a", "b", "u"), "g", 1, 1, categorical = "u", covariance_model = "full", n_starts = 1)
  information <- parameter_inference(fit, dat)
  covariance <- cov(dat[c("a", "b")]) * (nrow(dat) - 1) / nrow(dat)
  probability <- mean(dat$u == 1)
  covariance_errors <- sqrt((covariance^2 + outer(diag(covariance),
                                                diag(covariance))) / nrow(dat))
  expected <- c(sqrt(diag(covariance) / nrow(dat)),
                covariance_errors[lower.tri(covariance, diag = TRUE)],
                rep(sqrt(probability * (1 - probability) / nrow(dat)), 2), 0, 0)
  expect_equal(information$standard_error, unname(expected), tolerance = 1e-7)
  expect_length(coef(fit, scale = "unconstrained"), fit$n_parameters)
  expect_false(any(grepl("NA", names(coef(fit)))))
  expect_true(all(is.na(information$p_value[c(3, 5, 6, 7, 8, 9)])))
  expect_true(is.finite(information$p_value[4]))

  # Closed-form cluster influence contributions for empirical Gaussian moments
  # and the binomial proportion give the entire robust covariance independently.
  residual <- sweep(as.matrix(dat[c("a", "b")]), 2L, colMeans(dat[c("a", "b")]))
  influence <- cbind(residual,
    residual[, 1]^2 - covariance[1, 1],
    residual[, 1] * residual[, 2] - covariance[1, 2],
    residual[, 2]^2 - covariance[2, 2],
    (dat$u == 1) - probability, (dat$u == 2) - (1 - probability), 0, 0)
  expected_robust <- crossprod(rowsum(influence, dat$g)) / nrow(dat)^2
  robust <- parameter_inference(fit, dat, vcov_type = "robust")
  expect_equal(unname(attr(robust, "covariance")), unname(expected_robust),
               tolerance = 1e-7)
})

test_that("all-categorical full covariance inference has no Gaussian parameters", {
  set.seed(7602)
  dat <- data.frame(g = rep(seq_len(30), each = 8), u = sample(1:3, 240, TRUE))
  fits <- lapply(c("diagonal", "full"), function(covariance_model) {
    multilpa(dat, "u", "g", 1, 1, categorical = "u", covariance_model = covariance_model, n_starts = 1)
  })
  expect_equal(coef(fits[[1]]), coef(fits[[2]]))
  invisible(lapply(c("observed", "robust"), function(vcov_type) {
    expect_equal(parameter_inference(fits[[1]], dat, vcov_type = vcov_type),
                 parameter_inference(fits[[2]], dat, vcov_type = vcov_type))
  }))
  renamed <- dat
  renamed$u <- letters[dat$u]
  expect_error(parameter_inference(fits[[1]], renamed),
               class = "latents_bad_inference_data")
})

test_that("full covariance inference is invariant to indicator units", {
  set.seed(7605)
  dat <- data.frame(g = rep(seq_len(30), each = 8), a = rnorm(240), b = rnorm(240))
  dat$b <- dat$b + 0.4 * dat$a
  fit <- multilpa(dat, c("a", "b"), "g", 1, 1,
                 covariance_model = "full", n_starts = 1)
  scaled_data <- transform(dat, a = a * 1e6, b = b * 1e4)
  scaled <- multilpa(scaled_data, c("a", "b"), "g", 1, 1,
                    covariance_model = "full", n_starts = 1)
  factors <- c(1e6, 1e4, 1e12, 1e10, 1e8, 1, 1)
  invisible(lapply(c("observed", "robust"), function(vcov_type) {
    original <- parameter_inference(fit, dat, vcov_type = vcov_type)
    rescaled <- parameter_inference(scaled, scaled_data, vcov_type = vcov_type)
    expect_equal(rescaled$standard_error / factors, original$standard_error,
                 tolerance = 1e-7)
  }))
})

test_that("Wald inference refuses active categorical probability constraints", {
  dat <- data.frame(g = rep(seq_len(30), each = 8),
                    u = c(rep(1, 139), 2, rep(3, 100)))
  fit <- multilpa(dat, "u", "g", 1, 1, categorical = "u",
                 min_probability = 0.01, n_starts = 1)
  expect_equal(min(fit$response_probabilities$u), 0.01, tolerance = 1e-10)
  expect_error(parameter_inference(fit, dat), "bound-active")
})

test_that("covariate inference supports empty membership coefficient blocks", {
  set.seed(7603)
  dat <- data.frame(g = rep(seq_len(30), each = 8), y = rnorm(240, 3, 2))
  # The covariate estimator with no covariates, which `multilpa()` cannot be
  # asked for: naming no covariate requests the covariate-free model.
  fit <- latents:::.multilpa_fit_covariates(dat, "y", "g", 1, 1, n_starts = 1)
  information <- parameter_inference(fit, dat)
  variance <- mean((dat$y - mean(dat$y))^2)
  expect_equal(information$estimate, c(mean(dat$y), variance))
  expect_equal(information$standard_error,
               c(sqrt(variance / nrow(dat)), variance * sqrt(2 / nrow(dat))),
               tolerance = 1e-7)
  expect_equal(unname(coef(fit)), information$estimate)
  expect_equal(unname(confint(fit, data = dat)),
               unname(cbind(information$conf_low, information$conf_high)))
  expect_equal(fit$n_informative, nrow(dat))
  residual <- dat$y - mean(dat$y)
  influence <- cbind(residual, residual^2 - variance)
  robust <- parameter_inference(fit, dat, vcov_type = "robust")
  expect_equal(robust$standard_error,
               sqrt(diag(crossprod(rowsum(influence, dat$g)))) / nrow(dat),
               ignore_attr = TRUE, tolerance = 1e-7)

  dat$y <- dat$y * 1e6
  scaled_fit <- latents:::.multilpa_fit_covariates(dat, "y", "g", 1, 1, n_starts = 1)
  scaled <- parameter_inference(scaled_fit, dat)
  expect_equal(scaled$standard_error / c(1e6, 1e12), information$standard_error,
               tolerance = 1e-7)
  scaled_robust <- parameter_inference(scaled_fit, dat, vcov_type = "robust")
  expect_equal(scaled_robust$standard_error / c(1e6, 1e12), robust$standard_error,
               tolerance = 1e-7)
})

test_that("covariate inference refuses nonregular fits and altered fitting data", {
  set.seed(7604)
  dat <- data.frame(g = rep(seq_len(30), each = 8), z = rnorm(240))
  profile <- rbinom(240, 1, plogis(dat$z))
  dat$y <- 6 * profile + rnorm(240, sd = 0.5)
  fit <- multilpa(dat, "y", "g", 2, 1, "z", n_starts = 1, tol = 1e-12)
  expect_equal(nrow(parameter_inference(fit, dat)), fit$n_parameters)
  unconverged <- fit
  unconverged$converged <- FALSE
  expect_error(parameter_inference(unconverged, dat), "converged")
  boundary <- fit
  boundary$boundary <- TRUE
  expect_error(parameter_inference(boundary, dat), "bound-active")
  few_groups <- fit
  few_groups$n_groups <- fit$n_parameters
  expect_error(parameter_inference(few_groups, dat, vcov_type = "robust"),
               class = "latents_too_few_groups")
  wrong_group <- dat
  wrong_group$g[1] <- 2L
  wrong_covariate <- dat
  wrong_covariate$z[1] <- wrong_covariate$z[1] + 1
  wrong_indicator <- dat
  wrong_indicator$y[1:2] <- rev(wrong_indicator$y[1:2])
  invisible(lapply(list(wrong_group, wrong_covariate, wrong_indicator, dat[-1, ]),
    function(wrong) {
      expect_error(parameter_inference(fit, wrong), class = "latents_bad_inference_data")
    }))

  # Changing predictor units rescales only that coefficient and its error.
  scaled <- fit
  scaled_data <- dat
  scaled_data$z <- dat$z * 1e6
  scaled$profile_design <- lapply(fit$profile_design, function(design) {
    design[, ncol(design)] <- design[, ncol(design)] * 1e6
    design
  })
  scaled$profile_coefficients["z", ] <- fit$profile_coefficients["z", ] / 1e6
  information <- parameter_inference(fit, dat)
  scaled_information <- parameter_inference(scaled, scaled_data)
  factor <- ifelse(information$term == "z", 1e6, 1)
  expect_equal(scaled_information$standard_error * factor,
               information$standard_error, tolerance = 1e-6)
})

test_that("covariate posteriors preserve priors under huge common density offsets", {
  x <- matrix(c(-1e8, 1e8), 2, 1)
  parameters <- list(means = matrix(0, 2, 1), variances = matrix(1, 2, 1))
  design <- list(matrix(rep(c(1, 0), 2), 2, 2, byrow = TRUE),
                 matrix(rep(c(0, 1), 2), 2, 2, byrow = TRUE))
  result <- .multilpa_cov_expectation(x, c(1L, 2L), parameters, design,
    matrix(1, 2, 1), matrix(qlogis(c(0.8, 0.3)), 2, 1), matrix(qlogis(0.6), 1, 1))
  expect_equal(result$group_posteriors, matrix(c(0.6, 0.4), 2, 2, byrow = TRUE),
               tolerance = 1e-14)
  expect_equal(result$subject_posteriors, matrix(c(0.6, 0.4), 2, 2, byrow = TRUE),
               tolerance = 1e-14)
  expect_equal(result$log_likelihood, sum(dnorm(x, log = TRUE)))
})

test_that("covariate likelihood retains finite log priors below probability underflow", {
  parameters <- list(means = matrix(c(0, 50), 2, 1), variances = matrix(1, 2, 1))
  result <- .multilpa_cov_expectation(matrix(0, 1, 1), 1L, parameters,
    list(matrix(1, 1, 1)), matrix(1, 1, 1), matrix(-1000, 1, 1),
    matrix(numeric(0), 1, 0))
  expect_equal(result$log_likelihood, -1000 + dnorm(0, log = TRUE), tolerance = 1e-12)
  expect_equal(result$subject_posteriors[1, 1], 1, tolerance = 1e-12)
})

test_that("large common logits and group evidence preserve multinomial normalization", {
  design <- matrix(1, 2, 1)
  coefficients <- matrix(c(1e16, 1e16), 1, 2)
  expect_equal(.multilpa_softmax(design, coefficients),
               matrix(c(0.5, 0.5, 0), 2, 3, byrow = TRUE))
  expect_equal(.multilpa_log_softmax(design, coefficients)[, 1:2],
               matrix(-log(2), 2, 2))
  parameters <- list(means = matrix(c(0, 2e8), 2, 1), variances = matrix(1, 2, 1))
  profile_design <- list(matrix(c(1, 0), 1, 2), matrix(c(0, 1), 1, 2))
  result <- .multilpa_cov_expectation(matrix(0, 1, 1), 1L, parameters,
    profile_design, matrix(1, 1, 1), matrix(-1e16, 2, 1), matrix(qlogis(0.25), 1, 1))
  expect_equal(result$group_posteriors, matrix(c(0.25, 0.75), 1, 2), tolerance = 1e-14)
  expect_equal(result$subject_posteriors, matrix(c(1, 0), 1, 2), tolerance = 1e-14)
})

test_that("observed information rejects a negative definite Hessian", {
  expect_error(.multilpa_observed_hessian(
    function(displacement) -sum(displacement^2),
    function(displacement) -2 * displacement, c(1, 1), 1e-4),
    "not positive definite")
})
