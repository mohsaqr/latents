test_that("observed information reproduces analytic Gaussian standard errors", {
  set.seed(329)
  dat <- data.frame(group = rep(seq_len(20), each = 10),
                    y1 = rnorm(200, 2, 0.7), y2 = rnorm(200, -1, 3))
  fit <- multilpa(dat, c("y1", "y2"), "group", 1, 1, n_starts = 1)
  information <- parameter_inference(fit, dat)
  expected <- c(sqrt(fit$variances / 200), sqrt(2 / 200) * fit$variances, 0, 0)
  expect_equal(unname(information$standard_error), unname(expected), tolerance = 1e-7)
  expect_equal(attr(information, "covariance"), t(attr(information, "covariance")), tolerance = 1e-12)
  expect_lt(attr(information, "scaled_score"), 1e-7)
  expect_gt(attr(information, "condition_ratio"), 0.1)
  expect_equal(unname(coef(fit)), information$estimate)
  expect_length(coef(fit, scale = "unconstrained"), fit$n_parameters)
  expect_equal(vcov(fit, data = dat), attr(information, "covariance"))
  expect_equal(vcov(fit, data = dat, scale = "unconstrained"), attr(information, "covariance_unconstrained"))
  expect_equal(unname(confint(fit, data = dat)),
               unname(cbind(information$conf_low, information$conf_high)))
  expect_equal(confint(fit, parm = 1, level = 0.9, data = dat),
    confint(fit, parm = names(coef(fit))[1], level = 0.9, data = dat))
  expect_equal(as.numeric(confint(fit, parm = 1, level = 0.9, data = dat)),
    coef(fit)[1] + c(-1, 1) * qnorm(0.95) * expected[1], ignore_attr = TRUE)
  expect_error(confint(fit, parm = "unknown", data = dat), "parm")
  expect_error(confint(fit, parm = -1, data = dat), "indices")
  expect_error(confint(fit, level = 1.2, data = dat))
})

test_that("score and probability Jacobian agree with independent central differences", {
  reference <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-synthetic-varying.rds"))
  fit <- multilpa(reference$data, c("y1", "y2"), "clus", 2, 2,
    n_starts = 1, seed = 12, tol = 1e-12)
  theta <- coef(fit, scale = "unconstrained")
  x <- as.matrix(reference$data[, c("y1", "y2")])
  perturbations <- diag(length(theta)) * 1e-5
  objective <- function(parameters) {
    stopifnot(is.numeric(parameters))
    -.multilpa_expectation(x, fit$group_index, .multilpa_decode(parameters, fit))$log_likelihood
  }
  point <- theta + seq_along(theta) / 100
  finite_score <- vapply(seq_along(point), function(index) {
    stopifnot(is.numeric(index))
    (objective(point + perturbations[, index]) - objective(point - perturbations[, index])) / 2e-5
  }, numeric(1))
  expect_equal(.multilpa_score(point, x, fit), finite_score, tolerance = 1e-7)
  decoded <- .multilpa_decode(theta, fit)
  expect_equal(decoded$means, unname(fit$means))
  expect_equal(decoded$variances, unname(fit$variances))
  expect_equal(decoded$profile_probabilities, unname(fit$profile_probabilities))
  expect_equal(decoded$group_probabilities, unname(fit$group_probabilities))
  natural_transform <- function(parameters) {
    stopifnot(is.numeric(parameters))
    transformed_fit <- fit
    transformed_fit[names(decoded)] <- .multilpa_decode(parameters, fit)
    .multilpa_coefficients(transformed_fit, "natural")
  }
  numerical_jacobian <- vapply(seq_along(theta), function(index) {
    stopifnot(is.numeric(index))
    (natural_transform(theta + perturbations[, index]) -
       natural_transform(theta - perturbations[, index])) / 2e-5
  }, numeric(length(coef(fit))))
  expect_equal(unname(.multilpa_inference_jacobian(fit)),
               unname(numerical_jacobian), tolerance = 1e-8)
})

test_that("inference refuses data mismatches and nonregular fits", {
  set.seed(910)
  dat <- data.frame(g = rep(seq_len(10), each = 10), y = rnorm(100))
  fit <- multilpa(dat, "y", "g", 1, 1, n_starts = 1)
  # `vcov()` no longer demands data: the fit carries the columns it was built
  # from, and falling back to them gives the same answer as passing them.
  expect_equal(vcov(fit), vcov(fit, dat))
  expect_equal(vcov(fit), vcov(fit, data = get_results(fit, "data")))
  expect_error(parameter_inference(fit, dat[-1, ]),
               class = "multilpa_bad_inference_data")
  altered <- dat
  altered$y[1] <- altered$y[1] + 1
  expect_error(parameter_inference(fit, altered),
               class = "multilpa_bad_inference_data")
  if (!is.null(fit$indicator_data)) {
    altered <- dat
    altered$y[1:2] <- rev(altered$y[1:2])
    expect_error(parameter_inference(fit, altered),
                 class = "multilpa_bad_inference_data")
  }
  older_fit <- fit
  older_fit$indicator_data <- NULL
  expect_equal(parameter_inference(older_fit, dat)$standard_error,
    parameter_inference(fit, dat)$standard_error)
  expect_error(parameter_inference(fit, dat[100:1, ]),
               class = "multilpa_bad_inference_data")
  altered <- dat
  altered$y[1] <- Inf
  expect_error(parameter_inference(fit, altered),
               class = "multilpa_bad_inference_data")
  altered$y[1] <- NA_real_
  expect_error(parameter_inference(fit, altered),
               class = "multilpa_bad_inference_data")
  altered_fit <- fit
  altered_fit$converged <- FALSE
  expect_error(parameter_inference(altered_fit, dat),
               class = "multilpa_no_converge")
  altered_fit <- fit
  altered_fit$boundary <- TRUE
  expect_error(parameter_inference(altered_fit, dat),
               class = "multilpa_boundary_fit")
  expect_error(parameter_inference(fit, dat, step = 0))
  expect_error(parameter_inference(fit, dat, level = 0))
  # Identical components are an intentionally unidentified mixture.
  repeated_fit <- fit
  repeated_fit$n_profiles <- 2L
  repeated_fit$n_parameters <- 5L
  repeated_fit$means <- rbind(fit$means, fit$means)
  repeated_fit$variances <- rbind(fit$variances, fit$variances)
  repeated_fit$profile_probabilities <- matrix(c(0.5, 0.5), 1L)
  expect_error(parameter_inference(repeated_fit, dat),
               class = "multilpa_singular_information")
})

test_that("full-covariance inference matches analytic multivariate Gaussian information", {
  set.seed(954)
  dat <- data.frame(g = rep(seq_len(20), each = 20), y1 = rnorm(400), y2 = rnorm(400))
  dat$y2 <- dat$y2 + 0.8 * dat$y1
  fit <- multilpa(dat, c("y1", "y2"), "g", 1, 1,
    n_starts = 1, covariance_model = "full")
  covariance <- fit$covariances[, , 1]
  expected_covariance_se <- sqrt((covariance^2 + outer(diag(covariance), diag(covariance))) / nrow(dat))
  expected <- c(sqrt(diag(covariance) / nrow(dat)),
    expected_covariance_se[lower.tri(covariance, diag = TRUE)], 0, 0)
  information <- parameter_inference(fit, dat)
  expect_equal(unname(information$standard_error), unname(expected), tolerance = 1e-7)
  expect_equal(unname(attr(information, "covariance")[1:2, 1:2]), unname(covariance / nrow(dat)), tolerance = 1e-7)
  expect_lt(attr(information, "scaled_score"), 1e-6)
  expect_length(coef(fit, scale = "unconstrained"), fit$n_parameters)
})

test_that("diagonal FIML standard errors use each indicator's observed sample size", {
  set.seed(962)
  dat <- data.frame(g = rep(seq_len(20), each = 10), y1 = rnorm(200), y2 = rnorm(200, 2, 3))
  dat$y1[seq(1, 200, by = 5)] <- NA_real_
  dat$y2[seq(2, 200, by = 3)] <- NA_real_
  fit <- multilpa(dat, c("y1", "y2"), "g", 1, 1,
    n_starts = 1, missing = "fiml", tol = 1e-13)
  counts <- colSums(!is.na(dat[, c("y1", "y2")]))
  expected <- c(sqrt(fit$variances[1, ] / counts), sqrt(2 / counts) * fit$variances[1, ], 0, 0)
  information <- parameter_inference(fit, dat)
  expect_equal(unname(information$standard_error), unname(expected), tolerance = 1e-6)
  expect_lt(attr(information, "scaled_score"), 1e-3)
})

test_that("full-covariance score and Jacobian match numerical derivatives with missingness", {
  set.seed(481)
  group <- rep(seq_len(40), each = 12)
  group_type <- rep(seq_len(40) %% 2, each = 12)
  profile <- rbinom(length(group), 1, ifelse(group_type == 1, 0.85, 0.15))
  residual <- rnorm(length(group))
  dat <- data.frame(g = group, y1 = 5 * profile + residual,
    y2 = 4 * profile + residual / 2 + rnorm(length(group)))
  dat$y1[seq(1, nrow(dat), by = 11)] <- NA_real_
  dat$y2[seq(3, nrow(dat), by = 7)] <- NA_real_
  invisible(lapply(c("varying", "equal"), function(variance_model) {
    stopifnot(is.character(variance_model))
    fit <- multilpa(dat, c("y1", "y2"), "g", 2, 2, n_starts = 2,
      covariance_model = "full", missing = "fiml", variance_model = variance_model,
      seed = 791, tol = 1e-12)
    theta <- coef(fit, scale = "unconstrained")
    decoded <- .multilpa_decode(theta, fit)
    expect_equal(decoded$covariances, unname(fit$covariances), tolerance = 1e-12)
    expect_equal(decoded$means, unname(fit$means), tolerance = 1e-12)
    x <- as.matrix(dat[, c("y1", "y2")])
    perturbations <- diag(length(theta)) * 1e-5
    objective <- function(parameters) {
      stopifnot(is.numeric(parameters))
      -.multilpa_expectation(x, fit$group_index, .multilpa_decode(parameters, fit))$log_likelihood
    }
    point <- theta + seq_along(theta) / 1000
    finite_score <- vapply(seq_along(point), function(index) {
      stopifnot(is.numeric(index))
      (objective(point + perturbations[, index]) - objective(point - perturbations[, index])) / 2e-5
    }, numeric(1))
    expect_equal(.multilpa_score(point, x, fit), finite_score, tolerance = 1e-7)
    transform_natural <- function(parameters) {
      stopifnot(is.numeric(parameters))
      model <- fit
      model[names(decoded)] <- .multilpa_decode(parameters, model)
      .multilpa_coefficients(model, "natural")
    }
    numerical_jacobian <- vapply(seq_along(theta), function(index) {
      stopifnot(is.numeric(index))
      (transform_natural(theta + perturbations[, index]) -
        transform_natural(theta - perturbations[, index])) / 2e-5
    }, numeric(length(coef(fit))))
    expect_equal(unname(.multilpa_inference_jacobian(fit)), unname(numerical_jacobian), tolerance = 1e-8)
    information <- parameter_inference(fit, dat)
    expect_true(all(is.finite(information$standard_error)))
    expect_gt(attr(information, "condition_ratio"), 1e-5)
  }))
})

test_that("centered inference is stable for indicators with large location offsets", {
  set.seed(667)
  dat <- data.frame(g = rep(seq_len(20), each = 10), y = rnorm(200))
  fit <- multilpa(dat, "y", "g", 1, 1, n_starts = 1)
  shifted <- transform(dat, y = y + 1e9)
  shifted_fit <- multilpa(shifted, "y", "g", 1, 1, n_starts = 1)
  expect_equal(unname(parameter_inference(fit, dat)$standard_error),
    unname(parameter_inference(shifted_fit, shifted)$standard_error), tolerance = 1e-7)
  scaled <- transform(dat, y = y * 1e8)
  scaled_fit <- multilpa(scaled, "y", "g", 1, 1, n_starts = 1)
  expect_equal(unname(parameter_inference(scaled_fit, scaled)$standard_error[1:2]) / c(1e8, 1e16),
    unname(parameter_inference(fit, dat)$standard_error[1:2]), tolerance = 1e-7)
})
