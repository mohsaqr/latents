test_that("Gaussian missing-pattern moments agree with conditional normal formulas", {
  x <- rbind(c(3, NA), c(NA, 5), c(NA, NA), c(3, 5))
  parameters <- list(means = matrix(c(1, 2), 1L),
                     variances = matrix(c(4, 9), 1L),
                     covariances = array(c(4, 1, 1, 9), c(2L, 2L, 1L)))
  moments <- .multilpa_gaussian_moments(x, parameters)
  expect_equal(dim(moments$log_density), c(4L, 1L))
  expect_equal(moments$moments[[1L]]$expected,
               rbind(c(3, 2.5), c(4 / 3, 5), c(1, 2), c(3, 5)))
  adjustment <- moments$moments[[1L]]$adjustments
  expect_equal(adjustment$TRUEFALSE$covariance[2L, 2L], 8.75)
  expect_equal(adjustment$FALSETRUE$covariance[1L, 1L], 4 - 1 / 9)
  expect_equal(adjustment$FALSEFALSE$covariance, matrix(c(4, 1, 1, 9), 2L))
  covariance <- matrix(c(4, 1, 1, 9), 2L)
  complete_log_density <- -log(2 * pi) - log(det(covariance)) / 2 -
    as.numeric(crossprod(c(2, 3), solve(covariance, c(2, 3)))) / 2
  expect_equal(as.numeric(moments$log_density),
               c(dnorm(3, 1, 2, log = TRUE), dnorm(5, 2, 3, log = TRUE),
                 0, complete_log_density), tolerance = 1e-12)
  diagonal <- parameters
  diagonal$covariances <- NULL
  independent <- .multilpa_gaussian_moments(x, diagonal)
  expect_equal(independent$moments[[1L]]$expected,
               rbind(c(3, 2), c(1, 5), c(1, 2), c(3, 5)))
  expect_error(.multilpa_gaussian_moments(matrix(Inf, 1L), parameters))
})

test_that("diagonal FIML single-profile estimates equal available-case Gaussian ML", {
  set.seed(244)
  data <- data.frame(group = rep(seq_len(20), each = 10),
                     x = rnorm(200, 2, 1.4), y = rnorm(200, -1, 0.7))
  data$x[c(1:35, 72)] <- NA_real_
  data$y[c(1, 91:140)] <- NA_real_
  fit <- multilpa(data, c("x", "y"), "group", 1, 1,
                    missing = "fiml", n_starts = 1, seed = 8, tol = 1e-13)
  means <- vapply(data[c("x", "y")], mean, numeric(1), na.rm = TRUE)
  variances <- vapply(data[c("x", "y")], function(values) {
    mean((values - mean(values, na.rm = TRUE))^2, na.rm = TRUE)
  }, numeric(1))
  expect_equal(as.numeric(fit$means), unname(means), tolerance = 1e-6)
  expect_equal(as.numeric(fit$variances), unname(variances), tolerance = 1e-6)
  expected_log_likelihood <- sum(dnorm(data$x, means["x"], sqrt(variances["x"]),
                                     log = TRUE), na.rm = TRUE) +
    sum(dnorm(data$y, means["y"], sqrt(variances["y"]), log = TRUE), na.rm = TRUE)
  expect_equal(fit$log_likelihood, expected_log_likelihood, tolerance = 1e-10)
  expect_equal(fit$n_observations, nrow(data))
  expect_equal(fit$n_observed_by_indicator, c(x = 164L, y = 149L))
  expect_true(all(diff(fit$log_likelihood_history) >= -1e-9))
  expect_equal(fit$subject_posteriors[1L, 1L], 1)
  expect_error(multilpa(data, c("x", "y"), "group", 1, 1), "missing")
})

test_that("full covariance complete-data Gaussian ML has analytic means and covariance", {
  set.seed(415)
  data <- data.frame(group = rep(seq_len(10), each = 12), x = rnorm(120))
  data$y <- 1.5 + 0.8 * data$x + rnorm(120, sd = .6)
  fit <- multilpa(data, c("x", "y"), "group", 1, 1,
                    covariance_model = "full", n_starts = 1)
  x <- as.matrix(data[c("x", "y")])
  expect_equal(as.numeric(fit$means), colMeans(x), ignore_attr = TRUE)
  expected_covariance <- crossprod(sweep(x, 2, colMeans(x), "-")) / nrow(x)
  expect_equal(unname(fit$covariances[, , 1L]), unname(expected_covariance))
  expect_equal(as.numeric(fit$variances), diag(expected_covariance), ignore_attr = TRUE)
  expect_equal(fit$n_parameters, 5)
  expect_true(fit$converged)
  expect_false(fit$boundary)
  one <- multilpa(data, "x", "group", 1, 1,
                    covariance_model = "full", n_starts = 1)
  expect_equal(dim(one$covariances), c(1L, 1L, 1L))
  expect_equal(as.numeric(one$covariances), mean((data$x - mean(data$x))^2))
})

test_that("full covariance FIML agrees with independent direct Gaussian likelihood", {
  set.seed(148)
  data <- data.frame(group = rep(seq_len(20), each = 10), x = rnorm(200))
  data$y <- 2 + .75 * data$x + rnorm(200, sd = .8)
  data$x[1:40] <- NA_real_
  data$y[c(1:4, 80:120)] <- NA_real_
  fit <- multilpa(data, c("x", "y"), "group", 1, 1,
                    covariance_model = "full", missing = "fiml", n_starts = 1,
                    tol = 1e-13)
  # Independent bivariate density plus univariate marginals; no package helpers.
  objective <- function(theta) {
    stopifnot(is.numeric(theta), length(theta) == 5L)
    sx <- exp(theta[3]); sy <- exp(theta[4]); rho <- tanh(theta[5])
    complete <- !is.na(data$x) & !is.na(data$y)
    x_only <- !is.na(data$x) & is.na(data$y)
    y_only <- is.na(data$x) & !is.na(data$y)
    zx <- (data$x[complete] - theta[1]) / sx
    zy <- (data$y[complete] - theta[2]) / sy
    complete_log <- -log(2 * pi * sx * sy) - log(1 - rho^2) / 2 -
      (zx^2 + zy^2 - 2 * rho * zx * zy) / (2 * (1 - rho^2))
    -sum(complete_log) - sum(dnorm(data$x[x_only], theta[1], sx, log = TRUE)) -
      sum(dnorm(data$y[y_only], theta[2], sy, log = TRUE))
  }
  reference <- optim(c(0, 2, 0, 0, 0), objective, method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 2000))
  expect_equal(reference$convergence, 0L)
  expect_equal(fit$log_likelihood, -reference$value, tolerance = 1e-9)
  expect_equal(as.numeric(fit$means), reference$par[1:2], tolerance = 1e-5)
  expected_covariance <- matrix(c(exp(2 * reference$par[3]),
    tanh(reference$par[5]) * exp(sum(reference$par[3:4])),
    tanh(reference$par[5]) * exp(sum(reference$par[3:4])),
    exp(2 * reference$par[4])), 2L)
  expect_equal(unname(fit$covariances[, , 1]), expected_covariance, tolerance = 1e-5)
  expect_true(all(diff(fit$log_likelihood_history) >= -1e-9))
})

test_that("FIML multilevel rows use their observed group context", {
  set.seed(123)
  group <- rep(seq_len(40), each = 12)
  group_type <- rep(c(1, 2), each = 20)
  profile <- 1 + (runif(length(group)) > c(.9, .15)[group_type[group]])
  data <- data.frame(group = group, x = rnorm(length(group), c(-3, 3)[profile]))
  data$y <- .4 * data$x + rnorm(nrow(data), c(-1, 2)[profile], .8)
  data$x[seq(1, nrow(data), 7)] <- NA
  data$y[seq(1, nrow(data), 9)] <- NA
  invisible(lapply(c("equal", "varying"), function(variance_model) {
    fit <- multilpa(data, c("x", "y"), "group", 2, 2, missing = "fiml",
                      covariance_model = "full", variance_model = variance_model,
                      n_starts = 3, seed = 23, tol = 1e-10)
    expect_true(fit$converged)
    expect_true(all(diff(fit$log_likelihood_history) >= -1e-8))
    expected <- as.numeric(fit$group_posteriors[fit$group_index[1], ] %*%
                            fit$profile_probabilities)
    expect_equal(as.numeric(fit$subject_posteriors[1, ]), expected, tolerance = 1e-12)
    expect_equal(unname(rowSums(fit$subject_posteriors)), rep(1, nrow(data)))
    expect_equal(unname(rowSums(fit$group_posteriors)), rep(1, 40))
    if (variance_model == "equal") {
      expect_equal(unname(fit$covariances[, , 1]), unname(fit$covariances[, , 2]))
    }
  }))
})

test_that("covariance eigenvalue bounds and invalid missing-data inputs are explicit", {
  covariance <- matrix(c(1, .99, .99, 1), 2)
  bounded <- .multilpa_bound_covariance(covariance, .1)
  expect_equal(eigen(bounded, symmetric = TRUE)$values, c(1.99, .1))
  expect_equal(.multilpa_bound_covariance(covariance, .001), covariance)
  expect_true(.multilpa_covariance_boundary(list(means = matrix(0, 1, 2),
    covariances = array(bounded, c(2, 2, 1))), .1))
  expect_error(.multilpa_bound_covariance(matrix(NA_real_, 1), .1))
  data <- data.frame(group = rep(1:5, each = 4), x = seq_len(20), y = seq_len(20))
  expect_warning(fit <- multilpa(data, c("x", "y"), "group", 1, 1,
    n_starts = 1, covariance_model = "full", min_variance = .1), "bound-active")
  expect_true(fit$boundary)
  expect_true(fit$starts$boundary)
  data$x <- NA_real_
  expect_error(multilpa(data, c("x", "y"), "group", 1, 1,
                          missing = "fiml"), "observed values")
  data$x <- c(Inf, rep(NA_real_, 19))
  expect_error(multilpa(data, c("x", "y"), "group", 1, 1,
                          missing = "fiml"), "non-finite")
  data$x <- c(1, rep(NA_real_, 19))
  expect_error(multilpa(data, c("x", "y"), "group", 1, 1,
                          missing = "fiml"), "Constant")
  start <- list(means = matrix(c(0, 0), 1), covariances = array(diag(2), c(2, 2, 1)),
                profile_probabilities = matrix(1), group_probabilities = 1)
  valid <- .multilpa_validate_start(start, 1, 1, 2, "varying", .1, "full")
  expect_equal(valid$variances, matrix(c(1, 1), 1))
  start$covariances[1, 2, 1] <- .2
  expect_error(.multilpa_validate_start(start, 1, 1, 2, "varying", .1, "full"), "symmetric")
  start$covariances <- array(diag(2), c(2, 2, 1))
  start$variances <- matrix(c(2, 2), 1)
  expect_error(.multilpa_validate_start(start, 1, 1, 2, "varying", .1, "full"), "diagonals")
  start$covariances <- diag(2)
  expect_error(.multilpa_validate_start(start, 1, 1, 2, "varying", .1, "full"), "array")
})
test_that("FIML rejects NaN consistently with inference", {
  data <- data.frame(group = rep(1:3, each = 2), y = c(1, 2, NaN, 4, 5, 6))
  expect_error(multilpa(data, "y", "group", 1, 1, missing = "fiml"), "non-finite")
})
