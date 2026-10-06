
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
