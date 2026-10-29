# Comparisons with independent references moved from tests/testthat/test-review-1026-inference.R,
# so that they run in CI's equivalence job rather than under R CMD check.

review_two_profile_data <- function(seed = 11L, n_groups = 30L, per_group = 6L) {
  set.seed(seed)
  data <- data.frame(g = rep(seq_len(n_groups), each = per_group))
  membership <- stats::rbinom(nrow(data), 1L, 0.4)
  data$a <- stats::rnorm(nrow(data), ifelse(membership == 1L, 2, 0),
                         ifelse(membership == 1L, 0.7, 1))
  data$b <- stats::rnorm(nrow(data), ifelse(membership == 1L, -1, 0.5), 1)
  data
}

# Independent observed-data log likelihood of a two-profile, one-group-class
# diagonal (VVI) mixture in the package's unconstrained order: means by
# profile, log variances by profile, the profile-1 logit against profile 2.
review_vvi_loglik <- function(theta, x, rows = seq_len(nrow(x))) {
  means <- matrix(theta[1:4], 2L, byrow = TRUE)
  log_variances <- matrix(theta[5:8], 2L, byrow = TRUE)
  weights <- exp(c(theta[9], 0) - max(theta[9], 0))
  weights <- weights / sum(weights)
  density <- vapply(1:2, function(k) {
    weights[k] *
      stats::dnorm(x[rows, 1L], means[k, 1L], sqrt(exp(log_variances[k, 1L]))) *
      stats::dnorm(x[rows, 2L], means[k, 2L], sqrt(exp(log_variances[k, 2L])))
  }, numeric(length(rows)))
  sum(log(rowSums(matrix(density, ncol = 2L))))
}

test_that("observed, robust (CR1) and OPG covariances match independent numDeriv computations", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  data <- review_two_profile_data()
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          variance_model = "varying", n_starts = 3, seed = 1, tol = 1e-10))
  x <- as.matrix(data[, c("a", "b")])
  theta <- unname(coef(fit, scale = "unconstrained"))
  expect_equal(review_vvi_loglik(theta, x), fit$log_likelihood, tolerance = 1e-10)

  observed <- parameter_inference(fit)
  hessian <- numDeriv::hessian(function(t) -review_vvi_loglik(t, x), theta)
  expect_lt(max(abs(hessian - unname(attr(observed, "hessian")))) / max(abs(hessian)), 1e-6)
  covariance <- solve(hessian)
  natural <- function(t) {
    weights <- exp(c(t[9], 0)) / sum(exp(c(t[9], 0)))
    c(t[1:4], exp(t[5:8]), weights, 1)
  }
  jacobian <- numDeriv::jacobian(natural, theta)
  delta_se <- sqrt(diag(jacobian %*% covariance %*% t(jacobian)))
  expect_equal(observed$standard_error, delta_se, tolerance = 1e-6)

  # Per-group scores of the independent likelihood, then the CR1 sandwich
  # G / (G - 1) * A^-1 B A^-1 (Cameron and Miller 2015) and the OPG B^-1.
  n_groups <- max(data$g)
  scores <- t(vapply(seq_len(n_groups), function(group) {
    numDeriv::grad(function(t) review_vvi_loglik(t, x, which(data$g == group)), theta)
  }, numeric(length(theta))))
  meat <- crossprod(scores)
  robust <- parameter_inference(fit, vcov_type = "robust")
  sandwich_covariance <- covariance %*% meat %*% covariance * n_groups / (n_groups - 1)
  expect_lt(max(abs(sandwich_covariance - unname(attr(robust, "covariance_unconstrained")))) /
              max(abs(sandwich_covariance)), 1e-6)
  expect_equal(attr(robust, "scaling_correction"),
               sum(diag(covariance %*% meat)) / length(theta), tolerance = 1e-6)
  opg <- parameter_inference(fit, vcov_type = "opg")
  expect_lt(max(abs(solve(meat) - unname(attr(opg, "covariance_unconstrained")))) /
              max(abs(solve(meat))), 1e-6)

  # Intervals: symmetric for means, log scale for variances, logit scale for
  # interior probabilities.
  critical <- stats::qnorm(0.975)
  mean_rows <- observed$parameter == "mean"
  expect_equal(observed$conf_low[mean_rows],
               observed$estimate[mean_rows] - critical * observed$standard_error[mean_rows])
  variance_rows <- observed$parameter == "variance"
  expect_equal(observed$conf_low[variance_rows],
               observed$estimate[variance_rows] *
                 exp(-critical * observed$standard_error[variance_rows] /
                       observed$estimate[variance_rows]))
  interior <- observed$parameter == "probability" & observed$estimate < 1
  p <- observed$estimate[interior]
  expect_equal(observed$conf_high[interior],
               stats::plogis(stats::qlogis(p) +
                               critical * observed$standard_error[interior] / (p * (1 - p))))
})

test_that("boundary = 'fix' inverts the Hessian of the free coordinates only", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  set.seed(31)
  n <- 200L
  state <- sample.int(2L, n, TRUE)
  data <- data.frame(g = rep(1:40, each = 5), y = c(-2, 2)[state] + stats::rnorm(n))
  data$c <- factor(ifelse(state == 1L, sample(c("b", "c"), n, TRUE),
                          sample(c("a", "b", "c"), n, TRUE)))
  fit <- quietly(multilpa(data, c("y", "c"), "g", n_profiles = 2, n_group_classes = 1,
                          categorical = "c", n_starts = 3, seed = 1, tol = 1e-12,
                          max_iter = 5000))
  expect_error(parameter_inference(fit), class = "latents_boundary_fit")
  held <- parameter_inference(fit, boundary = "fix")
  expect_identical(attr(held, "fixed_at_bound"), "measurement.response.profile_1.c:a")

  theta <- unname(coef(fit, scale = "unconstrained"))
  free <- which(!grepl("profile_1\\.c:a$", names(coef(fit, scale = "unconstrained"))))
  code <- match(as.character(data$c), colnames(fit$response_probabilities$c))
  loglik <- function(t) {
    response <- t(apply(cbind(matrix(t[5:8], 2L, byrow = TRUE), 0), 1L,
                        function(z) exp(z - max(z)) / sum(exp(z - max(z)))))
    weights <- exp(c(t[9], 0)) / sum(exp(c(t[9], 0)))
    sum(log(rowSums(vapply(1:2, function(k) {
      weights[k] * stats::dnorm(data$y, t[k], sqrt(exp(t[2 + k]))) * response[k, code]
    }, numeric(n)))))
  }
  expect_equal(loglik(theta), fit$log_likelihood, tolerance = 1e-10)
  hessian <- numDeriv::hessian(function(t) {
    full <- theta
    full[free] <- t
    -loglik(full)
  }, theta[free])
  covariance <- solve(hessian)
  expect_lt(max(abs(covariance - unname(attr(held, "covariance_unconstrained")))) /
              max(abs(covariance)), 1e-5)
})

test_that("Rubin's rules match hand computation and mice::pool.scalar", {
  set.seed(4)
  m <- 5L
  estimates <- matrix(stats::rnorm(3 * m, c(1, 2, 0.5), 0.2), 3L)
  variances <- matrix(stats::runif(3 * m, 0.01, 0.05), 3L)
  estimates[3, ] <- 1
  variances[3, ] <- 0
  tables <- lapply(seq_len(m), function(i) {
    data.frame(level = "measurement", outcome = "profile_1", term = c("a", "b", "c"),
               parameter = "mean", estimate = estimates[, i],
               standard_error = sqrt(variances[, i]), statistic = c(1, NA, NA))
  })
  pooled <- .multilpa_rubin(tables, 0.95)
  q_bar <- mean(estimates[1, ])
  u_bar <- mean(variances[1, ])
  b <- stats::var(estimates[1, ])
  total <- u_bar + (1 + 1 / m) * b
  riv <- (1 + 1 / m) * b / u_bar
  df <- (m - 1) * (1 + 1 / riv)^2
  expect_equal(pooled$estimate[1], q_bar)
  expect_equal(pooled$standard_error[1], sqrt(total))
  expect_equal(pooled$df[1], df)
  expect_equal(pooled$fmi[1], (riv + 2 / (df + 3)) / (riv + 1))
  expect_equal(pooled$conf_low[1], q_bar - stats::qt(0.975, df) * sqrt(total))
  expect_equal(pooled$p_value[1], 2 * stats::pt(-abs(q_bar / sqrt(total)), df))
  expect_true(is.na(pooled$p_value[2]))
  # A structural constant: no variance, no missing information, infinite df.
  expect_identical(c(pooled$riv[3], pooled$fmi[3], pooled$df[3]), c(0, 0, Inf))
  skip_if_not_installed("mice")
  reference <- lapply(1:2, function(j) {
    mice::pool.scalar(estimates[j, ], variances[j, ], n = Inf, k = 1)
  })
  expect_equal(pooled$standard_error[1:2]^2, vapply(reference, `[[`, numeric(1), "t"))
  expect_equal(pooled$df[1:2], vapply(reference, `[[`, numeric(1), "df"))
  expect_equal(pooled$riv[1:2], vapply(reference, `[[`, numeric(1), "r"))
  expect_equal(pooled$fmi[1:2], vapply(reference, `[[`, numeric(1), "fmi"))
})
