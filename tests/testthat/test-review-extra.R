# Regression tests from past audits; CI runs them on every platform.
skip_on_cran()

test_that("extra indicator columns are finite category/count vectors", {
  invisible(lapply(list(c(1, 2, Inf), c(1, 2, -Inf), matrix(1:3, 3)), function(value) {
    data <- data.frame(o = I(value))
    expect_error(.latents_prepare_extra(data, ordinal = "o"), class = "latents_bad_data")
  }))
  data <- data.frame(k = I(matrix(0:2, 3)))
  expect_error(.latents_prepare_extra(data, count = "k"), class = "latents_bad_data")
  object <- list(extra_data = .latents_prepare_extra(data.frame(o = 1:3, k = 0:2), "o", "k"))
  expect_error(.latents_prepare_extra_like(object, data.frame(o = 1:3, k = I(matrix(0:2, 3)))),
               class = "latents_bad_data")
  expect_error(.latents_prepare_extra_like(object, data.frame(o = I(matrix(1:3, 3)), k = 0:2)),
               class = "latents_bad_data")
})

test_that("NB scores and curvature stay accurate near the Poisson limit", {
  y <- 0:30
  means <- c(0.1, 3, 20)
  invisible(lapply(c(1e-8, 1e-6, 1e-5, 1e-4, 0.01, 0.2), function(alpha) {
    actual <- .latents_negative_binomial_dispersion_derivatives(y, means, rep(alpha, 3))
    reference <- lapply(means, function(mu) {
      t(vapply(y, function(count) {
        j <- seq_len(count) - 1
        x <- alpha * mu
        difference <- log1p(x) / alpha - mu / (1 + x)
        c(score = sum(j * alpha / (1 + j * alpha)) + difference - count * x / (1 + x),
          curvature = sum(j * alpha / (1 + j * alpha)^2) - difference +
            mu * x / (1 + x)^2 - count * x / (1 + x)^2)
      }, numeric(2)))
    })
    expect_lt(max(abs(actual$score - do.call(cbind, lapply(reference, function(x) x[, 1])))), 1e-9)
    expect_lt(max(abs(actual$curvature - do.call(cbind, lapply(reference, function(x) x[, 2])))), 1e-9)
    if (alpha == 1e-8) {
      expect_lt(max(abs(actual$score - do.call(cbind, lapply(reference, function(x) x[, 1])))), 1e-13)
      expect_lt(max(abs(actual$curvature - do.call(cbind, lapply(reference, function(x) x[, 2])))), 1e-13)
    }
  }))
})

test_that("covariate ordinal inference refuses changed category identities", {
  set.seed(876)
  n <- 300
  x <- rnorm(n)
  profile <- rbinom(n, 1, plogis(x)) + 1L
  data <- data.frame(x = x, y = rnorm(n, c(-2, 2)[profile], 0.8),
                     o = sample(1:3, n, replace = TRUE))
  fit <- quietly(lpa(data, c("y", "o"), 2, ordinal = "o", profile_covariates = "x",
                     n_starts = 1, seed = 1, tol = 1e-12))
  changed <- data
  changed$o <- 10 * changed$o
  expect_error(vcov(fit, changed), class = "latents_bad_inference_data")
  expect_equal(vcov(fit), vcov(fit, data), tolerance = 1e-12)
})

test_that("weighted mixed NB covariate scores and natural maps match independent derivatives", {
  set.seed(878)
  n <- 400L
  group <- rep(seq_len(40L), each = 10L)
  group_covariate <- rnorm(40L)
  group_class <- rbinom(40L, 1, plogis(group_covariate))
  x <- rnorm(n)
  profile <- rbinom(n, 1, plogis(-1.5 + 3 * group_class[group] + 0.6 * x)) + 1L
  data <- data.frame(g = group, x = x, v = group_covariate[group],
                     sw = rep(runif(40L, 0.5, 2), each = 10L),
                     y1 = rnorm(n, c(-2, 2)[profile], 0.8),
                     y2 = rnorm(n, c(-1, 1)[profile], 0.8),
                     o = sample(1:3, n, replace = TRUE),
                     k = rnbinom(n, mu = c(2, 8)[profile], size = 2))
  data$o[5] <- NA
  data$k[9] <- NA
  data$y2[11] <- NA
  invisible(lapply(c("equal", "varying"), function(dispersion) {
    fit <- quietly(multilpa(data, c("y1", "y2", "o", "k"), "g", 2, 2,
                            profile_covariates = "x", group_covariates = "v",
                            profile_slopes = "group_class", ordinal = "o", count = "k",
                            count_model = "negative_binomial", count_dispersion = dispersion,
                            covariance_model = "full", weights = "sw", missing = "fiml",
                            n_starts = 1, seed = 1, tol = 1e-10))
    theta <- .multilpa_cov_encode(fit) + rep_len(c(-0.03, 0.02), fit$n_parameters)
    indicators <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
    objective <- function(point) {
      pieces <- .multilpa_cov_decode(point, fit)
      .multilpa_cov_expectation(indicators, fit$group_index, pieces$parameters,
                                fit$profile_design, fit$group_design, pieces$beta,
                                pieces$gamma, sampling_weights = unname(fit$sampling_weights),
                                extra = fit$extra_data)$log_likelihood
    }
    gradient <- vapply(seq_along(theta), function(index) {
      up <- down <- theta
      up[index] <- up[index] + 1e-5
      down[index] <- down[index] - 1e-5
      (objective(up) - objective(down)) / 2e-5
    }, numeric(1))
    expect_lt(max(abs(colSums(.multilpa_cov_group_scores(theta, indicators, fit)) - gradient)), 1e-5)
    labels <- .multilpa_cov_labels(fit)
    natural <- function(point) .multilpa_cov_natural_estimate(fit, point, labels)
    numerical_jacobian <- vapply(seq_along(theta), function(index) {
      up <- down <- theta
      up[index] <- up[index] + 1e-5
      down[index] <- down[index] - 1e-5
      (natural(up) - natural(down)) / 2e-5
    }, numeric(nrow(labels)))
    expect_equal(.multilpa_cov_natural_jacobian(fit, theta, labels), numerical_jacobian,
                 tolerance = 1e-7)
    expect_identical(length(theta), as.integer(fit$n_parameters))
  }))
})

test_that("Poisson log densities retain precision for large count means", {
  extra <- .latents_prepare_extra(data.frame(k = c(1e8, 1e15, NA_real_)), count = "k", missing = "fiml")
  parameters <- list(count_means = matrix(c(1e8, 1e15), 2L))
  actual <- .latents_extra_log_density(extra, parameters, 3L, 2L)
  expected <- rbind(c(dpois(1e8, 1e8, log = TRUE), dpois(1e8, 1e15, log = TRUE)),
                     c(dpois(1e15, 1e8, log = TRUE), dpois(1e15, 1e15, log = TRUE)), c(0, 0))
  expect_identical(actual, expected)
})
