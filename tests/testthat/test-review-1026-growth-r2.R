# Review 2026-10-26, round 2, growth mixtures: the fit is free of the units
# and origin of the predictors as well as of the outcome (GR2-1 to GR2-4),
# plus the class table under zero weights (GR2-5), the multilevel ICL (GR2-6),
# the held residual-SD rows (GR2-7) and the shared EM/finish bounds (GR2-8).
# A growth model with time t' = c t + o is the model with time t
# reparametrized: b' = T b with T = [[1, -o/c], [0, 1/c]], so every reported
# quantity must equal the exact image of the fit in t.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

r2_growth_data <- function(n = 90L, seed = 3L) {
  withr::with_seed(seed, {
    class <- stats::rbinom(n, 1L, 0.45) + 1L
    covariance <- list(matrix(c(1, 0.3, 0.3, 0.25), 2L),
                       matrix(c(0.6, -0.1, -0.1, 0.15), 2L))
    effects <- t(vapply(seq_len(n), function(i) {
      as.vector(t(chol(covariance[[class[i]]])) %*% stats::rnorm(2L))
    }, numeric(2L)))
    out <- data.frame(id = rep(seq_len(n), each = 5L), time = rep(0:4, n))
    out$x <- stats::rnorm(n)[out$id]
    person_class <- class[out$id]
    out$y <- c(2, 4.5)[person_class] + c(0.5, -0.3)[person_class] * out$time +
      effects[out$id, 1L] + effects[out$id, 2L] * out$time +
      stats::rnorm(nrow(out), 0, c(0.7, 0.5)[person_class])
    out[stats::runif(nrow(out)) > 0.15, ]
  })
}

r2_fit <- function(data, formula = y ~ time, random = ~ 1 + time, ...) {
  suppressWarnings(mixture_regression(formula, data, n_classes = 2L, id = "id",
                                      class_level = "group", random = random,
                                      membership = ~ x, n_starts = 2L, seed = 7L, ...))
}

# The coefficient-and-random-effect map of time t' = c t + o, applied to the
# estimates of the fit in t, with standard errors by central differences
# (steps of 1e-4 standard errors) of that map.
r2_expected <- function(base, c, o) {
  inference <- latents:::.growth_resolve_inference(base)
  theta <- inference$theta
  vcov <- inference$vcov
  spec <- base$spec
  transform <- matrix(c(1, 0, -o / c, 1 / c), 2L)
  natural <- function(v) {
    params <- latents:::.growth_unpack(spec, v, base$params)
    coefficients <- unlist(lapply(seq_len(spec$n_classes), function(k) {
      beta <- params$beta[, k]
      c(beta[1L] - beta[2L] * o / c, beta[2L] / c)
    }))
    random <- unlist(lapply(latents:::.growth_covariances(spec, params), function(g) {
      mapped <- transform %*% g %*% t(transform)
      c(diag(mapped), sqrt(diag(mapped)), mapped[1L, 2L] / sqrt(prod(diag(mapped))))
    }))
    c(coefficients, random)
  }
  estimate <- natural(theta)
  se_theta <- sqrt(diag(vcov))
  jacobian <- vapply(seq_along(theta), function(j) {
    h <- 1e-4 * se_theta[j]
    (natural(replace(theta, j, theta[j] + h)) -
       natural(replace(theta, j, theta[j] - h))) / (2 * h)
  }, numeric(length(estimate)))
  list(estimate = estimate,
       std_error = sqrt(rowSums((jacobian %*% vcov) * jacobian)))
}

# Every table of the fit in t' against the image of the fit in t, in standard
# errors; flags identical.
r2_expect_mapped <- function(base, moved, c, o, tolerance = 1e-6) {
  expected <- r2_expected(base, c, o)
  coefficients <- as.data.frame(get_results(moved, "coefficients"))
  random <- as.data.frame(get_results(moved, "random"))
  # Random rows come in the order variance (2), sd (2), correlation, by class.
  obtained <- c(coefficients$estimate, random$estimate)
  obtained_se <- c(coefficients$std_error, random$std_error)
  expect_lt(max(abs(obtained - expected$estimate) / expected$std_error), tolerance)
  # A correlation within 1e-4 of -1 (the calendar-year chart, r = -0.9999999)
  # has a delta standard error that numerical differences determine only to
  # about 1e-4 (central differences and numDeriv::grad of the same map differ
  # by 3.4e-5; the package's own steps by up to 2.2e-4). Its estimate is
  # still held to 1e-6 standard errors above.
  near_one <- c(rep(FALSE, nrow(coefficients)),
                random$parameter == "correlation" & abs(random$estimate) > 0.9999)
  expect_lt(max(abs(obtained_se / expected$std_error - 1)[!near_one]), tolerance)
  expect_lt(max(c(0, abs(obtained_se / expected$std_error - 1)[near_one])), 1e-3)
  same <- function(what, columns) {
    a <- as.data.frame(get_results(base, what))[columns]
    b <- as.data.frame(get_results(moved, what))[columns]
    expect_equal(b, a, tolerance = tolerance)
  }
  same("classes", c("share", "share_std_error", "expected_persons", "assigned_persons",
                    "average_posterior", "residual_sd", "residual_sd_std_error"))
  same("membership", c("estimate", "std_error"))
  same("fit", c("log_likelihood", "bic", "icl", "entropy", "converged", "random_boundary"))
  same("assignments", c("posterior_class_1", "posterior_class_2"))
  same("trajectories", c("estimate", "std_error", "conf_low", "conf_high"))
  same("individual", c("class_mean", "predicted"))
  expect_equal(moved$starts$log_likelihood, base$starts$log_likelihood, tolerance = 1e-9)
  effects <- as.data.frame(get_results(moved, "random_effects"))
  base_effects <- as.matrix(as.data.frame(get_results(base, "random_effects"))[3:4])
  transform <- matrix(c(1, 0, -o / c, 1 / c), 2L)
  expect_equal(unname(as.matrix(effects[3:4])),
               unname(base_effects %*% t(transform)), tolerance = 1e-6)
  expect_identical(moved$converged, base$converged)
  expect_identical(moved$degenerate_classes, base$degenerate_classes)
  expect_identical(moved$random_boundary, base$random_boundary)
}

test_that("time in days fits as time in years (GR2-1)", {
  data <- r2_growth_data()
  years <- r2_fit(data)
  days <- data
  days$time <- days$time * 365
  # Before: the slope variance was held at the variance floor, 11 units of
  # log likelihood low, reported converged and without a boundary warning.
  r2_expect_mapped(years, r2_fit(days), 365, 0)
})

test_that("time in hours and in units of 1e5 fit as time in years (GR2-3, GR2-4)", {
  data <- r2_growth_data()
  years <- r2_fit(data)
  hours <- data
  hours$time <- hours$time * 8760
  r2_expect_mapped(years, r2_fit(hours), 8760, 0)
  large <- data
  large$time <- large$time * 1e5
  r2_expect_mapped(years, r2_fit(large), 1e5, 0)
})

test_that("a calendar-year random slope fits as time from zero (GR2-2)", {
  data <- r2_growth_data()
  base <- r2_fit(data)
  calendar <- data
  calendar$time <- calendar$time + 2000
  # Before: "EM likelihood decreased ... This is a defect".
  r2_expect_mapped(base, r2_fit(calendar), 1, 2000)
})

test_that("a covariance with an ill-conditioned chart is factorized exactly (GR2-2)", {
  # A calendar-year random-effect covariance: variances 1.6e6 and 0.39, smallest
  # eigenvalue about 5e-7. A jitter of 1e-12 of the largest variance was three
  # times that eigenvalue.
  transform <- matrix(c(1, 0, -2000, 1), 2L)
  covariance <- transform %*% matrix(c(0.71, 0.28, 0.28, 0.39), 2L) %*% t(transform)
  factor <- latents:::.growth_factor(covariance)
  expect_equal(factor, t(chol(covariance)), tolerance = 1e-12)
  expect_equal(tcrossprod(factor), covariance, tolerance = 1e-12)
})

test_that("the class table counts persons of positive weight only (GR2-5)", {
  data <- r2_growth_data()
  zero <- c(3L, 10L, 44L, 45L, 70L)
  data$w <- withr::with_seed(2L, stats::runif(90L, 0.5, 2))[data$id]
  data$w[data$id %in% zero] <- 0
  fit_args <- function(d) suppressWarnings(mixture_regression(
    y ~ time, d, n_classes = 2L, id = "id", class_level = "group",
    random = ~ 1 + time, weights = "w", n_starts = 2L, seed = 4L))
  all <- as.data.frame(get_results(fit_args(data), "classes"))
  dropped <- as.data.frame(get_results(fit_args(data[!data$id %in% zero, ]), "classes"))
  expect_identical(all$assigned_persons, dropped$assigned_persons)
  expect_equal(all$average_posterior, dropped$average_posterior, tolerance = 1e-6)
  expect_equal(all$expected_persons, dropped$expected_persons, tolerance = 1e-6)
})

test_that("the multilevel ICL charges the joint classification entropy (GR2-6)", {
  data <- r2_growth_data(n = 120L, seed = 13L)
  data$school <- (data$id - 1L) %/% 6L + 1L
  fit <- suppressWarnings(mixture_regression(
    y ~ time, data, n_classes = 2L, id = "id", class_level = "group",
    random = ~ 1 + time, cluster = "school", n_group_classes = 2L,
    n_starts = 2L, seed = 4L))
  # By hand from the person log densities: rho_h prod_i p(k_i | h, y_i).
  softmax <- function(v) exp(v - max(v)) / sum(exp(v - max(v)))
  within <- t(apply(fit$params$class_logits[1:2, , drop = FALSE], 1L, softmax))
  shares <- softmax(fit$params$delta[1L, ])
  log_density <- fit$expectation$log_density
  cluster <- fit$spec$cluster_index
  entropy <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  joint <- sum(vapply(seq_len(fit$spec$n_clusters), function(j) {
    members <- which(cluster == j)
    conditional <- lapply(1:2, function(h) {
      t(vapply(members, function(i) softmax(log(within[h, ]) + log_density[i, ]),
               numeric(2L)))
    })
    log_marginal <- vapply(1:2, function(h) {
      log(shares[h]) + sum(vapply(members, function(i) {
        v <- log(within[h, ]) + log_density[i, ]
        max(v) + log(sum(exp(v - max(v))))
      }, numeric(1)))
    }, numeric(1))
    rho <- softmax(log_marginal)
    entropy(rho) + sum(rho * vapply(conditional, entropy, numeric(1)))
  }, numeric(1)))
  table <- as.data.frame(get_results(fit, "fit"))
  expect_equal(table$icl, table$bic + 2 * joint, tolerance = 1e-8)
})

test_that("only the floored class's residual SD loses its interval (GR2-7)", {
  fit <- r2_fit(r2_growth_data())
  fit$residual_boundary <- "log_sigma.class_2"
  natural <- latents:::.growth_natural_parameters(fit)
  held <- natural$degenerate[natural$labels$block == "residual_sd"]
  expect_identical(unname(held), c(FALSE, TRUE))
})

test_that("EM floors conditional variances, as the finish bounds them (GR2-8)", {
  covariance <- matrix(c(1, 0.99999, 0.99999, 1), 2L)
  floored <- latents:::.growth_floor_covariance(covariance, 1e-3)
  expect_gte(min(diag(t(chol(floored)))^2), 1e-3 * (1 - 1e-12))
  expect_equal(floored[1L, ], covariance[1L, ])
  # Untouched above the floor.
  expect_identical(latents:::.growth_floor_covariance(diag(2), 1e-3), diag(2))
  # A fit's standard parameters satisfy the finish's bounds.
  fit <- r2_fit(r2_growth_data())
  floor <- fit$spec$min_variance / fit$standard$scale^2
  conditional <- unlist(lapply(fit$standard$params$random_covariance, function(g) {
    diag(t(chol(g)))^2
  }))
  expect_true(all(conditional >= floor * (1 - 1e-12)))
})
