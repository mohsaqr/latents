# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-review-1026-growth.R,
# so that they run in CI's equivalence job rather than under R CMD check.

review_growth_data <- function(n = 70L, seed = 21L, classes = 2L) {
  withr::with_seed(seed, {
    times <- 0:4
    class <- sample(seq_len(classes), n, TRUE)
    effects <- t(t(chol(matrix(c(1, 0.2, 0.2, 0.15), 2L))) %*%
                   matrix(stats::rnorm(2L * n), 2L))
    out <- data.frame(id = rep(seq_len(n), each = length(times)),
                      time = rep(times, n), x = stats::rnorm(n * length(times)))
    person_class <- class[out$id]
    out$y <- c(2, 6)[person_class] + c(0.5, -0.4)[person_class] * out$time +
      0.3 * out$x + effects[out$id, 1L] + effects[out$id, 2L] * out$time +
      stats::rnorm(nrow(out), 0, 0.7)
    out[stats::runif(nrow(out)) > 0.2, ]
  })
}

# Dense multivariate normal log density by a hand Cholesky factorization.
review_dense_log_density <- function(y, mean, covariance) {
  factor <- chol(covariance)
  z <- backsolve(factor, y - mean, transpose = TRUE)
  -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(z^2))
}

# The random-effect covariance of each class from the estimates, as the
# parametrisation defines it: one matrix per class ("varying"), one shared
# ("equal"), or one shared and scaled by a class factor ("proportional").
review_class_covariances <- function(spec, params) {
  lapply(seq_len(spec$n_classes), function(k) {
    switch(spec$random_covariance,
           varying = params$random_covariance[[k]],
           equal = params$random_covariance[[1L]],
           proportional = params$random_covariance[[1L]] * params$random_scale[k]^2)
  })
}

# Class k's regression coefficients: its own, then the common ones.
review_class_coefficients <- function(params, k) c(params$beta[, k], params$common)

# Every person's class log density from the parameters, densely.
review_dense_terms <- function(fit, params = fit$params) {
  spec <- fit$spec
  covariances <- review_class_covariances(spec, params)
  rows <- split(seq_len(spec$n), spec$group_index)
  design <- cbind(spec$x, spec$z)
  random <- spec$random_design
  vapply(seq_len(spec$n_classes), function(k) {
    beta <- review_class_coefficients(params, k)
    vapply(rows, function(i) {
      z <- random[i, , drop = FALSE]
      covariance <- z %*% covariances[[k]] %*% t(z) + diag(params$sigma2[k], length(i))
      review_dense_log_density(spec$y[i], as.vector(design[i, , drop = FALSE] %*% beta),
                               covariance)
    }, numeric(1))
  }, numeric(length(rows))) |> matrix(length(rows), spec$n_classes)
}

review_growth_fit <- function(data, n_classes = 2, ...) {
  mixture_regression(y ~ time + x, data, n_classes = n_classes, id = "id",
                     class_level = "group", common = ~ x, random = ~ 1 + time,
                     n_starts = 2, seed = 1, vcov_type = "none", ...)
}

test_that("the per-person log density and posteriors equal a dense computation", {
  skip_on_cran()
  fit <- review_growth_fit(review_growth_data())
  dense <- review_dense_terms(fit)
  independent <- "independent: dense per-person growth densities and empirical Bayes moments, base R"
  # 1e-8, the original tolerance: a Cholesky solve per person.
  expect_reference(unname(fit$expectation$log_density), unname(dense),
                   reference = independent, covers = "mixture_regression()",
                   quantity = "per-person class log densities", tolerance = 1e-8,
                   measure = "relative",
                   note = "Package side: the E step's per-person class log densities of mixture_regression(random = ); no exported accessor returns them.")
  # No membership covariates: every person's prior is the softmax of the
  # intercept logits.
  logits <- fit$params$gamma[1L, ]
  shares <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
  prior <- matrix(shares, nrow(dense), length(shares), byrow = TRUE)
  expect_reference(as.numeric(logLik(fit)), sum(log(rowSums(prior * exp(dense)))),
                   reference = independent,
                   covers = c("mixture_regression()", "logLik.latents_growth_mixture()"),
                   quantity = "log likelihood", tolerance = 1e-10, measure = "relative")
  # Empirical Bayes: E(b | y, k) = G Z' V^-1 r, Var = (G^-1 + Z'Z / sigma^2)^-1.
  spec <- fit$spec
  covariances <- review_class_covariances(spec, fit$params)
  rows <- split(seq_len(spec$n), spec$group_index)
  invisible(lapply(seq_len(spec$n_classes), function(k) {
    beta <- review_class_coefficients(fit$params, k)
    sigma2 <- fit$params$sigma2[k]
    expected <- t(vapply(rows, function(i) {
      z <- spec$random_design[i, , drop = FALSE]
      v <- z %*% covariances[[k]] %*% t(z) + diag(sigma2, length(i))
      r <- spec$y[i] - cbind(spec$x, spec$z)[i, , drop = FALSE] %*% beta
      as.vector(covariances[[k]] %*% t(z) %*% solve(v, r))
    }, numeric(2L)))
    # Scaled: the effects are deviations centred at zero, so some are near
    # zero, where a relative measure would judge rounding noise. 1e-8, the original.
    expect_reference(unname(fit$expectation$means[[k]]), unname(expected),
                     reference = independent, covers = "mixture_regression()",
                     quantity = sprintf("empirical Bayes means, class %d", k), tolerance = 1e-8,
                     measure = "scaled",
                     note = "Package side: the E step's posterior random-effect means; no exported accessor returns them.")
    person <- 3L
    z <- spec$random_design[rows[[person]], , drop = FALSE]
    expect_reference(unname(fit$expectation$covariances[[k]][, , fit$stats$pattern[person]]),
                     unname(solve(solve(covariances[[k]]) + crossprod(z) / sigma2)),
                     reference = independent, covers = "mixture_regression()",
                     quantity = sprintf("empirical Bayes covariance of person 3, class %d", k),
                     tolerance = 1e-8, measure = "relative",
                     note = "Package side: the E step's posterior random-effect covariance; no exported accessor returns it.")
  }))
})

test_that("the correlated-residual log density equals a dense computation", {
  skip_on_cran()
  data <- withr::with_seed(5L, {
    times <- c(0, 0.5, 1.5, 2, 3)
    n <- 50L
    out <- data.frame(id = rep(seq_len(n), each = 5L), time = rep(times, n))
    process <- as.vector(t(chol(0.5 * exp(-1.2 * abs(outer(times, times, "-"))))) %*%
                           matrix(stats::rnorm(5L * n), 5L))
    out$y <- 1 + 0.4 * out$time + stats::rnorm(n)[out$id] + process +
      stats::rnorm(nrow(out), 0, 0.4)
    out[stats::runif(nrow(out)) > 0.1, ]
  })
  invisible(lapply(c("ar", "brownian"), function(kind) {
    fit <- quietly(mixture_regression(  # boundary warnings are not under test
      y ~ time, data, n_classes = 1, id = "id", class_level = "group", random = ~ 1,
      n_starts = 1, seed = 1, vcov_type = "none", residual_correlation = kind,
      residual_time = "time"), c(.multilpa_expected_warnings, "latents_random_boundary"))
    spec <- fit$spec
    p <- fit$params
    rows <- split(seq_len(spec$n), spec$group_index)
    dense <- vapply(rows, function(i) {
      t <- spec$residual_times[i]
      process <- if (identical(kind, "ar")) {
        p$correlation[["sd"]]^2 * exp(-p$correlation[["rate"]] * abs(outer(t, t, "-")))
      } else p$correlation[["sd"]]^2 * outer(t, t, pmin)
      covariance <- p$random_covariance[[1L]][1L, 1L] + process +
        diag(p$sigma2[1L], length(i))
      review_dense_log_density(spec$y[i], as.vector(spec$x[i, , drop = FALSE] %*% p$beta[, 1L]),
                               covariance)
    }, numeric(1))
    expect_reference(fit$expectation$log_density[, 1L], unname(dense),
                     reference = "independent: dense log density with a correlated residual process, base R",
                     covers = "mixture_regression()",
                     quantity = sprintf("per-person log densities, residual correlation %s", kind),
                     tolerance = 1e-10, measure = "relative",
                     note = "Package side: the E step's per-person log densities of mixture_regression(residual_correlation = ); no exported accessor returns them.")
  }))
})

test_that("a random-effect covariance in small units is factorized without an absolute jitter", {
  # Regression: .growth_factor() added 1e-12 * max(1, scale) to G, an
  # absolute floor that swamps a covariance near 1e-14 and biased the density.
  covariance <- matrix(c(4, 1, 1, 2), 2L) * 1e-14
  factor <- latents:::.growth_factor(covariance)
  # The factor reproduces the matrix it was given (a unit check of the input).
  # Scaled to order one: waldo compares values below the tolerance absolutely.
  expect_equal(tcrossprod(factor) * 1e14, covariance * 1e14, tolerance = 1e-10)
  stats <- list(n = 3L, yy = 0, xy = matrix(0, 1L, 1L), zy = matrix(0, 1L, 2L),
                xx = array(0, c(1L, 1L, 1L)), zz = array(0, c(2L, 2L, 1L)),
                xz = array(0, c(1L, 2L, 1L)))
  z <- cbind(1, c(0, 1, 2))
  y <- c(1, -2, 0.5) * 1e-7
  x <- matrix(1, 3L, 1L)
  stats$yy <- sum(y^2)
  stats$xy[] <- crossprod(x, y)
  stats$zy[] <- crossprod(z, y)
  stats$xx[] <- crossprod(x)
  stats$zz[] <- crossprod(z)
  stats$xz[] <- crossprod(x, z)
  stats <- latents:::.growth_patterns(stats)
  sigma2 <- 1e-14
  terms <- latents:::.growth_class_terms(stats, 0, factor, sigma2)
  expect_reference(terms$log_density,
                   review_dense_log_density(y, rep(0, 3L),
                                            z %*% covariance %*% t(z) + diag(sigma2, 3L)),
                   reference = "independent: dense multivariate normal log density by Cholesky, base R",
                   covers = "mixture_regression()",
                   quantity = "log density with a covariance near 1e-14", tolerance = 1e-9,
                   measure = "relative",
                   note = "Package side: the internal class density (.growth_class_terms) of the growth engine, from its own factor (.growth_factor) of the covariance.")
})
