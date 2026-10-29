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

# Every person's class log density from the parameters, densely.
review_dense_terms <- function(fit, params = fit$params) {
  spec <- fit$spec
  covariances <- latents:::.growth_covariances(spec, params)
  rows <- split(seq_len(spec$n), spec$group_index)
  design <- cbind(spec$x, spec$z)
  random <- spec$random_design
  vapply(seq_len(spec$n_classes), function(k) {
    beta <- latents:::.growth_coefficients(params, k)
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
  expect_equal(fit$expectation$log_density, dense, tolerance = 1e-8)
  prior <- latents:::.growth_class_prior(fit$spec, fit$params)
  expect_equal(fit$log_likelihood, sum(log(rowSums(prior * exp(dense)))),
               tolerance = 1e-10)
  # Empirical Bayes: E(b | y, k) = G Z' V^-1 r, Var = (G^-1 + Z'Z / sigma^2)^-1.
  spec <- fit$spec
  covariances <- latents:::.growth_covariances(spec, fit$params)
  rows <- split(seq_len(spec$n), spec$group_index)
  invisible(lapply(seq_len(spec$n_classes), function(k) {
    beta <- latents:::.growth_coefficients(fit$params, k)
    sigma2 <- fit$params$sigma2[k]
    expected <- t(vapply(rows, function(i) {
      z <- spec$random_design[i, , drop = FALSE]
      v <- z %*% covariances[[k]] %*% t(z) + diag(sigma2, length(i))
      r <- spec$y[i] - cbind(spec$x, spec$z)[i, , drop = FALSE] %*% beta
      as.vector(covariances[[k]] %*% t(z) %*% solve(v, r))
    }, numeric(2L)))
    expect_equal(unname(fit$expectation$means[[k]]), unname(expected), tolerance = 1e-8)
    person <- 3L
    z <- spec$random_design[rows[[person]], , drop = FALSE]
    expect_equal(fit$expectation$covariances[[k]][, , fit$stats$pattern[person]],
                 unname(solve(solve(covariances[[k]]) + crossprod(z) / sigma2)),
                 tolerance = 1e-8)
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
    fit <- suppressWarnings(mixture_regression(  # boundary warnings are not under test
      y ~ time, data, n_classes = 1, id = "id", class_level = "group", random = ~ 1,
      n_starts = 1, seed = 1, vcov_type = "none", residual_correlation = kind,
      residual_time = "time"))
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
    expect_equal(fit$expectation$log_density[, 1L], unname(dense), tolerance = 1e-10)
  }))
})

test_that("a random-effect covariance in small units is factorized without an absolute jitter", {
  # Regression: .growth_factor() added 1e-12 * max(1, scale) to G, an
  # absolute floor that swamps a covariance near 1e-14 and biased the density.
  covariance <- matrix(c(4, 1, 1, 2), 2L) * 1e-14
  factor <- latents:::.growth_factor(covariance)
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
  expect_equal(terms$log_density,
               review_dense_log_density(y, rep(0, 3L),
                                        z %*% covariance %*% t(z) + diag(sigma2, 3L)),
               tolerance = 1e-9)
})
