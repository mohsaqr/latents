# Review 2026-10-26, growth mixture models: calibration of the marginal
# likelihood y_i | k ~ N(X_i beta_k, Z_i G_k Z_i' + sigma_k^2 I) (Verbeke and
# Lesaffre 1996; Muthen and Shedden 1999) against a dense computation and
# nlme::lme(method = "ML"), the correlated-residual covariances of
# lcmm::hlme(cor = AR / BM) (Proust-Lima, Philipps and Liquet 2017), and
# invariants (class relabelling, time origin, outcome scale).

# More than 1.5 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

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

test_that("analytic scores equal central differences of the log likelihood", {
  skip_on_cran()
  fit <- review_growth_fit(review_growth_data())
  spec <- fit$spec
  theta <- latents:::.growth_pack(spec, fit$params) +
    withr::with_seed(2L, stats::rnorm(fit$n_parameters, 0, 0.05))
  log_likelihood <- function(v) {
    latents:::.growth_expectation(spec, fit$stats,
                                  latents:::.growth_unpack(spec, v, fit$params))$log_likelihood
  }
  numeric_gradient <- vapply(seq_along(theta), function(j) {
    h <- 1e-5 * (1 + abs(theta[[j]]))
    (log_likelihood(replace(theta, j, theta[j] + h)) -
       log_likelihood(replace(theta, j, theta[j] - h))) / (2 * h)
  }, numeric(1))
  analytic <- colSums(latents:::.growth_scores(
    spec, fit$stats, latents:::.growth_unpack(spec, theta, fit$params)))
  expect_equal(unname(analytic), numeric_gradient, tolerance = 1e-6)
})

test_that("AR and Brownian residual covariances follow their formulas and are positive definite", {
  spec_ar <- list(residual_correlation = "ar")
  spec_bm <- list(residual_correlation = "brownian")
  times <- c(0, 0.5, 1.5, 2, 3.25)
  ar <- latents:::.growth_process_covariance(spec_ar, c(rate = 1.3, sd = 0.7), times)
  bm <- latents:::.growth_process_covariance(spec_bm, c(sd = 0.7), times)
  expect_equal(ar[2L, 4L], 0.49 * exp(-1.3 * 1.5))
  expect_equal(diag(ar), rep(0.49, 5L))
  expect_equal(bm[3L, 5L], 0.49 * 1.5)
  expect_equal(bm[5L, 2L], 0.49 * 0.5)
  expect_true(isSymmetric(ar) && isSymmetric(bm))
  # exp(-rho |t - s|) is positive definite on distinct times; min(t, s) is
  # positive semi-definite on nonnegative times (singular at t = 0).
  expect_gt(min(eigen(ar, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_gte(min(eigen(bm, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
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

test_that("relabelling the classes leaves the likelihood unchanged and swaps posteriors", {
  skip_on_cran()
  fit <- review_growth_fit(review_growth_data())
  swapped <- latents:::.growth_permute_params(fit$spec, fit$params, 2:1)
  expectation <- latents:::.growth_expectation(fit$spec, fit$stats, swapped)
  expect_equal(expectation$log_likelihood, fit$log_likelihood, tolerance = 1e-12)
  expect_equal(expectation$posterior[, 2:1], fit$expectation$posterior, tolerance = 1e-12)
})

test_that("shifting the time origin leaves an intercept-and-slope model's likelihood unchanged", {
  skip_on_cran()
  data <- review_growth_data()
  fit <- review_growth_fit(data)
  shifted_data <- transform(data, time = time - 2.5)
  shifted <- review_growth_fit(shifted_data)
  # Exactly: b_new = T b with T the change of origin, so beta and G map by T.
  transform_matrix <- matrix(c(1, 0, 2.5, 1), 2L)
  params <- fit$params
  params$beta[c("(Intercept)", "time"), ] <- transform_matrix %*%
    params$beta[c("(Intercept)", "time"), ]
  params$random_covariance <- lapply(params$random_covariance, function(g) {
    transform_matrix %*% g %*% t(transform_matrix)
  })
  mapped <- latents:::.growth_expectation(shifted$spec, shifted$stats, params)
  expect_equal(mapped$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  # And the refitted maximum agrees.
  expect_equal(shifted$log_likelihood, fit$log_likelihood, tolerance = 1e-6)
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

test_that("an outcome in small units fits like the same outcome in ordinary units", {
  skip_on_cran()
  data <- review_growth_data(n = 80L, seed = 11L, classes = 1L)
  ordinary <- review_growth_fit(data, n_classes = 1)
  small <- data
  small$y <- small$y * 1e-5
  # Before the relative Cholesky jitter, EM stopped with "likelihood
  # decreased". Scale equivariance: logLik(c y) = logLik(y) - N log(c). The
  # quasi-Newton finish is not scale-free (reported), so only 0.05 here.
  fit <- review_growth_fit(small, n_classes = 1)
  expect_lt(abs(fit$log_likelihood - (ordinary$log_likelihood - nrow(data) * log(1e-5))),
            0.05)
})

test_that("a growth fit is equivariant in the outcome's units (G2)", {
  withr::local_seed(11)
  n <- 80L
  covariance <- matrix(c(1, 0.3, 0.3, 0.25), 2L)
  effects <- t(t(chol(covariance)) %*% matrix(stats::rnorm(2L * n), 2L))
  data <- data.frame(id = rep(seq_len(n), each = 5L), time = rep(0:4, n))
  data$y <- 2 + 0.5 * data$time + effects[data$id, 1L] + effects[data$id, 2L] * data$time +
    stats::rnorm(nrow(data), 0, 0.7)
  data <- data[stats::runif(nrow(data)) > 0.15, ]
  fit_at <- function(scale) {
    mixture_regression(y ~ time, transform(data, y = y * scale), n_classes = 1, id = "id",
                       class_level = "group", random = ~ 1 + time, n_starts = 1, seed = 1)
  }
  table_of <- function(fit) subset(as.data.frame(parameter_inference(fit)), block != "share")
  base <- fit_at(1)
  reference <- table_of(base)
  power <- c(coefficient = 1, residual_sd = 1, random_variance = 2, random_sd = 1,
             random_correlation = 0)[reference$block] |> unname()
  invisible(lapply(c(1e-5, 1e5), function(scale) {
    fit <- fit_at(scale)
    expect_true(fit$converged)
    expect_equal(fit$log_likelihood, base$log_likelihood - nrow(data) * log(scale),
                 tolerance = 1e-10, info = format(scale))
    table <- table_of(fit)
    expected_se <- reference$standard_error * scale^power
    # Within a thousandth of a standard error; before the fix 0.06 to 0.085.
    expect_lt(max(abs(table$estimate - reference$estimate * scale^power) / expected_se), 1e-3)
    # Before the fix one standard error was 6.2 times too large at 1e-5.
    expect_equal(table$standard_error, expected_se, tolerance = 1e-4, info = format(scale))
  }))
})

test_that("positive_rows is the exact support count, not an underflow count (GD-1)", {
  withr::local_seed(909)
  data <- data.frame(y1 = c(stats::rnorm(140L, 0, 1), stats::rnorm(100L, 3, 1),
                            stats::rnorm(60L, 1.5, 0.01)))
  data$unit <- seq_len(nrow(data))
  fit <- quietly(multilpa(data, vars = "y1", id = "unit", n_profiles = 3L,
                          n_group_classes = 1L, n_starts = 10L, seed = 1L, tol = 1e-12))
  # Before the fix the narrow component reported 87: the rows whose posterior
  # had not underflowed.
  expect_identical(get_results(fit, "gaussian_diagnostics")$positive_rows, rep(300L, 3L))
  data$y1[c(5L, 9L)] <- NA
  refit <- quietly(multilpa(data, vars = "y1", id = "unit", n_profiles = 3L,
                            n_group_classes = 1L, n_starts = 10L, seed = 1L,
                            missing = "fiml"))
  expect_identical(get_results(refit, "gaussian_diagnostics")$positive_rows, rep(298L, 3L))
})

test_that("a trajectory inference names why a weight column cannot be read", {
  data <- subset(growth_scores, student <= 40)
  data$w <- 1
  fit <- quietly(mixture_regression(score ~ wave, data, n_classes = 1, id = "student",
                                    class_level = "group", random = ~ 1, weights = "w",
                                    n_starts = 0, seed = 1, vcov_type = "none"))
  broken <- data
  broken$w <- NULL
  # three_step() and r3step() refuse weighted trajectory fits before this
  # check, so the helper is exercised directly: a missing weight column used
  # to be swallowed into a generic refusal.
  failure <- expect_error(.trajectory_check_three_step(fit, broken),
                          class = "latents_bad_inference_data")
  expect_match(conditionMessage(failure), "weights")
})
