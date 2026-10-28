library(testthat)
library(latents)

log_sum <- function(x) max(x) + log(sum(exp(x - max(x))))

test_that("one-profile Gaussian fit and inference match independent MLE formulas", {
  set.seed(20261020)
  d <- data.frame(a = rnorm(150, 3, 1.5), b = rnorm(150, -2, 0.7))
  f <- lpa(d, c("a", "b"), n_profiles = 1, model = "VVV",
           n_starts = 1, seed = 7, tol = 1e-12)
  x <- as.matrix(d)
  mu <- colMeans(x)
  centered <- sweep(x, 2, mu)
  sigma <- crossprod(centered) / nrow(x)
  ll <- -nrow(x) / 2 * (2 * log(2 * pi) +
    as.numeric(determinant(sigma, logarithm = TRUE)$modulus) + 2)
  expect_equal(as.numeric(f$means), mu, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(unname(f$covariances[, , 1]), unname(sigma), tolerance = 1e-10)
  expect_equal(f$log_likelihood, ll, tolerance = 1e-10)
  info <- parameter_inference(f)
  mean_rows <- info$parameter == "mean"
  expect_equal(info$standard_error[mean_rows], sqrt(diag(sigma) / nrow(x)),
               ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("weighted Gaussian fit matches closed form weighted MLE", {
  set.seed(20261021)
  d <- data.frame(a = rnorm(180), b = rnorm(180), w = exp(rnorm(180)))
  f <- lpa(d, c("a", "b"), n_profiles = 1, model = "VVV", weights = "w",
           n_starts = 1, seed = 7, tol = 1e-12)
  x <- as.matrix(d[c("a", "b")])
  w <- d$w / sum(d$w) * nrow(x)
  mu <- colSums(x * w) / sum(w)
  r <- sweep(x, 2, mu)
  sigma <- crossprod(r, r * w) / sum(w)
  expect_equal(as.numeric(f$means), mu, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(unname(f$covariances[, , 1]), unname(sigma), tolerance = 1e-10)
  info <- parameter_inference(f)
  expect_true(all(is.finite(info$standard_error)))
  expect_identical(attr(info, "vcov_type"), "robust")
})

test_that("nested mixture evaluate-only likelihood matches exact marginalization", {
  d <- data.frame(id = rep(1:4, c(2, 3, 2, 3)),
                  y = c(-2, 1, 0, 1, 2, -1, -2, 2, 1, -1))
  start <- list(means = matrix(c(-1, 2), 2, 1),
                variances = matrix(c(0.8, 1.2), 2, 1),
                profile_probabilities = rbind(c(0.8, 0.2), c(0.3, 0.7)),
                group_probabilities = c(0.4, 0.6))
  f <- multilpa(d, "y", "id", n_profiles = 2, n_group_classes = 2,
                start = start, max_iter = 0, n_starts = 1)
  group_ll <- vapply(split(d$y, d$id), function(y) {
    group_terms <- vapply(1:2, function(g) {
      log(start$group_probabilities[g]) + sum(vapply(y, function(value) {
        log_sum(log(start$profile_probabilities[g, ]) +
          dnorm(value, start$means[, 1], sqrt(start$variances[, 1]), log = TRUE))
      }, numeric(1)))
    }, numeric(1))
    log_sum(group_terms)
  }, numeric(1))
  expect_equal(f$log_likelihood, sum(group_ll), tolerance = 1e-12)
  expect_equal(rowSums(f$subject_posteriors), rep(1, nrow(d)),
               ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(rowSums(f$group_posteriors), rep(1, 4),
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("FIML likelihood equals independent observed-margin Gaussian densities", {
  set.seed(20261022)
  d <- data.frame(a = rnorm(160), b = rnorm(160))
  d$a[seq(2, 160, 4)] <- NA_real_
  d$b[seq(3, 160, 4)] <- NA_real_
  f <- lpa(d, c("a", "b"), n_profiles = 1, model = "VVV", missing = "fiml",
           n_starts = 1, seed = 7, tol = 1e-12)
  sigma <- f$covariances[, , 1]
  densities <- vapply(seq_len(nrow(d)), function(i) {
    observed <- which(!is.na(d[i, ]))
    r <- as.numeric(d[i, observed]) - f$means[1, observed]
    s <- sigma[observed, observed, drop = FALSE]
    -0.5 * (length(observed) * log(2 * pi) +
      as.numeric(determinant(s, logarithm = TRUE)$modulus) +
      as.numeric(crossprod(r, solve(s, r))))
  }, numeric(1))
  expect_equal(f$log_likelihood, sum(densities), tolerance = 1e-11)
  expect_true(f$converged)
})
