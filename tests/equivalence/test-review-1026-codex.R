# Comparisons with independent references moved from tests/testthat/test-review-1026-codex.R,
# so that they run in CI's equivalence job rather than under R CMD check.

codex_growth_weight_data <- function() {
  .latents_local_seed(42L)
  data <- data.frame(id = rep(seq_len(40L), each = 4L),
                     time = rep(0:3, 40L), weight = 1)
  data$y <- 2 + 0.4 * data$time + stats::rnorm(40L)[data$id] +
    stats::rnorm(nrow(data), sd = 0.5)
  data
}

test_that("two-level regression ICL charges the brute-force joint entropy", {
  .latents_local_seed(60)
  data <- data.frame(id = rep(seq_len(40), each = 3), x = stats::rnorm(120))
  group_class <- rep(stats::rbinom(40, 1, 0.5), each = 3)
  row_class <- stats::rbinom(120, 1, ifelse(group_class == 1, 0.75, 0.25))
  data$y <- ifelse(row_class == 1, 1 + 0.2 * data$x, -1 - 0.2 * data$x) +
    stats::rnorm(120, sd = 1.1)
  fit <- suppressWarnings(mixture_regression(
    y ~ x, data, 2, id = "id", n_group_classes = 2, n_starts = 3, seed = 1,
    vcov_type = "none"))
  # Enumerate all 2^3 row-class paths and both group classes; no package
  # entropy helper or stored posterior enters the independent reference.
  paths <- as.matrix(expand.grid(rep(list(1:2), 3)))
  prior <- exp(fit$params$class_logits)
  prior <- prior / rowSums(prior)
  group_prior <- exp(fit$params$delta[1, ])
  group_prior <- group_prior / sum(group_prior)
  density <- vapply(seq_len(2), function(k) {
    stats::dnorm(data$y, fit$spec$x %*% fit$params$beta[, k],
                 sqrt(fit$params$sigma2[k]))
  }, numeric(nrow(data)))
  entropy <- sum(vapply(seq_len(40), function(j) {
    members <- which(data$id == j)
    joint <- unlist(lapply(seq_len(2), function(h) {
      apply(paths, 1, function(path) {
        group_prior[h] * prod(prior[h, path] * density[cbind(members, path)])
      })
    }))
    posterior <- joint / sum(joint)
    -sum(posterior[posterior > 0] * log(posterior[posterior > 0]))
  }, numeric(1)))
  table <- get_results(fit, "fit")
  expect_equal(table$icl, table$bic + 2 * entropy, tolerance = 1e-10)
})

test_that("zero-weight outcomes cannot change Gaussian residual floors", {
  data <- codex_growth_weight_data()
  absent <- data.frame(id = rep(41:44, each = 4L), time = rep(0:3, 4L),
                       weight = 0, y = rep(c(-1e5, 1e5), each = 8L))
  fit <- function(d, random = NULL) mixture_regression(
    y ~ time, d, n_classes = 1L, id = "id", class_level = "group",
    random = random, weights = "weight", n_starts = 0L, seed = 1L,
    vcov_type = "none")
  reference <- stats::lm(y ~ time, data)
  ordinary <- fit(data)
  padded <- fit(rbind(data, absent))
  # Gaussian ML, independently from OLS: beta and sigma^2 = SSE / N.
  expect_equal(unname(ordinary$params$beta[, 1L]), unname(stats::coef(reference)),
               tolerance = 1e-10)
  expect_equal(unname(ordinary$params$sigma2), mean(stats::residuals(reference)^2),
               tolerance = 1e-10)
  expect_equal(padded$log_likelihood, ordinary$log_likelihood, tolerance = 1e-9)
  expect_equal(padded$spec$min_variance, ordinary$spec$min_variance)
  growth <- fit(data, ~ 1)
  growth_padded <- fit(rbind(data, absent), ~ 1)
  expect_equal(growth_padded$log_likelihood, growth$log_likelihood, tolerance = 1e-7)
  expect_equal(growth_padded$params$sigma2, growth$params$sigma2, tolerance = 1e-7)
  # The transform is estimated from the carrying rows, by its hand definition.
  expect_equal(growth_padded$standard$center, mean(data$y))
  expect_equal(growth_padded$standard$scale, stats::sd(data$y))
})

test_that("Brownian growth residuals preserve their model when time units change", {
  .latents_local_seed(12L)
  data <- local({
    n <- 50L
    out <- data.frame(id = rep(seq_len(n), each = 5L), time = rep(0:4, n))
    process <- unlist(lapply(seq_len(n), function(i) {
      c(0, cumsum(stats::rnorm(4L, sd = 0.8)))
    }))
    out$y <- 2 + 0.3 * out$time + rep(stats::rnorm(n, sd = 0.6), each = 5L) +
      process + stats::rnorm(nrow(out), sd = 0.3)
    out
  })
  fit <- function(d) mixture_regression(
    y ~ time, d, n_classes = 1L, id = "id", class_level = "group",
    random = ~ 1, residual_correlation = "brownian", n_starts = 0L,
    seed = 1L, vcov_type = "none")
  base <- fit(data)
  moved <- fit(transform(data, time = time * 1e8))
  # Independent dense Gaussian likelihood, W_ij = diffusion_sd^2 min(t_i,t_j).
  direct <- sum(vapply(split(data, data$id), function(person) {
    covariance <- base$params$random_covariance[[1L]][1L, 1L] +
      base$params$correlation[["sd"]]^2 * outer(person$time, person$time, pmin) +
      diag(base$params$sigma2, nrow(person))
    residual <- person$y - base$params$beta[1L, 1L] -
      base$params$beta[2L, 1L] * person$time
    -0.5 * (nrow(person) * log(2 * pi) +
      as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
      sum(residual * solve(covariance, residual)))
  }, numeric(1)))
  expect_equal(base$log_likelihood, direct, tolerance = 1e-9)
  expect_equal(moved$log_likelihood, direct, tolerance = 1e-8)
  expect_equal(moved$params$correlation[["sd"]] * 1e4,
               base$params$correlation[["sd"]], tolerance = 1e-7)
  expect_equal(moved$params$sigma2, base$params$sigma2, tolerance = 1e-7)
  expect_length(moved$residual_boundary, 0L)
  expect_equal(predict(moved, transform(data, time = time * 1e8),
                       type = "conditional")$fitted,
               predict(base, data, type = "conditional")$fitted, tolerance = 1e-7)
  data$weight <- 1
  absent <- data.frame(id = rep(51:54, each = 5L), time = rep(0:4, 4L) * 1e12,
                       weight = 0, y = rep(c(-1e5, 1e5), each = 10L))
  padded <- mixture_regression(
    y ~ time, rbind(data, absent), n_classes = 1L, id = "id", class_level = "group",
    random = ~ 1, residual_correlation = "brownian", weights = "weight",
    n_starts = 0L, seed = 1L, vcov_type = "none")
  expect_equal(padded$log_likelihood, direct, tolerance = 1e-8)
  expect_equal(padded$params$correlation[["sd"]], base$params$correlation[["sd"]],
               tolerance = 1e-7)
})
