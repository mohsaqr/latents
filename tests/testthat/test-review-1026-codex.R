# Adversarial interaction regressions, 2026-10-26 (round 3).

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

test_that("ordinary category labels cannot become partial joint observations", {
  data <- data.frame(
    a = factor(rep(c("a:x", "b:x", "a::"), c(40, 30, 30))),
    b = rep(c("u", "v"), 50))
  fit <- lca(data, names(data), 1, n_starts = 1, missing = "fiml", seed = 1)
  # With one class each item is a multinomial with its empirical proportions.
  expected <- 40 * log(.4) + 60 * log(.3) + 100 * log(.5)
  expect_equal(fit$categorical_levels$a, levels(data$a))
  expect_equal(unname(fit$response_probabilities$a[1, ]), c(.3, .4, .3))
  expect_equal(as.numeric(logLik(fit)), expected, tolerance = 1e-12)
  expect_null(attr(fit$categorical_data, "latents_partial"))
  expect_equal(sum(predict(fit, data, type = "density")$log_density), expected,
               tolerance = 1e-12)
})

test_that("prediction accepts partial markers only for declared joint items", {
  data <- data.frame(a = factor(rep(c("a:x", "b:x"), each = 20)))
  fit <- lca(data, "a", 1, n_starts = 1, missing = "fiml", seed = 1)
  expect_error(predict(fit, data.frame(a = factor("a::")), type = "density"),
               class = "latents_bad_data")
  # Genuine paired items still sum out their missing item and survive rebuilding.
  paired <- data.frame(a = c(rep("a", 4), rep("b", 4), "a", NA),
                       b = c("x", "x", "x", "y", "x", "y", "y", "y", NA, "x"))
  joint <- lca(paired, c("a", "b"), 1, local_dependence = "a:b",
               missing = "fiml", n_starts = 1, seed = 1, tol = 1e-12)
  probability <- drop(joint$response_probabilities[["a:b"]])
  half <- c(sum(probability[c("a:x", "a:y")]),
             sum(probability[c("a:x", "b:x")]))
  expect_equal(tail(predict(joint, paired, type = "density")$log_density, 2),
               log(half), tolerance = 1e-12)
  expect_equal(predict(joint, type = "density"),
               predict(joint, paired, type = "density"), tolerance = 1e-12)
  rebuilt <- get_results(joint, "data")
  again <- lca(rebuilt, "a:b", 1, missing = "fiml", n_starts = 1, seed = 1,
                tol = 1e-12)
  expect_equal(as.numeric(logLik(again)), as.numeric(logLik(joint)), tolerance = 1e-12)
  rebuilt$id <- seq_len(nrow(rebuilt))
  sampled <- .multilpa_resample_groups(rebuilt, as.list(rebuilt$id), "id", rebuilt$id)
  expect_identical(attr(sampled[["a:b"]], "latents_pair"), c("a", "b"))
  resampled_fit <- lca(sampled, "a:b", 1, missing = "fiml", n_starts = 1,
                       seed = 1, tol = 1e-12)
  expect_equal(as.numeric(logLik(resampled_fit)), as.numeric(logLik(joint)),
               tolerance = 1e-12)
  # Saved fits from before pair provenance was retained still identify their
  # joint indicators in local_dependence.
  legacy <- joint
  attr(legacy$categorical_values[["a:b"]], "latents_pair") <- NULL
  expect_equal(vcov(legacy), vcov(joint), tolerance = 1e-12)
})

test_that("compare_models retains released model labels that now name controls", {
  data <- data.frame(y = c(-2, -1, 0, 1, 2, 3))
  fit <- mixture_regression(y ~ 1, data, 1, n_starts = 0, vcov_type = "none")
  # All three names were arbitrary model labels in 0.9.21's dots-only API.
  comparison <- compare_models(n_profiles = fit, other = fit,
                                cutoff = fit, n_group_classes = fit)
  expect_identical(comparison$model,
                    c("n_profiles", "other", "cutoff", "n_group_classes"))
  expect_equal(comparison$log_likelihood, rep(as.numeric(logLik(fit)), 4))
  expect_equal(comparison$bic_weight, rep(.25, 4))
  # Numeric controls retain their current grid-selection behavior.
  grid <- enumerate_lpa(data, "y", 1:2, n_starts = 1, seed = 1)
  selected <- compare_models(grid, n_profiles = 1:2, cutoff = .7)
  expect_setequal(selected$n_profiles, 1:2)
  expect_error(compare_models(other = fit, n_profiles = 2),
               class = "latents_bad_argument")
})

test_that("zero-weight rows do not increase observed indicator counts under FIML", {
  .latents_local_seed(312)
  data <- data.frame(
    group = rep(seq_len(12), each = 4), x = stats::rnorm(48),
    y = stats::rnorm(48), item = factor(sample(1:3, 48, TRUE)), weight = 1)
  data$x[c(2, 9)] <- NA_real_
  data$y[c(8, 33, 39)] <- NA_real_
  data$item[c(1, 17)] <- NA
  absent <- data.frame(
    group = rep(13:15, each = 4), x = c(NA_real_, seq_len(11)),
    y = seq_len(12), item = factor(rep(1:3, 4), levels = levels(data$item)),
    weight = 0)
  padded <- rbind(data, absent)
  fit <- function(frame, drop_empty) suppressWarnings(multilpa(
    frame, c("x", "y", "item"), "group", n_profiles = 1,
    n_group_classes = 1, categorical = "item", missing = "fiml",
    weights = "weight", n_starts = 1, seed = 3, drop_empty = drop_empty))
  expected <- vapply(padded[padded$weight > 0, c("x", "y", "item")],
                     function(column) as.numeric(sum(!is.na(column))), numeric(1))
  invisible(lapply(c(TRUE, FALSE), function(drop_empty) {
    reference <- fit(data, drop_empty)
    actual <- fit(padded, drop_empty)
    expect_equal(actual$n_observed_by_indicator, expected)
    expect_equal(actual$n_observed_by_indicator, reference$n_observed_by_indicator)
    expect_equal(actual$log_likelihood, reference$log_likelihood, tolerance = 1e-12)
  }))
})

test_that("positive sampling weights stay present when their ratio underflows", {
  data <- data.frame(x = c(100, -1, 1), weight = c(1e-300, 1e300, 1e300))
  start <- list(means = matrix(0), variances = matrix(1),
                profile_probabilities = matrix(1), group_probabilities = 1)
  actual <- suppressMessages(multilpa(
    data, "x", id = NULL, n_profiles = 1, weights = "weight", start = start,
    n_starts = 1, max_iter = 0))
  # At double precision the tiny first contribution is zero, but all three
  # raw-positive units set the normalization. Each large weight is 3 / 2.
  expected <- 1.5 * sum(stats::dnorm(c(-1, 1), log = TRUE))
  expect_true(all(actual$sampling_weights > 0))
  expect_equal(sum(actual$sampling_weights), 3)
  expect_equal(nobs(actual), 3)
  expect_equal(actual$log_likelihood, expected, tolerance = 1e-12)
  expect_equal(actual$bic, -2 * expected + 2 * log(3), tolerance = 1e-12)
  # Finite weights whose sum overflows must still normalize on their ratios.
  overflow <- data.frame(weight = rep(1e308, 3))
  expect_equal(.latents_sampling_weights(overflow, "weight", 1:3, 3), rep(1, 3))
})

# Cross-family propagation of the sampling-weight and complete-entropy rules.
codex_criteria_entropy <- function(posterior, weights) {
  cells <- ifelse(posterior > 0, -posterior * log(posterior), 0)
  sum(weights * rowSums(cells))
}

codex_criteria_avepp <- function(posterior, weights) {
  assigned <- max.col(posterior, ties.method = "first")
  vapply(seq_len(ncol(posterior)), function(k) {
    keep <- assigned == k & weights > 0
    if (any(keep)) sum(weights[keep] * posterior[keep, k]) / sum(weights[keep]) else NA_real_
  }, numeric(1))
}

codex_criteria_regression_data <- function() {
  .latents_local_seed(42)
  n <- 180L
  data <- data.frame(x = stats::rnorm(n), w = sample(c(1, 5), n, TRUE))
  latent <- stats::rbinom(n, 1, 0.5)
  data$y <- ifelse(latent == 1, 1 + 0.4 * data$x, -1 - 0.4 * data$x) + stats::rnorm(n)
  data
}

codex_criteria_group_data <- function() {
  .latents_local_seed(91)
  data <- data.frame(id = rep(seq_len(65), each = 4), time = rep(0:3, 65))
  data$w <- rep(sample(c(1, 5), 65, TRUE), each = 4)
  data$y <- rep(stats::rnorm(65, rep(c(-0.8, 0.8), length.out = 65), 0.6), each = 4) +
    0.3 * data$time + stats::rnorm(nrow(data), sd = 0.9)
  data
}

test_that("regression classification estimates follow the sampling weights", {
  fit <- mixture_regression(y ~ x, codex_criteria_regression_data(), 2, weights = "w",
                            n_starts = 3, seed = 2, vcov_type = "none")
  posterior <- fit$expectation$tau
  weights <- fit$spec$row_weights
  entropy <- codex_criteria_entropy(posterior, weights)
  table <- get_results(fit, "fit")
  expect_equal(table$icl, table$bic + 2 * entropy, tolerance = 1e-12)
  expect_equal(table$entropy, 1 - entropy / (sum(weights) * log(2)), tolerance = 1e-12)
  expect_equal(get_results(fit, "classes")$mean_posterior,
               codex_criteria_avepp(posterior, weights), tolerance = 1e-12)
  classification <- get_results(fit, "classification")
  diagonal <- classification$assigned == classification$class
  expect_equal(classification$mean_posterior[diagonal],
               codex_criteria_avepp(posterior, weights), tolerance = 1e-12)
})

test_that("growth classification estimates follow the sampling weights", {
  fit <- suppressWarnings(mixture_regression(
    y ~ time, codex_criteria_group_data(), 2, id = "id", class_level = "group",
    random = ~ 1, weights = "w", n_starts = 2, seed = 3, vcov_type = "none"))
  posterior <- fit$expectation$posterior
  weights <- fit$spec$sampling_weights
  entropy <- codex_criteria_entropy(posterior, weights)
  table <- get_results(fit, "fit")
  expect_equal(table$icl, table$bic + 2 * entropy, tolerance = 1e-12)
  expect_equal(table$entropy, 1 - entropy / (sum(weights) * log(2)), tolerance = 1e-12)
  expect_equal(get_results(fit, "classes")$average_posterior,
               codex_criteria_avepp(posterior, weights), tolerance = 1e-12)
})

test_that("additive classification estimates follow the sampling weights", {
  fit <- suppressWarnings(multilpa(
    codex_criteria_group_data(), "y", "id", n_group_classes = 2, family = "additive",
    weights = "w", n_starts = 2, seed = 3))
  posterior <- fit$group_posteriors
  weights <- fit$sampling_weights
  entropy <- codex_criteria_entropy(posterior, weights)
  expect_equal(get_results(fit, "fit")$entropy,
               1 - entropy / (sum(weights) * log(2)), tolerance = 1e-12)
  expect_equal(get_results(fit, "group_classes")$mean_posterior,
               codex_criteria_avepp(posterior, weights), tolerance = 1e-12)
  classification <- get_results(fit, "classification")
  diagonal <- classification$assigned == classification$group_class
  expect_equal(classification$mean_posterior[diagonal],
               codex_criteria_avepp(posterior, weights), tolerance = 1e-12)
})

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

codex_partial_covariate_data <- function() {
  .latents_local_seed(91)
  n <- 300L
  latent <- sample(1:2, n, TRUE)
  data <- data.frame(a = sample(c("a", "b"), n, TRUE),
                     b = sample(c("x", "y"), n, TRUE),
                     y = stats::rnorm(n, 3 * (latent - 1)), z = stats::rnorm(n))
  data$b[seq(1, 100, 3)] <- NA
  data
}

test_that("paired FIML covariate inference retains each half-observed item", {
  data <- codex_partial_covariate_data()
  fit <- lpa(data, c("a", "b", "y"), 2, categorical = c("a", "b"),
              local_dependence = "a:b", profile_covariates = "z", missing = "fiml",
              n_starts = 2, seed = 1)
  rebuilt <- .multilpa_cov_stored_data(fit)
  half <- is.na(data$b)
  expect_identical(as.character(rebuilt[["a:b"]])[half], paste0(data$a[half], "::"))
  expect_identical(attr(rebuilt[["a:b"]], "latents_pair"), c("a", "b"))
  stored <- vcov(fit, scale = "unconstrained")
  supplied <- vcov(fit, data = data, scale = "unconstrained")
  expect_equal(stored, supplied, tolerance = 1e-12)
  expect_true(all(is.finite(stored)))
  legacy <- fit
  attr(legacy$categorical_values[["a:b"]], "latents_pair") <- NULL
  expect_equal(vcov(legacy, scale = "unconstrained"), stored, tolerance = 1e-12)
  changed <- data
  changed$a[1] <- if (data$a[1] == "a") "b" else "a"
  expect_error(vcov(fit, data = changed), class = "latents_bad_inference_data")
})

# Round-three cross-family review: absent units must not set Gaussian floors
# or growth coordinates; Brownian diffusion SD has units y/sqrt(time).

codex_growth_weight_data <- function() {
  .latents_local_seed(42L)
  data <- data.frame(id = rep(seq_len(40L), each = 4L),
                     time = rep(0:3, 40L), weight = 1)
  data$y <- 2 + 0.4 * data$time + stats::rnorm(40L)[data$id] +
    stats::rnorm(nrow(data), sd = 0.5)
  data
}

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

test_that("a constant Gaussian outcome on carrying rows is refused by class", {
  data <- codex_growth_weight_data()
  data$weight[data$id > 3L] <- 0
  data$y[data$weight > 0] <- 0.1
  expect_error(mixture_regression(
    y ~ time, data, n_classes = 1L, id = "id", weights = "weight",
    n_starts = 0L, seed = 1L), class = "latents_bad_data")
})

test_that("relative regression variance floors remain relative below machine epsilon", {
  data <- codex_growth_weight_data()
  fit <- function(d, random = NULL, weights = NULL) mixture_regression(
    y ~ time, d, n_classes = 1L, id = "id", class_level = "group",
    random = random, weights = weights, n_starts = 0L, seed = 1L,
    vcov_type = "none")
  ordinary <- fit(data)
  tiny <- fit(transform(data, y = y * 1e-12))
  # Density under y -> c*y changes by the exact Jacobian -N*log(c).
  expect_equal(tiny$log_likelihood,
               ordinary$log_likelihood - nrow(data) * log(1e-12), tolerance = 1e-9)
  expect_equal(tiny$spec$min_variance / 1e-24, ordinary$spec$min_variance)
  expect_equal(tiny$params$sigma2 / 1e-24, ordinary$params$sigma2, tolerance = 1e-9)
  weighted <- fit(transform(data, y = y * 1e-12), weights = "weight")
  expect_equal(weighted$spec$min_variance, tiny$spec$min_variance)
  growth <- fit(data, ~ 1)
  tiny_growth <- fit(transform(data, y = y * 1e-12), ~ 1)
  expect_equal(tiny_growth$log_likelihood,
               growth$log_likelihood - nrow(data) * log(1e-12), tolerance = 1e-8)
  expect_equal(tiny_growth$params$sigma2 / 1e-24, growth$params$sigma2, tolerance = 1e-8)
})

test_that("singleton persons cannot identify a random intercept and residual", {
  .latents_local_seed(3L)
  data <- data.frame(id = seq_len(30L), x = stats::rnorm(30L),
                     y = stats::rnorm(30L), weight = 1)
  # The one-row marginal covariance is G_00 + sigma^2, unchanged along this
  # nonzero parameter direction: there are three identified parameters (OLS).
  expect_equal(0.2 + 0.8, 0.3 + 0.7)
  expect_equal(attr(stats::logLik(stats::lm(y ~ x, data)), "df"), 3)
  expect_error(mixture_regression(
    y ~ x, data, n_classes = 1L, id = "id", class_level = "group",
    random = ~ 1, n_starts = 0L, seed = 1L), class = "latents_not_identified")
  absent <- data[1:4, ]
  absent$id <- 31L
  absent$weight <- 0
  expect_error(mixture_regression(
    y ~ x, rbind(data, absent), n_classes = 1L, id = "id", class_level = "group",
    weights = "weight", random = ~ 1, n_starts = 0L, seed = 1L),
    class = "latents_not_identified")
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
