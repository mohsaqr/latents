test_that("empty groups do not evade robust validity or alter CR1", {
  set.seed(792)
  make_fit <- function(data, ...) quietly(multilpa(data, "y", "g", 1, 1,
    missing = "fiml", n_starts = 1, seed = 3, tol = 1e-12, ...))
  invisible(lapply(1:2, function(n_groups) {
    data <- data.frame(g = rep(seq_len(n_groups), each = 50),
                       y = stats::rnorm(50 * n_groups))
    padded <- rbind(data, data.frame(g = n_groups + seq_len(20), y = NA_real_))
    invisible(lapply(c("robust", "opg"), function(type) {
      expect_error(parameter_inference(make_fit(data), vcov_type = type),
                   class = "latents_too_few_groups")
      expect_error(parameter_inference(make_fit(padded), vcov_type = type),
                   class = "latents_too_few_groups")
      expect_error(parameter_inference(make_fit(padded, drop_empty = FALSE), vcov_type = type),
                   class = "latents_too_few_groups")
    }))
  }))
  data <- data.frame(g = rep(1:6, each = 30), y = stats::rnorm(180))
  padded <- rbind(data, data.frame(g = 7:26, y = NA_real_))
  invisible(lapply(c("robust", "opg"), function(type) {
    a <- parameter_inference(make_fit(data), vcov_type = type)
    b <- parameter_inference(make_fit(padded), vcov_type = type)
    expect_equal(attr(b, "covariance"), attr(a, "covariance"), tolerance = 1e-8)
    expect_equal(attr(b, "group_scores")[1:6, ], attr(a, "group_scores"),
                 tolerance = 1e-8)
    expect_lt(max(abs(attr(b, "group_scores")[-(1:6), ])), 1e-12)
  }))
  # Legacy counting intentionally uses every empty group in CR1.
  legacy <- parameter_inference(make_fit(padded, drop_empty = FALSE), vcov_type = "robust")
  current <- parameter_inference(make_fit(padded), vcov_type = "robust")
  expect_equal(attr(legacy, "covariance"),
               attr(current, "covariance") * (26 / 25) / (6 / 5), tolerance = 1e-8)
  # Normalizing weights changes the pseudo likelihood, but cancels in the sandwich.
  data$weight <- rep(c(1, 2, 1, 3, 2, 1), each = 30)
  padded <- rbind(data, data.frame(g = 7:26, y = NA_real_, weight = 2))
  zero <- rbind(data, data.frame(g = 7:26, y = stats::rnorm(20), weight = 0))
  weighted <- function(d) parameter_inference(make_fit(d, weights = "weight"))
  expect_equal(attr(weighted(padded), "covariance"), attr(weighted(data), "covariance"),
               tolerance = 1e-8)
  expect_equal(attr(weighted(zero), "covariance"), attr(weighted(data), "covariance"),
               tolerance = 1e-8)
  # One-dimensional scores have full numerical rank even when the lone
  # informative unit's score is only roundoff. Legacy CR1 cannot authorize
  # an otherwise invalid covariance, including with fixed measurement.
  one <- subset(data, g == 1)
  start <- starting_values(make_fit(one))
  one_padded <- rbind(one, data.frame(g = 7:26, y = NA_real_, weight = 1))
  fixed <- make_fit(one_padded, start = start, fixed = "means", drop_empty = FALSE)
  expect_error(parameter_inference(fixed, vcov_type = "robust"),
               class = "latents_too_few_groups")
})

test_that("direct-effect inference transforms the full covariance across predictor origins", {
  set.seed(802)
  data <- data.frame(g = rep(1:20, each = 20), z = stats::rnorm(400),
                     w = stats::rnorm(400))
  data$y <- 0.8 + 0.7 * data$z - 0.4 * data$w + stats::rnorm(400)
  fit <- function(d, ...) quietly(multilpa(d, "y", "g", 1, 1,
    direct_effects = "y ~ z + w", n_starts = 1, seed = 3, tol = 1e-12, ...))
  original <- fit(data)
  shifted_data <- transform(data, z = z + 2026, w = w - 1000)
  shifted <- fit(shifted_data)
  # Coordinates are intercept, log variance, z slope, w slope.
  transform <- diag(4)
  transform[1, 3:4] <- c(-2026, 1000)
  inverse_transform <- diag(4)
  inverse_transform[1, 3:4] <- -transform[1, 3:4]
  invisible(lapply(c("observed", "robust", "opg"), function(type) {
    a <- parameter_inference(original, vcov_type = type)
    b <- parameter_inference(shifted, vcov_type = type)
    expected <- transform %*% attr(a, "covariance_unconstrained") %*% t(transform)
    expect_equal(unname(attr(b, "covariance_unconstrained")), unname(expected),
                 tolerance = 2e-7)
    expect_equal(b$standard_error[b$parameter == "direct_effect"],
                 a$standard_error[a$parameter == "direct_effect"], tolerance = 2e-7)
    expect_equal(b$p_value[b$parameter == "direct_effect"],
                 a$p_value[a$parameter == "direct_effect"], tolerance = 2e-7)
    if (type != "observed") {
      expect_equal(unname(attr(b, "group_scores")),
                   unname(attr(a, "group_scores") %*% inverse_transform), tolerance = 2e-6)
    }
    # The cached observed Hessian must work for a second covariance type too.
    expect_equal(unname(attr(b, "hessian")),
                 unname(t(inverse_transform) %*% attr(a, "hessian") %*%
                          inverse_transform), tolerance = 2e-7)
  }))
  reference <- stats::lm(y ~ z + w, shifted_data)
  mle_variance <- sum(stats::residuals(reference)^2) / nrow(data)
  reference_covariance <- summary(reference)$cov.unscaled * mle_variance
  inferred <- parameter_inference(shifted)
  expect_equal(unname(attr(inferred, "covariance")[c(1, 3, 4), c(1, 3, 4)]),
               unname(reference_covariance), tolerance = 2e-7)
  expect_equal(inferred$estimate[c(1, 3, 4)], unname(stats::coef(reference)), tolerance = 2e-7)
  # Independently form the normal-model sandwich in the original predictor
  # coordinates, then apply the known intercept transformation.
  reference <- stats::lm(y ~ z + w, data)
  residual <- stats::residuals(reference)
  variance <- mean(residual^2)
  design <- stats::model.matrix(reference)
  unit_scores <- rowsum(cbind(design * (residual / variance),
                              (residual^2 / variance - 1) / 2), data$g)
  inverse <- matrix(0, 4, 4)
  inverse[1:3, 1:3] <- summary(reference)$cov.unscaled * variance
  inverse[4, 4] <- 2 / nrow(data)
  reference_covariance <- inverse %*% crossprod(unit_scores) %*% inverse * 20 / 19
  order <- c(1, 4, 2, 3)
  reference_covariance <- transform %*% reference_covariance[order, order] %*% t(transform)
  expect_equal(unname(attr(parameter_inference(shifted, vcov_type = "robust"),
                            "covariance_unconstrained")),
               unname(reference_covariance), tolerance = 2e-7)
})

test_that("direct-effect translation stays regular with missing indicators and sampling weights", {
  set.seed(215)
  n <- 400
  data <- data.frame(g = rep(1:40, each = 10), z = stats::rnorm(n),
                     w = stats::rnorm(n), weight = rep(c(1, 2), each = 200))
  class <- sample(c(-3, 3), n, TRUE)
  data$a <- class + 0.6 * data$z + stats::rnorm(n)
  data$b <- class - 0.4 * data$z + 0.3 * data$w + stats::rnorm(n)
  data$a[seq(5, n, 11)] <- NA_real_
  data$b[seq(8, n, 13)] <- NA_real_
  fit <- function(d) quietly(multilpa(d, c("a", "b"), "g", 2, 1,
    direct_effects = c("a ~ z", "b ~ z + w"), weights = "weight", missing = "fiml",
    n_starts = 1, seed = 1, tol = 1e-12))
  a <- fit(data)
  b <- fit(transform(data, z = z + 2026, w = w - 1000))
  expect_equal(a$log_likelihood, b$log_likelihood, tolerance = 1e-8)
  ia <- parameter_inference(a)
  ib <- parameter_inference(b)
  effects <- ia$parameter == "direct_effect"
  expect_equal(ib$estimate[effects], ia$estimate[effects], tolerance = 2e-6)
  expect_equal(ib$standard_error[effects], ia$standard_error[effects], tolerance = 2e-6)
  expect_equal(ib$p_value[effects], ia$p_value[effects], tolerance = 2e-6)
  expect_gt(attr(ib, "condition_ratio"), 1e-5)
})

test_that("irrelevant direct covariates cannot destabilize fitting or inference", {
  set.seed(803)
  data <- data.frame(g = rep(1:20, each = 20), z = stats::rnorm(400),
                     w = stats::rnorm(400), weight = 1)
  data$y <- 0.8 + 0.7 * data$z - 0.4 * data$w + stats::rnorm(400)
  fit <- function(d, ...) quietly(multilpa(d, "y", "g", 1, 1,
    direct_effects = "y ~ z + w", missing = "fiml", n_starts = 1, seed = 1,
    tol = 1e-12, ...))
  base <- fit(data)
  padding <- data.frame(g = 21:30, z = 1e5, w = -1e5, weight = 1, y = NA_real_)
  padded <- fit(rbind(data, padding))
  expect_true(padded$converged)
  expect_equal(padded$log_likelihood, base$log_likelihood, tolerance = 1e-9)
  expect_equal(parameter_inference(padded, vcov_type = "robust")$standard_error,
               parameter_inference(base, vcov_type = "robust")$standard_error,
               tolerance = 2e-7)
  padding$y <- stats::rnorm(10)
  padding$weight <- 0
  weighted <- fit(data, weights = "weight")
  zero_padded <- fit(rbind(data, padding), weights = "weight")
  expect_true(zero_padded$converged)
  expect_equal(zero_padded$direct_effects, weighted$direct_effects, tolerance = 2e-7)
  expect_equal(parameter_inference(zero_padded)$standard_error,
               parameter_inference(weighted)$standard_error, tolerance = 2e-7)
})

test_that("direct-effect charts use target-specific observations and preserve public covariates", {
  set.seed(804)
  data <- data.frame(g = rep(1:40, each = 10), z = stats::rnorm(400),
                     w = stats::rnorm(400))
  data$z[201:400] <- data$z[201:400] + 1e5
  data$a <- 0.8 + 0.7 * data$z + stats::rnorm(400)
  data$a[201:400] <- NA_real_
  data$b <- -0.3 - 0.4 * data$z + 0.2 * data$w + stats::rnorm(400)
  fit <- quietly(multilpa(data, c("a", "b"), "g", 1, 1,
    direct_effects = c("a ~ z", "b ~ z + w"), missing = "fiml", n_starts = 1,
    seed = 1, tol = 1e-12))
  expect_true(fit$converged)
  inference <- parameter_inference(fit)
  references <- list(a = stats::lm(a ~ z, data), b = stats::lm(b ~ z + w, data))
  invisible(lapply(names(references), function(indicator) {
    reference <- references[[indicator]]
    variance <- mean(stats::residuals(reference)^2)
    rows <- which((inference$parameter == "mean" & inference$term == indicator) |
      (inference$parameter == "direct_effect" & startsWith(inference$term, paste0(indicator, ":"))))
    expect_equal(inference$estimate[rows], unname(stats::coef(reference)), tolerance = 2e-6)
    expect_equal(unname(attr(inference, "covariance")[rows, rows]),
                 unname(summary(reference)$cov.unscaled * variance), tolerance = 2e-6)
  }))
  expect_identical(colnames(fit$extra_data$direct$covariates), c("z", "w"))
  expect_identical(fit$extra_data$direct$effects$covariate_column, c(1L, 1L, 2L))
  expect_identical(.multilpa_refit_arguments(fit)$direct_effects, c("a ~ z", "b ~ z + w"))
  expect_equal(predict(fit, newdata = data, type = "posterior"),
               predict(fit, type = "posterior"), tolerance = 1e-12)
})

test_that("covariate and both transition engines preserve robust inference under empty padding", {
  set.seed(34)
  data <- data.frame(g = rep(1:40, each = 5), time = rep(1:5, 40),
                     z = stats::rnorm(200))
  data$y <- stats::rnorm(200, sample(c(-3, 3), 200, TRUE))
  padding <- data.frame(g = rep(41:60, each = 5), time = rep(1:5, 20),
                        z = stats::rnorm(100), y = NA_real_)
  fitters <- list(
    covariate = function(d, ...) quietly(multilpa(d, "y", "g", 2, 1,
      profile_covariates = "z", missing = "fiml", n_starts = 1, seed = 1,
      tol = 1e-12, ...)),
    homogeneous = function(d, ...) quietly(lta(d, "y", "g", 2, time = "time",
      missing = "fiml", n_starts = 1, seed = 1, tol = 1e-12, ...)),
    general = function(d, ...) quietly(lta(d, "y", "g", 2, time = "time",
      transitions = "occasion", missing = "fiml", n_starts = 1, seed = 1,
      tol = 1e-12, ...)))
  covariance <- function(fit, type) {
    if (inherits(fit, "multilpa_lta")) vcov(fit, type = type) else
      vcov(fit, vcov_type = type)
  }
  invisible(lapply(fitters, function(fitter) {
    base <- fitter(data)
    padded <- fitter(rbind(data, padding))
    expect_true(base$converged && padded$converged)
    expect_equal(base$log_likelihood, padded$log_likelihood, tolerance = 1e-9)
    invisible(lapply(c("robust", "opg"), function(type) {
      expect_equal(covariance(padded, type), covariance(base, type), tolerance = 3e-6)
    }))
    legacy <- fitter(rbind(data, padding), drop_empty = FALSE)
    expect_equal(covariance(legacy, "robust"), covariance(padded, "robust") *
                   (60 / 59) / (40 / 39), tolerance = 3e-6)
    data$weight <- rep(c(1, 2), each = 100)
    padding$weight <- 2
    weighted_base <- fitter(data, weights = "weight")
    weighted_pad <- fitter(rbind(data, padding), weights = "weight")
    expect_equal(covariance(weighted_pad, "robust"),
                 covariance(weighted_base, "robust"), tolerance = 3e-6)
  }))
})
