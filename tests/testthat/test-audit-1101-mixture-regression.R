skip_on_cran()

# Audit 2026-11-01, slice `mixture-regression`: regression tests for the fixes
# and independent checks of the regression-mixture engine (M-step, ordinal
# likelihood, random intercept, predictions, standard errors).

audit_regression_data <- function(n_groups = 60L, per_group = 5L, seed = 4L) {
  set.seed(seed)
  n <- n_groups * per_group
  data <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                     id = rep(seq_len(n_groups), each = per_group))
  row_class <- stats::rbinom(n, 1L, stats::plogis(0.3 + 0.8 * data$z)) + 1L
  group_class <- rep(stats::rbinom(n_groups, 1L, 0.5), each = per_group) + 1L
  data$y <- ifelse(row_class == 1L, 1 + 2 * data$x, -1 - data$x) +
    stats::rnorm(n, sd = 0.6)
  data$yg <- ifelse(group_class == 1L, 1 + 2 * data$x, -1 - data$x) +
    stats::rnorm(n, sd = 0.6)
  latent <- ifelse(group_class == 1L, 1 + 1.5 * data$x, -1 - data$x) +
    stats::rlogis(n)
  data$o <- as.integer(cut(latent, c(-Inf, -1, 0.5, 2, Inf)))
  data$e <- stats::runif(n)
  latent_common <- ifelse(group_class == 1L, 1, -1) + 0.7 * data$x + data$e +
    stats::rlogis(n)
  data$oc <- as.integer(cut(latent_common, c(-Inf, -0.5, 1, 2.5, Inf)))
  data
}

# Observed-data log likelihood of a fit as a function of its packed vector.
audit_log_likelihood <- function(fit) {
  function(theta) {
    params <- latents:::.mixture_unpack(fit$spec, theta, fit$params)
    latents:::.mixture_expectation(fit$spec, params)$log_likelihood
  }
}

audit_hessian_errors <- function(fit) {
  theta <- latents:::.mixture_pack(fit$spec, fit$params)
  hessian <- numDeriv::hessian(audit_log_likelihood(fit), theta)
  sqrt(diag(solve(-hessian)))
}

test_that("an ordinal mixture with every slope shared (or none) is fitted and inferred", {
  # Fix: .mixture_pack() named a zero-length class block with one stray name,
  # so `common` holding every predictor, or `o ~ 1`, stopped with an
  # unclassed base error once inference ran.
  data <- audit_regression_data()
  fit <- mixture_regression(oc ~ x + offset(e), data, 2, family = "ordinal",
                            id = "id", class_level = "group", common = "x",
                            n_starts = 2, seed = 1, tol = 1e-12)
  table <- as.data.frame(fit)
  expect_identical(table$term,
                   c(rep(c("threshold:1|2", "threshold:2|3", "threshold:3|4"), 2), "x"))
  expect_true(all(is.finite(table$std_error)))
  # The likelihood written from the model: per class, cumulative logits with
  # the shared slope and the offset; groups share one class.
  params <- fit$params
  row_log <- vapply(1:2, function(k) {
    eta <- params$common[["x"]] * data$x + data$e
    cumulative <- cbind(0, stats::plogis(outer(eta, params$thresholds[, k],
                                               function(e, t) t - e)), 1)
    log(cumulative[cbind(seq_len(nrow(data)), data$oc + 1L)] -
          cumulative[cbind(seq_len(nrow(data)), data$oc)])
  }, numeric(nrow(data)))
  prior <- exp(params$gamma[1, ]) / sum(exp(params$gamma[1, ]))
  by_group <- rowsum(row_log, data$id)
  expect_equal(fit$log_likelihood, sum(log(exp(by_group) %*% prior)),
               tolerance = 1e-10)
  skip_if_not_installed("numDeriv")
  expect_equal(unname(sqrt(diag(fit$inference$vcov))), audit_hessian_errors(fit),
               tolerance = 1e-5)

  no_predictor <- mixture_regression(o ~ 1, data, 2, family = "ordinal", id = "id",
                                     class_level = "group", n_starts = 2, seed = 1)
  expect_identical(nrow(as.data.frame(no_predictor)), 6L)
  expect_false(any(startsWith(names(coef(no_predictor)), "coefficient.")))
})

test_that("bootstrap replicates do not repeat boundary warnings", {
  set.seed(2)
  data <- data.frame(x = stats::rnorm(120))
  data$count <- stats::rpois(120, exp(0.5 + 0.4 * data$x))
  boundary_count <- function(bootstrap) {
    count <- 0L
    withCallingHandlers(
      enumerate_regressions(count ~ x, data, n_classes = 1:2,
                            family = "negative_binomial", n_starts = 2,
                            seed = 1, bootstrap = bootstrap),
      latents_boundary = function(w) {
        count <<- count + 1L
        invokeRestart("muffleWarning")
      },
      latents_failed_replicates = function(w) invokeRestart("muffleWarning"))
    count
  }
  observed_only <- boundary_count(0L)
  expect_gt(observed_only, 0L)
  expect_identical(boundary_count(3L), observed_only)
})

test_that("the Gaussian M-step is weighted least squares per class", {
  data <- audit_regression_data()
  data$f <- factor(rep(c("a", "b", "c"), length.out = nrow(data)))
  spec <- latents:::.mixture_spec(y ~ x + f, data, 2L, "gaussian", NULL, "observation",
                                  1L, NULL, ~1, ~1, "varying", 1e-6, "error")
  set.seed(3)
  first <- stats::runif(nrow(data))
  tau <- cbind(first, 1 - first)
  updated <- latents:::.mixture_update_regression(spec, latents:::.mixture_null_params(spec), tau)
  updated <- latents:::.mixture_update_variance(spec, updated, tau)
  fits <- lapply(1:2, function(k) stats::lm(y ~ x + f, data, weights = tau[, k]))
  # One step from zero carries the 1e-10 ridge (about 1e-9 here) ...
  expect_equal(unname(updated$beta), unname(vapply(fits, stats::coef, numeric(4))),
               tolerance = 1e-8)
  expect_equal(unname(updated$sigma2), vapply(1:2, function(k) {
    sum(tau[, k] * stats::residuals(fits[[k]])^2) / sum(tau[, k])
  }, numeric(1)), tolerance = 1e-10)
  # ... and the weighted least-squares solution is its exact fixed point.
  at_solution <- updated
  at_solution$beta[] <- vapply(fits, stats::coef, numeric(4))
  again <- latents:::.mixture_update_regression(spec, at_solution, tau)
  expect_equal(again$beta, at_solution$beta, tolerance = 1e-10)
})

test_that("a one-class ordinal mixture is MASS::polr", {
  skip_if_not_installed("MASS")
  data <- audit_regression_data()
  data$ordered <- factor(data$o, ordered = TRUE)
  fit <- mixture_regression(ordered ~ x + z, data, 1, family = "ordinal",
                            n_starts = 1, seed = 1, tol = 1e-12)
  reference <- MASS::polr(ordered ~ x + z, data, Hess = TRUE)
  expect_equal(fit$log_likelihood, as.numeric(stats::logLik(reference)), tolerance = 1e-8)
  table <- as.data.frame(fit)
  expect_equal(table$estimate, unname(c(reference$zeta, stats::coef(reference))),
               tolerance = 1e-4)
  expect_equal(table$std_error, unname(sqrt(diag(stats::vcov(reference)))[c(3:5, 1:2)]),
               tolerance = 1e-3)
})

test_that("a one-class random intercept is nlme::lme by maximum likelihood", {
  skip_if_not_installed("nlme")
  set.seed(3)
  data <- data.frame(id = rep(1:40, each = 5), t = rep(0:4, 40))
  data$y <- 2 + 0.5 * data$t + rep(stats::rnorm(40, sd = 0.8), each = 5) +
    stats::rnorm(200, sd = 0.6)
  fit <- mixture_regression(y ~ t, data, 1, id = "id", class_level = "group",
                            random = "intercept", n_starts = 1, seed = 1)
  reference <- nlme::lme(y ~ t, random = ~ 1 | id, data = data, method = "ML")
  expect_equal(fit$log_likelihood, as.numeric(stats::logLik(reference)), tolerance = 1e-8)
})

test_that("predictions are the design times the coefficients, mixed by the priors", {
  data <- audit_regression_data()
  data$f <- factor(rep(c("a", "b"), length.out = nrow(data)))
  fit <- mixture_regression(y ~ x + f, data, 2, membership = ~z, common = "f",
                            n_starts = 3, seed = 1)
  new <- data.frame(x = c(-1, 0.5, 2), z = c(0, 1, -1), f = factor(c("a", "b", "b")))
  linear <- cbind(1, new$x) %*% fit$params$beta + (new$f == "b") * fit$params$common[["fb"]]
  logits <- cbind(1, new$z) %*% fit$params$gamma
  prior <- exp(logits) / rowSums(exp(logits))
  by_class <- predict(fit, new, type = "class_response")
  expect_equal(by_class$fitted, as.vector(linear), tolerance = 1e-12)
  expect_equal(by_class$prior, as.vector(prior), tolerance = 1e-12)
  expect_equal(predict(fit, new)$fitted, rowSums(prior * linear), tolerance = 1e-12)

  set.seed(5)
  counts <- data.frame(x = stats::rnorm(200), exposure = stats::runif(200, 1, 3))
  counts$y <- stats::rpois(200, counts$exposure *
                             exp(ifelse(seq_len(200) %% 2 == 0, 0.3 + 0.8 * counts$x, 1.8 - 0.4 * counts$x)))
  poisson <- mixture_regression(y ~ x + offset(log(exposure)), counts, 2,
                                family = "poisson", n_starts = 2, seed = 1)
  new_counts <- counts[1:4, ]
  expect_equal(predict(poisson, new_counts, type = "class_response")$fitted,
               as.vector(exp(cbind(1, new_counts$x) %*% poisson$params$beta +
                               log(new_counts$exposure))), tolerance = 1e-12)
})

test_that("standard errors are the inverse numerical Hessian in every variant", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("MASS")
  data <- audit_regression_data()
  set.seed(9)
  data$count <- stats::rpois(nrow(data), exp(ifelse(data$y > 0, 0.2 + 0.7 * data$x, 1.5 - 0.3 * data$x)))
  data$nb <- MASS::rnegbin(nrow(data), exp(ifelse(data$y > 0, 0.2 + 0.7 * data$x, 1.5 - 0.3 * data$x)), theta = 3)
  data$trials <- rep(4:9, length.out = nrow(data))
  data$successes <- stats::rbinom(nrow(data), data$trials, stats::plogis(ifelse(data$y > 0, -1 + data$x, 1 - data$x)))
  fits <- list(
    gaussian = mixture_regression(y ~ x, data, 2, membership = ~z, n_starts = 3, seed = 1, tol = 1e-12),
    group = mixture_regression(yg ~ x, data, 2, id = "id", class_level = "group", n_starts = 3, seed = 1, tol = 1e-12),
    two_level = mixture_regression(y ~ x, data, 2, id = "id", n_group_classes = 2, n_starts = 3, seed = 1, tol = 1e-12),
    poisson = mixture_regression(count ~ x, data, 2, family = "poisson", n_starts = 3, seed = 1, tol = 1e-12),
    negative_binomial = mixture_regression(nb ~ x, data, 2, family = "negative_binomial", n_starts = 3, seed = 1, tol = 1e-12),
    binomial = mixture_regression(cbind(successes, trials - successes) ~ x, data, 2, family = "binomial", n_starts = 3, seed = 1, tol = 1e-12),
    ordinal = mixture_regression(o ~ x, data, 2, family = "ordinal", id = "id", class_level = "group", n_starts = 3, seed = 1, tol = 1e-12))
  invisible(lapply(names(fits), function(name) {
    fit <- fits[[name]]
    expect_equal(unname(sqrt(diag(fit$inference$vcov))), audit_hessian_errors(fit),
                 tolerance = 1e-4, label = name)
  }))
})
