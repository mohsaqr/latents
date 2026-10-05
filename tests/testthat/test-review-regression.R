# Regression tests from past audits; CI runs them on every platform.
skip_on_cran()

review_regression_rows <- function() {
  set.seed(31)
  data <- data.frame(x = rnorm(160), g = rep(1:40, each = 4))
  class <- rep(rep(1:2, each = 20), each = 4)
  data$y <- c(-3, 3)[class] + .5 * data$x + rnorm(160, sd = .4)
  data$w <- c(1, 3)[class]
  data
}

test_that("weighted regression class sizes use weighted posteriors", {
  data <- review_regression_rows()
  lapply(c("observation", "group"), function(level) {
    fit <- mixture_regression(y ~ x, data, 2, weights = "w", n_starts = 1,
                               seed = 1, class_level = level,
                               id = if (level == "group") "g" else NULL)
    table <- get_results(fit, "classes")
    tau <- if (level == "group") fit$expectation$group_tau else fit$expectation$tau
    weights <- if (level == "group") fit$spec$sampling_weights else fit$spec$row_weights
    expected <- colSums(tau * weights)
    expect_equal(table$count, unname(expected))
    expect_equal(table$share, unname(expected / sum(expected)))
    expect_equal(get_results(fit, "fit")$smallest_share, min(table$share))
    expect_equal(unname(rowSums(tau)), rep(1, nrow(tau)))
  })
})

test_that("two-level regression group shares and composition use sampling weights", {
  data <- review_regression_rows()
  fit <- quietly(mixture_regression(y ~ x, data, 2, id = "g", n_group_classes = 2,
                                     weights = "w", n_starts = 2, seed = 1,
                                     membership = ~ x, vcov_type = "none"),
                 classes = c(.multilpa_expected_warnings, "latents_separation"))
  # This table is defined without Wald inference when row covariates occur.
  table <- .mixture_group_class_table(fit, .95, NULL)
  weighted <- fit$expectation$rho * fit$spec$sampling_weights
  expect_equal(table$group_share, rep(colSums(weighted) / sum(weighted), each = 2),
               ignore_attr = TRUE)
  lapply(1:2, function(h) {
    prior <- exp(fit$expectation$log_prior_by_group_class[[h]])
    weight <- weighted[fit$spec$group_index, h]
    expect_equal(table$probability[table$group_class == paste0("group_class_", h)],
                 unname(colSums(prior * weight) / sum(weight)))
  })
})

test_that("regression counts and iteration controls reject fractional or overflowing values", {
  data <- review_regression_rows()
  lapply(c(.5, 2^31, Inf, NA_real_), function(value) {
    expect_error(mixture_regression(y ~ x, data, 1, n_starts = value), "n_starts")
    expect_error(mixture_regression(y ~ x, data, 1, max_iter = value), "max_iter")
    expect_error(enumerate_regressions(y ~ x, data, 1, bootstrap = value), "bootstrap")
    expect_error(enumerate_regressions(y ~ x, data, 1, bootstrap_starts = value),
                 "bootstrap_starts")
  })
  expect_error(mixture_regression(y ~ x, data, 1.000000001), "n_classes")
  expect_error(mixture_regression(y ~ x, data, 2^31), "n_classes")
  expect_error(mixture_regression(y ~ x, data, 1, seed = Inf), class = "latents_bad_argument")
  expect_error(enumerate_regressions(y ~ x, data, 1, seed = .1),
               class = "latents_bad_argument")
  fit <- quietly(mixture_regression(y ~ x, data, 1, n_starts = 0, max_iter = 0,
                                    vcov_type = "none"))
  expect_identical(fit$settings$n_starts, 0L)
})

review_regression_likelihood <- function(spec, parameters) {
  eta <- spec$x %*% parameters$beta + spec$offset
  eta <- eta + drop(spec$z %*% parameters$common)
  density <- vapply(1:2, function(k) {
    switch(spec$family,
      gaussian = dnorm(spec$y, eta[, k], sqrt(parameters$sigma2[k])),
      poisson = dpois(spec$y, exp(eta[, k])),
      binomial = dbinom(spec$y, spec$trials, plogis(eta[, k])))
  }, numeric(spec$n))
  softmax <- function(eta) exp(eta) / rowSums(exp(eta))
  if (spec$nesting == "observation") {
    likelihood <- rowSums(softmax(spec$w %*% parameters$gamma) * density)
    return(sum(spec$row_weights * log(likelihood)))
  }
  likelihood <- vapply(seq_len(spec$n_groups), function(j) {
    rows <- which(spec$group_index == j)
    if (spec$nesting == "group") {
      prior <- softmax(spec$w[j, , drop = FALSE] %*% parameters$gamma)
      return(sum(prior * apply(density[rows, , drop = FALSE], 2L, prod)))
    }
    group_prior <- softmax(spec$v[j, , drop = FALSE] %*% parameters$delta)
    sum(vapply(1:2, function(h) {
      group_design <- matrix(0, length(rows), 2L)
      group_design[, h] <- 1
      class_prior <- softmax(cbind(group_design, spec$w[rows, , drop = FALSE]) %*%
                               parameters$class_logits)
      group_prior[h] * prod(rowSums(class_prior * density[rows, , drop = FALSE]))
    }, numeric(1)))
  }, numeric(1))
  sum(spec$sampling_weights * log(likelihood))
}

test_that("weighted regression likelihoods and scores hold for all families and nestings", {
  set.seed(291)
  data <- data.frame(g = rep(1:6, each = 4), x = rnorm(24),
                     u = rep(rnorm(6), each = 4), w = rep(1:3, each = 8),
                     y = rnorm(24), k = rpois(24, 3), s = rbinom(24, 5, .4))
  data$f <- 5 - data$s
  grid <- expand.grid(family = c("gaussian", "poisson", "binomial"),
                      nesting = c("observation", "group", "two-level"),
                      stringsAsFactors = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    family <- grid$family[i]
    nesting <- grid$nesting[i]
    formula <- switch(family, gaussian = y ~ x + u, poisson = k ~ x + u,
                       binomial = cbind(s, f) ~ x + u)
    spec <- .mixture_spec(formula, data, 2L, family,
                          if (nesting == "observation") NULL else "g",
                          if (nesting == "group") "group" else "observation",
                          if (nesting == "two-level") 2L else 1L,
                          common = ~ u, membership = ~ u,
                          group_membership = if (nesting == "two-level") ~ u else ~1,
                          variance = "varying", min_variance = 1e-6, missing = "error")
    spec <- .mixture_weight_spec(spec, data, "w")
    template <- .mixture_null_params(spec)
    theta <- .mixture_pack(spec, template)
    theta[] <- sin(seq_along(theta)) / 4
    parameters <- .mixture_unpack(spec, theta, template)
    objective <- function(point) {
      review_regression_likelihood(spec, .mixture_unpack(spec, point, template))
    }
    actual <- .mixture_expectation(spec, parameters)
    expect_equal(actual$log_likelihood, objective(theta), tolerance = 1e-12,
                 info = paste(family, nesting))
    gradient <- vapply(seq_along(theta), function(j) {
      up <- down <- theta
      up[j] <- up[j] + 1e-5
      down[j] <- down[j] - 1e-5
      (objective(up) - objective(down)) / 2e-5
    }, numeric(1))
    score <- colSums(.mixture_unit_scores(spec, parameters, actual))
    expect_equal(unname(score), gradient, tolerance = 1e-6,
                 info = paste(family, nesting))
  })
})
