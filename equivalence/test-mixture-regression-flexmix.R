# mixture_regression() against flexmix (Leisch 2004; Gruen and Leisch 2008).
#
# Two checks per model. (1) The likelihood: latents' log likelihood evaluated
# at flexmix's own estimates equals flexmix's logLik(), which shows the two
# define the same model independently of how either optimizes. (2) The
# optimum: latents reaches at least flexmix's likelihood, and for the
# binomial and Poisson families, whose M-steps both packages solve exactly,
# the same estimates.
#
# flexmix's Gaussian M-step divides the residual sum of squares by
# (n - rank) rather than n, so its Gaussian estimates are not the maximum
# likelihood ones and latents' likelihood is slightly higher; that case is
# compared on (1) and on the inequality only.

skip_if_not_installed("flexmix")
flexmix_control <- list(tolerance = 1e-13, iter.max = 10000L)

best_flexmix <- function(...) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()),
          add = TRUE)
  set.seed(11)
  flexmix::stepFlexmix(..., nrep = 10, verbose = FALSE,
                       control = flexmix_control)
}

# Latents' parameter list from a flexmix fit, classes in flexmix's order.
params_from_flexmix <- function(fit, flexmix_fit) {
  spec <- fit$spec
  estimates <- flexmix::parameters(flexmix_fit)
  coefficient_rows <- startsWith(rownames(estimates), "coef.")
  params <- fit$params
  params$beta[] <- estimates[coefficient_rows, , drop = FALSE]
  if (identical(spec$family, "gaussian")) {
    params$sigma2 <- estimates["sigma", ]^2
  }
  # Without a concomitant model flexmix reports the class probabilities
  # themselves; with one, the logit coefficients with the first class as the
  # reference, which is latents' own layout.
  concomitant <- as.matrix(flexmix::parameters(flexmix_fit,
                                               which = "concomitant"))
  params$gamma <- if (identical(rownames(concomitant), "prior")) {
    matrix(log(concomitant) - log(concomitant[1L]), nrow = 1L)
  } else concomitant
  params
}

latents_at <- function(fit, params) {
  latents:::.mixture_expectation(fit$spec, params)$log_likelihood
}

data_for_comparison <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  set.seed(2)
  groups <- 80L
  size <- 6L
  id <- rep(seq_len(groups), each = size)
  group_class <- 1L + (stats::runif(groups) < 0.5)
  class <- group_class[id]
  x <- stats::rnorm(groups * size)
  w <- stats::rnorm(groups)[id]
  out <- data.frame(
    id = id, x = x, w = w, r = stats::rnorm(groups * size),
    yb = stats::rbinom(groups * size, 1L,
                       stats::plogis(ifelse(class == 1L, -1 + 1.5 * x, 1 - x))),
    yp = stats::rpois(groups * size,
                      exp(0.5 + ifelse(class == 1L, 0.8, -0.5) * x)),
    yg = ifelse(class == 1L, 2 + x, -1 - 2 * x) +
      stats::rnorm(groups * size))
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  out
})

compare <- function(ours, theirs, estimates_too) {
  params <- params_from_flexmix(ours, theirs)
  expect_equal(latents_at(ours, params), theirs@logLik,
               tolerance = 1e-9)
  expect_gte(ours$log_likelihood,
             theirs@logLik - 1e-7)
  if (estimates_too) {
    theirs_beta <- params$beta
    ours_beta <- ours$params$beta
    # Match classes by their coefficients, not by label.
    order_theirs <- order(theirs_beta[1L, ])
    order_ours <- order(ours_beta[1L, ])
    expect_equal(unname(ours_beta[, order_ours]),
                 unname(theirs_beta[, order_theirs]), tolerance = 1e-5)
  }
}

test_that("Gaussian rows", {
  ours <- mixture_regression(yg ~ x, data_for_comparison, 2, n_starts = 5, seed = 1,
                 tol = 1e-13, vcov_type = "none")
  theirs <- best_flexmix(yg ~ x, data = data_for_comparison, k = 2)
  compare(ours, theirs, estimates_too = FALSE)
})

test_that("Gaussian rows with a concomitant model", {
  ours <- mixture_regression(yg ~ x, data_for_comparison, 2, membership = ~ r,
                 n_starts = 5, seed = 1, tol = 1e-13, vcov_type = "none")
  theirs <- best_flexmix(yg ~ x, data = data_for_comparison, k = 2,
                         concomitant = flexmix::FLXPmultinom(~ r))
  compare(ours, theirs, estimates_too = FALSE)
})

test_that("Gaussian, one class per group", {
  ours <- mixture_regression(yg ~ x, data_for_comparison, 2, id = "id",
                 class_level = "group", n_starts = 5, seed = 1, tol = 1e-13,
                 vcov_type = "none")
  theirs <- best_flexmix(yg ~ x | id, data = data_for_comparison, k = 2)
  compare(ours, theirs, estimates_too = FALSE)
})

test_that("binomial, one class per group", {
  ours <- mixture_regression(yb ~ x, data_for_comparison, 2, family = "binomial",
                 id = "id", class_level = "group", n_starts = 5, seed = 1,
                 tol = 1e-13, vcov_type = "none")
  theirs <- best_flexmix(cbind(yb, 1 - yb) ~ x | id,
                         data = data_for_comparison, k = 2,
                         model = flexmix::FLXMRglm(family = "binomial"))
  compare(ours, theirs, estimates_too = TRUE)
})

test_that("binomial, one class per group, with a concomitant model", {
  ours <- mixture_regression(yb ~ x, data_for_comparison, 2, family = "binomial",
                 id = "id", class_level = "group", membership = ~ w,
                 n_starts = 5, seed = 1, tol = 1e-13, vcov_type = "none")
  theirs <- best_flexmix(cbind(yb, 1 - yb) ~ x | id,
                         data = data_for_comparison, k = 2,
                         model = flexmix::FLXMRglm(family = "binomial"),
                         concomitant = flexmix::FLXPmultinom(~ w))
  compare(ours, theirs, estimates_too = TRUE)
})

test_that("Poisson rows and Poisson per group", {
  ours <- mixture_regression(yp ~ x, data_for_comparison, 2, family = "poisson",
                 n_starts = 5, seed = 1, tol = 1e-13, vcov_type = "none")
  theirs <- best_flexmix(yp ~ x, data = data_for_comparison, k = 2,
                         model = flexmix::FLXMRglm(family = "poisson"))
  compare(ours, theirs, estimates_too = TRUE)
  ours <- mixture_regression(yp ~ x, data_for_comparison, 2, family = "poisson",
                 id = "id", class_level = "group", n_starts = 5, seed = 1,
                 tol = 1e-13, vcov_type = "none")
  theirs <- best_flexmix(yp ~ x | id, data = data_for_comparison, k = 2,
                         model = flexmix::FLXMRglm(family = "poisson"))
  compare(ours, theirs, estimates_too = TRUE)
})

test_that("Poisson with a coefficient shared across classes", {
  ours <- mixture_regression(yp ~ x + r, data_for_comparison, 2, family = "poisson",
                 common = ~ r, n_starts = 5, seed = 1, tol = 1e-13,
                 vcov_type = "none")
  theirs <- best_flexmix(yp ~ x, data = data_for_comparison, k = 2,
                         model = flexmix::FLXMRglmfix(family = "poisson",
                                                      fixed = ~ r))
  expect_equal(ours$log_likelihood, theirs@logLik,
               tolerance = 1e-7)
  shared <- flexmix::parameters(theirs)["coef.r", 1L]
  expect_equal(unname(ours$params$common), unname(shared), tolerance = 1e-5)
})

test_that("Gaussian with equal variances", {
  ours <- mixture_regression(yg ~ x, data_for_comparison, 2, variance = "equal",
                 n_starts = 5, seed = 1, tol = 1e-13, vcov_type = "none")
  theirs <- best_flexmix(yg ~ x, data = data_for_comparison, k = 2,
                         model = flexmix::FLXMRglmfix(varFix = TRUE))
  expect_gte(ours$log_likelihood, theirs@logLik - 1e-7)
})

test_that("the parameter count is flexmix's", {
  ours <- mixture_regression(yb ~ x, data_for_comparison, 2, family = "binomial",
                 id = "id", class_level = "group", membership = ~ w,
                 n_starts = 2, seed = 1, vcov_type = "none")
  theirs <- best_flexmix(cbind(yb, 1 - yb) ~ x | id,
                         data = data_for_comparison, k = 2,
                         model = flexmix::FLXMRglm(family = "binomial"),
                         concomitant = flexmix::FLXPmultinom(~ w))
  expect_identical(ours$n_parameters, as.integer(theirs@df))
})
