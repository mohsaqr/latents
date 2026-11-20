# mixture_regression(): likelihood, estimation, inference, accessors and refusals.

simulate_rows <- function(n = 300L, seed = 1L) {
  latents:::.mixture_with_seed(seed, {
    class <- 1L + (stats::runif(n) < 0.4)
    x <- stats::runif(n, 0, 10)
    w <- stats::rnorm(n)
    data.frame(
      x = x, w = w,
      y = ifelse(class == 1L, 1 + 2 * x, 10 - 0.5 * x) +
        stats::rnorm(n, sd = ifelse(class == 1L, 1, 2)),
      count = stats::rpois(n, exp(ifelse(class == 1L, 0.2 + 0.2 * x,
                                         1.5 - 0.1 * x))))
  })
}

# The likelihood written from its definition, one term per unit, sharing no
# code with the engine.
direct_log_likelihood <- function(fit) {
  spec <- fit$spec
  p <- fit$params
  eta <- spec$x %*% p$beta + spec$offset
  if (ncol(spec$z) > 0L) eta <- eta + drop(spec$z %*% p$common)
  density <- switch(spec$family,
    gaussian = vapply(seq_len(spec$n_classes), function(k) {
      stats::dnorm(spec$y, eta[, k], sqrt(p$sigma2[k]))
    }, numeric(spec$n)),
    poisson = vapply(seq_len(spec$n_classes), function(k) {
      stats::dpois(spec$y, exp(eta[, k]))
    }, numeric(spec$n)),
    binomial = vapply(seq_len(spec$n_classes), function(k) {
      stats::dbinom(spec$y, spec$trials, stats::plogis(eta[, k]))
    }, numeric(spec$n)))
  density <- matrix(density, spec$n)
  softmax <- function(scores) exp(scores) / rowSums(exp(scores))
  switch(spec$nesting,
    observation = sum(log(rowSums(softmax(spec$w %*% p$gamma) * density))),
    group = {
      prior <- softmax(spec$w %*% p$gamma)
      sum(vapply(seq_len(spec$n_groups), function(j) {
        rows <- spec$group_index == j
        log(sum(prior[j, ] * apply(density[rows, , drop = FALSE], 2L, prod)))
      }, numeric(1)))
    },
    "two-level" = {
      group_prior <- softmax(spec$v %*% p$delta)
      sum(vapply(seq_len(spec$n_groups), function(j) {
        rows <- which(spec$group_index == j)
        by_group_class <- vapply(seq_len(spec$n_group_classes), function(h) {
          indicator <- matrix(0, length(rows), spec$n_group_classes)
          indicator[, h] <- 1
          class_prior <- softmax(cbind(indicator, spec$w[rows, , drop = FALSE]) %*%
                                   p$class_logits)
          prod(rowSums(class_prior * density[rows, , drop = FALSE]))
        }, numeric(1))
        log(sum(group_prior[j, ] * by_group_class))
      }, numeric(1)))
    })
}

# One class's regression coefficients (or standard errors) from a vector named
# by coef(), renamed to the bare term labels lm() and glm() use.
class_terms <- function(values, k = 1L) {
  prefix <- sprintf("coefficient.class_%d.", k)
  kept <- values[startsWith(names(values), prefix)]
  stats::setNames(unname(kept), substring(names(kept), nchar(prefix) + 1L))
}

test_that("one class reduces to lm() and glm()", {
  rows <- simulate_rows()
  gaussian <- mixture_regression(y ~ x + w, rows, 1, n_starts = 0, vcov_type = "none")
  linear <- stats::lm(y ~ x + w, rows)
  expect_reference(as.numeric(logLik(gaussian)), as.numeric(stats::logLik(linear)),
                   reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one gaussian class", tolerance = 1e-10,
                   measure = "relative")
  expect_reference(class_terms(coef(gaussian)), stats::coef(linear), reference = "stats",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "coefficients, one gaussian class", tolerance = 1e-8,
                   measure = "relative")
  poisson <- mixture_regression(count ~ x, rows, 1, family = "poisson", n_starts = 0,
                                vcov_type = "none")
  reference <- stats::glm(count ~ x, stats::poisson(), rows)
  expect_reference(as.numeric(logLik(poisson)), as.numeric(stats::logLik(reference)),
                   reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one poisson class", tolerance = 1e-9,
                   measure = "relative")
  # glm()'s IRLS stops at its own deviance criterion (1e-8): 1e-7, the original tolerance.
  expect_reference(class_terms(coef(poisson)), stats::coef(reference), reference = "stats",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "coefficients, one poisson class", tolerance = 1e-7,
                   measure = "relative")
  # Standard errors of a one-class model are the GLM's. 1e-5, the original
  # tolerance: the observed information is differenced from analytic scores.
  fitted <- mixture_regression(count ~ x, rows, 1, family = "poisson", n_starts = 0)
  expect_reference(as.data.frame(get_results(fitted, "coefficients")),
                   data.frame(term = names(stats::coef(reference)),
                              std_error = unname(sqrt(diag(stats::vcov(reference))))),
                   reference = "stats",
                   covers = c("mixture_regression()", "get_results.latents_mixture_regression()"),
                   quantity = "standard errors, one poisson class", tolerance = 1e-5,
                   measure = "relative", key = "term")
})

# Moved from tests/testthat/test-mixture-regression.R: comparisons with independent references.

test_that("the likelihood is the mixture likelihood, in every nesting", {
  skip_on_cran()
  rows <- simulate_rows()
  # The same sum of logs computed twice: only the order of the additions differs.
  same_likelihood <- function(fit, setting) {
    expect_reference(as.numeric(logLik(fit)), direct_log_likelihood(fit),
                     reference = "independent: mixture likelihood from its definition, base R",
                     covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                     quantity = sprintf("log likelihood, %s", setting), tolerance = 1e-10,
                     measure = "relative")
  }
  single <- mixture_regression(y ~ x, rows, 2, n_starts = 1, seed = 1, vcov_type = "none")
  same_likelihood(single, "two gaussian classes")
  concomitant <- mixture_regression(y ~ x, rows, 2, membership = ~ w, n_starts = 1,
                                    seed = 1, vcov_type = "none")
  same_likelihood(concomitant, "two gaussian classes with a membership covariate")
  grouped <- mixture_regression(passed ~ hours, study_hours, 2, family = "binomial",
                                id = "student", class_level = "group", n_starts = 1,
                                seed = 1, vcov_type = "none")
  same_likelihood(grouped, "binomial, classes of groups")
  two_level <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                                  n_group_classes = 2, membership = ~ sleep,
                                  group_membership = ~ motivation, n_starts = 1,
                                  seed = 1, vcov_type = "none")
  same_likelihood(two_level, "two-level gaussian")
  counts <- mixture_regression(count ~ x, rows, 2, family = "poisson", n_starts = 1,
                               seed = 1, vcov_type = "none")
  same_likelihood(counts, "two poisson classes")
})

test_that("binomial counts with several trials are identified row by row", {
  rows <- latents:::.mixture_with_seed(5, {
    n <- 300
    class <- 1L + (stats::runif(n) < 0.5)
    x <- stats::rnorm(n)
    trials <- rep(8L, n)
    successes <- stats::rbinom(n, trials, stats::plogis(
      ifelse(class == 1L, -1 + 2 * x, 1 - x)))
    data.frame(x, successes, failures = trials - successes)
  })
  fit <- mixture_regression(cbind(successes, failures) ~ x, rows, 2, family = "binomial",
                            n_starts = 1, seed = 1)
  expect_reference(as.numeric(logLik(fit)), direct_log_likelihood(fit),
                   reference = "independent: mixture likelihood from its definition, base R",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, binomial counts of eight trials", tolerance = 1e-10,
                   measure = "relative")
  # Recovery of the simulated slopes: a unit check, not a reference comparison.
  expect_equal(sort(get_results(fit)$estimate[c(2, 4)]), c(-1, 2),
               tolerance = 0.25)
})
