# multilpa(family = "additive"): estimation. References are closed-form MLEs,
# numerical derivatives of the dense likelihood, and direct likelihoods that
# share no code with the EM.

additive_fit_quiet <- function(...) {
  withCallingHandlers(multilpa(...), latents_boundary = \(w) invokeRestart("muffleWarning"))
}

# Balanced one-class data built so the between variance is interior or on the
# boundary, with the closed-form maximum likelihood estimates.
additive_balanced <- function(shifts) {
  data <- data.frame(group = rep(seq_along(shifts), each = 4L),
                     y = rep(shifts, each = 4L) + c(-1, -0.3, 0.3, 1))
  averages <- tapply(data$y, data$group, mean)
  mu <- mean(averages)
  scatter <- sum(tapply(data$y, data$group, \(v) sum((v - mean(v))^2)))
  spread <- sum((averages - mu)^2)
  within <- scatter / (length(shifts) * 3)
  between <- spread / length(shifts) - within / 4
  if (between < 0) {
    between <- 0
    within <- (scatter + 4 * spread) / (length(shifts) * 4)
  }
  list(data = data, mean = mu, within = within, between = between)
}

test_that("one class, balanced groups: interior MLE matches the closed form", {
  case <- additive_balanced(seq(-2, 3, length.out = 6L))
  expect_gt(case$between, 0)
  fit <- multilpa(case$data, "y", "group", n_group_classes = 1,
                  family = "additive", tol = 1e-13, max_iter = 5000, seed = 1)
  expect_equal(unname(fit$means[1, 1]), case$mean, tolerance = 1e-8)
  expect_equal(unname(fit$within_variances[1, 1]), case$within, tolerance = 1e-6)
  expect_equal(unname(fit$between_variances[1, 1]), case$between, tolerance = 1e-6)
  expect_false(any(fit$between_zero))
  expect_true(fit$converged)
})

test_that("one class, balanced groups: boundary MLE is reached exactly", {
  case <- additive_balanced(seq(-0.05, 0.05, length.out = 6L))
  expect_identical(case$between, 0)
  expect_warning(
    fit <- multilpa(case$data, "y", "group", n_group_classes = 1,
                    family = "additive", tol = 1e-12, seed = 1),
    class = "latents_boundary")
  expect_identical(unname(fit$between_variances[1, 1]), 0)
  expect_true(fit$between_zero[1, 1])
  expect_true(fit$kkt)
  expect_equal(unname(fit$means[1, 1]), case$mean, tolerance = 1e-10)
  expect_equal(unname(fit$within_variances[1, 1]), case$within, tolerance = 1e-10)
  # KKT: the likelihood falls as the between variance leaves zero.
  stats <- .additive_prepare(case$data, "y", "group")
  point <- list(means = fit$means, between = fit$between_variances,
                within = fit$within_variances, weights = fit$group_probabilities)
  score <- .additive_between_score(stats, point,
                                   .additive_expectation(stats, point)$posterior)
  expect_lt(score[1, 1], 0)
  # A zero-between variance is still counted: estimated, not fixed.
  expect_identical(fit$n_parameters, 3L)
})

additive_small <- local({
  truth <- list(means = rbind(c(-1, 0.5), c(1.2, -0.7)),
                between = rbind(c(0.3, 0.4), c(0.5, 0.2)),
                within = c(0.8, 1.1), weights = c(0.45, 0.55))
  additive_draw(truth, rep(c(2L, 3L, 5L, 6L), 15L), seed = 20260930)
})

# Dense-covariance log likelihood on an unconstrained vector: means, log W,
# log T, baseline class logits. Shares no code with the package.
additive_dense_loglik <- function(theta, blocks, n_classes, d, between_variance) {
  n_between <- if (between_variance == "equal") 1L else n_classes
  means <- matrix(theta[seq_len(n_classes * d)], n_classes, d)
  within <- exp(theta[n_classes * d + seq_len(d)])
  between <- matrix(exp(theta[n_classes * d + d + seq_len(n_between * d)]),
                    n_between, d)[rep(seq_len(n_between), length.out = n_classes), ,
                                  drop = FALSE]
  # Class 1 is the baseline: logits are log(w_h / w_1) for h = 2..H.
  logits <- c(0, utils::tail(theta, n_classes - 1L))
  log_weights <- logits - log(sum(exp(logits)))
  sum(vapply(blocks, \(x) {
    joint <- vapply(seq_len(n_classes), \(h) log_weights[h] +
      additive_dense_density(x, means[h, ], between[h, ], within), numeric(1))
    max(joint) + log(sum(exp(joint - max(joint))))
  }, numeric(1)))
}

additive_theta <- function(fit) {
  n_classes <- nrow(fit$means)
  between <- if (fit$between_variance == "equal") fit$between_variances[1, ] else
    as.vector(fit$between_variances)
  c(as.vector(fit$means), log(fit$within_variances[1, ]), log(between),
    log(fit$group_probabilities[-1L] / fit$group_probabilities[1L]))
}

test_that("the fitted point is stationary for the dense likelihood", {
  skip_if_not_installed("numDeriv")
  blocks <- additive_blocks(additive_small, c("y1", "y2"), "group")
  lapply(c("varying", "equal"), \(between_variance) {
    fit <- additive_fit_quiet(additive_small, c("y1", "y2"), "group",
                              n_group_classes = 2, family = "additive",
                              between_variance = between_variance,
                              tol = 1e-12, max_iter = 5000, seed = 7)
    expect_false(any(fit$between_zero))
    theta <- additive_theta(fit)
    expect_equal(additive_dense_loglik(theta, blocks, 2L, 2L, between_variance),
                 fit$log_likelihood, tolerance = 1e-10)
    gradient <- numDeriv::grad(additive_dense_loglik, theta, blocks = blocks,
                               n_classes = 2L, d = 2L,
                               between_variance = between_variance)
    expect_lt(max(abs(gradient)), 1e-3)
  })
})

test_that("the between-variance score equals a numerical derivative", {
  skip_if_not_installed("numDeriv")
  stats <- .additive_prepare(additive_small, c("y1", "y2"), "group")
  point <- list(means = rbind(c(-0.5, 0.2), c(1, -0.4)),
                between = rbind(c(0.2, 0.6), c(0.4, 0.1)),
                within = c(0.9, 1.2), weights = c(0.3, 0.7))
  analytic <- .additive_between_score(
    stats, point, .additive_expectation(stats, point)$posterior)
  numeric_score <- matrix(numDeriv::grad(\(v) {
    moved <- point
    moved$between[] <- v
    .additive_expectation(stats, moved)$log_likelihood
  }, as.vector(point$between)), 2L)
  expect_equal(analytic, numeric_score, tolerance = 1e-7)
  expect_equal(.additive_between_score(
    stats, point, .additive_expectation(stats, point)$posterior, "equal"),
    colSums(numeric_score), tolerance = 1e-7)
})

test_that("with between variances held at zero, EM matches the group-level mixture", {
  skip_if_not_installed("numDeriv")
  stats <- .additive_prepare(additive_small, c("y1", "y2"), "group")
  held <- matrix(TRUE, 2L, 2L)
  structure <- .additive_structure("additive", "varying")
  start <- .additive_initialize(stats, 2L, structure, 1e-6, 1L)
  start$between[] <- 0
  fit <- .additive_em(stats, start, structure, 1e-6, 5000L, 1e-13, held)
  expect_true(fit$converged)
  expect_true(all(diff(fit$history) >= -1e-9))
  blocks <- additive_blocks(additive_small, c("y1", "y2"), "group")
  # Every rating of a group shares the group's class; ratings are otherwise
  # independent Gaussians with class means and shared variances.
  mixture <- function(theta) {
    means <- matrix(theta[1:4], 2L)
    sd <- sqrt(exp(theta[5:6]))
    log_weights <- log(c(1, exp(theta[7])) / (1 + exp(theta[7])))
    sum(vapply(blocks, \(x) {
      joint <- vapply(1:2, \(h) log_weights[h] + sum(stats::dnorm(
        x, rep(means[h, ], each = nrow(x)), rep(sd, each = nrow(x)), log = TRUE)),
        numeric(1))
      max(joint) + log(sum(exp(joint - max(joint))))
    }, numeric(1)))
  }
  theta <- c(as.vector(fit$parameters$means), log(fit$parameters$within[1, ]),
             log(fit$parameters$weights[2] / fit$parameters$weights[1]))
  expect_equal(mixture(theta), fit$expectation$log_likelihood, tolerance = 1e-10)
  expect_lt(max(abs(numDeriv::grad(mixture, theta))), 1e-3)
  posterior_reference <- t(vapply(blocks, \(x) {
    joint <- vapply(1:2, \(h) log(fit$parameters$weights[h]) + sum(stats::dnorm(
      x, rep(fit$parameters$means[h, ], each = nrow(x)),
      rep(sqrt(fit$parameters$within[1, ]), each = nrow(x)), log = TRUE)), numeric(1))
    exp(joint - max(joint)) / sum(exp(joint - max(joint)))
  }, numeric(2)))
  expect_equal(unname(fit$expectation$posterior), posterior_reference,
               tolerance = 1e-10)
})

test_that("EM is monotone and equal between variances are shared", {
  fit <- additive_fit_quiet(additive_small, c("y1", "y2"), "group",
                            n_group_classes = 2, family = "additive",
                            between_variance = "equal", seed = 3)
  expect_true(all(diff(fit$log_likelihood_history) >=
                    -1e-10 * (1 + abs(fit$log_likelihood))))
  expect_equal(fit$between_variances[1, ], fit$between_variances[2, ])
  expect_identical(fit$n_parameters, 9L)
  varying <- additive_fit_quiet(additive_small, c("y1", "y2"), "group",
                                n_group_classes = 2, family = "additive", seed = 3)
  expect_identical(varying$n_parameters, 11L)
  # Nested: the varying model cannot fit worse than the equal one.
  expect_gte(varying$log_likelihood, fit$log_likelihood - 1e-6)
  expect_equal(fit$bic, -2 * fit$log_likelihood + 9 * log(60))
  expect_equal(fit$bic_individual,
               -2 * fit$log_likelihood + 9 * log(nrow(additive_small)))
  expect_equal(unname(rowSums(fit$group_posteriors)), rep(1, 60), tolerance = 1e-12)
  expect_identical(fit$n_groups, 60L)
})

test_that("seeded fits reproduce and leave the caller's RNG untouched", {
  set.seed(99)
  before <- .Random.seed
  a <- additive_fit_quiet(additive_small, c("y1", "y2"), "group",
                          n_group_classes = 2, family = "additive",
                          n_starts = 4, seed = 11)
  expect_identical(.Random.seed, before)
  b <- additive_fit_quiet(additive_small, c("y1", "y2"), "group",
                          n_group_classes = 2, family = "additive",
                          n_starts = 4, seed = 11)
  expect_identical(a$log_likelihood, b$log_likelihood)
  expect_identical(a$starts, b$starts)
  expect_identical(nrow(a$starts), 4L)
  expect_identical(sum(a$starts$selected), 1L)
})

test_that("likelihood is invariant to row order", {
  a <- additive_fit_quiet(additive_small, c("y1", "y2"), "group",
                          n_group_classes = 2, family = "additive", seed = 5)
  shuffled <- additive_small[rev(seq_len(nrow(additive_small))), ]
  b <- additive_fit_quiet(shuffled, c("y1", "y2"), "group",
                          n_group_classes = 2, family = "additive", seed = 5)
  expect_equal(b$log_likelihood, a$log_likelihood, tolerance = 1e-8)
})

test_that("inapplicable arguments and unidentified requests are refused", {
  vars <- c("y1", "y2")
  expect_error(multilpa(additive_small, vars, "group", n_profiles = 2,
                        family = "additive"), class = "latents_bad_argument")
  expect_error(multilpa(additive_small, vars, "group", n_group_classes = 2,
                        family = "additive", categorical = "y1"),
               class = "latents_bad_argument")
  expect_error(multilpa(additive_small, vars, family = "additive"),
               class = "latents_bad_argument")
  expect_error(multilpa(additive_small, vars, NULL, family = "additive"),
               class = "latents_bad_argument")
  expect_error(lpa(additive_small, vars, n_profiles = 2, family = "additive"),
               class = "latents_bad_argument")
  expect_error(multilpa(additive_small, vars, "group", n_profiles = 2,
                        between_variance = "equal"),
               class = "latents_bad_argument")
  expect_error(enumerate_classes(additive_small, vars, "group", n_profiles = 2,
                                 family = "additive"),
               class = "latents_bad_argument")
  expect_error(multilpa(additive_small, vars, "group", n_group_classes = 61,
                        family = "additive"), class = "latents_unidentified")
  singles <- additive_small[!duplicated(additive_small$group), ]
  expect_error(multilpa(singles, vars, "group", n_group_classes = 2,
                        family = "additive"), class = "latents_unidentified")
  expect_error(multilpa(additive_small, vars, "group", n_group_classes = 0,
                        family = "additive"), class = "latents_bad_argument")
})

test_that("the print method is short and stable", {
  fit <- additive_fit_quiet(additive_small, c("y1", "y2"), "group",
                            n_group_classes = 2, family = "additive", seed = 3)
  expect_snapshot(print(fit))
  expect_s3_class(fit, "multilpa_additive")
  expect_false(inherits(fit, "multilpa"))
})
