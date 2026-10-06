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

test_that("the likelihood is the mixture likelihood, in every nesting", {
  skip_on_cran()
  rows <- simulate_rows()
  single <- mixture_regression(y ~ x, rows, 2, n_starts = 1, seed = 1, vcov_type = "none")
  expect_equal(single$log_likelihood, direct_log_likelihood(single),
               tolerance = 1e-10)
  concomitant <- mixture_regression(y ~ x, rows, 2, membership = ~ w, n_starts = 1,
                                    seed = 1, vcov_type = "none")
  expect_equal(concomitant$log_likelihood,
               direct_log_likelihood(concomitant), tolerance = 1e-10)
  grouped <- mixture_regression(passed ~ hours, study_hours, 2, family = "binomial",
                                id = "student", class_level = "group", n_starts = 1,
                                seed = 1, vcov_type = "none")
  expect_equal(grouped$log_likelihood, direct_log_likelihood(grouped),
               tolerance = 1e-10)
  two_level <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                                  n_group_classes = 2, membership = ~ sleep,
                                  group_membership = ~ motivation, n_starts = 1,
                                  seed = 1, vcov_type = "none")
  expect_equal(two_level$log_likelihood, direct_log_likelihood(two_level),
               tolerance = 1e-10)
  counts <- mixture_regression(count ~ x, rows, 2, family = "poisson", n_starts = 1,
                               seed = 1, vcov_type = "none")
  expect_equal(counts$log_likelihood, direct_log_likelihood(counts),
               tolerance = 1e-10)
})

test_that("a two-class Gaussian mixture recovers the generating regressions", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, n_starts = 1, seed = 1)
  table <- get_results(fit, "coefficients")
  truth <- c(35, 4.5, 55, 0.8)
  expect_true(all(abs(table$estimate - truth) < 3 * table$std_error + 0.5))
  classes <- get_results(fit, "classes")
  expect_equal(classes$sigma, c(5, 7), tolerance = 0.15)
})

test_that("standard errors match a numerical Hessian of the likelihood", {
  skip_on_cran()
  rows <- simulate_rows()
  check <- function(fit) {
    theta <- latents:::.mixture_pack(fit$spec, fit$params)
    log_lik <- function(t) {
      params <- latents:::.mixture_unpack(fit$spec, t, fit$params)
      latents:::.mixture_expectation(fit$spec, params)$log_likelihood
    }
    h <- 1e-4
    size <- length(theta)
    unit <- diag(h, size)
    pairs <- expand.grid(i = seq_len(size), j = seq_len(size))
    hessian <- matrix(vapply(seq_len(nrow(pairs)), function(cell) {
      e_i <- unit[, pairs$i[cell]]
      e_j <- unit[, pairs$j[cell]]
      (log_lik(theta + e_i + e_j) - log_lik(theta + e_i - e_j) -
         log_lik(theta - e_i + e_j) + log_lik(theta - e_i - e_j)) / (4 * h^2)
    }, numeric(1)), size)
    numerical <- sqrt(diag(solve(-hessian)))
    analytic <- sqrt(diag(fit$inference$vcov))
    expect_lt(max(abs(analytic / numerical - 1)), 1e-4)
  }
  check(mixture_regression(y ~ x, rows, 2, membership = ~ w, n_starts = 1, seed = 1,
                           tol = 1e-12))
  check(mixture_regression(y ~ x + w, rows, 2, common = ~ w, variance = "equal",
                           n_starts = 1, seed = 1, tol = 1e-12))
  check(mixture_regression(count ~ x, rows, 2, family = "poisson", n_starts = 1, seed = 1,
                           tol = 1e-12))
  check(mixture_regression(passed ~ hours, study_hours, 2, family = "binomial",
                           id = "student", class_level = "group",
                           membership = ~ motivation, n_starts = 1, seed = 1, tol = 1e-12))
  check(mixture_regression(score ~ hours, study_hours, 2, id = "student",
                           n_group_classes = 2, membership = ~ sleep,
                           group_membership = ~ motivation, n_starts = 1, seed = 1,
                           tol = 1e-12))
})

test_that("the score vanishes at the estimate", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                            n_group_classes = 2, n_starts = 1, seed = 1, tol = 1e-12)
  expect_lt(fit$inference$score_norm, 1e-2)
})

test_that("the robust covariance is the sandwich clustered on the unit", {
  fit <- mixture_regression(y ~ x, simulate_rows(), 2, n_starts = 1, seed = 1)
  robust <- vcov(fit, type = "robust")
  information <- fit$inference$information
  scores <- latents:::.mixture_unit_scores(fit$spec, fit$params,
                                          fit$expectation)
  bread <- solve(information)
  expect_equal(unname(robust), unname(bread %*% crossprod(scores) %*% bread),
               tolerance = 1e-8)
  # A single-level fit given `id` clusters on it, which changes the meat.
  clustered_rows <- cbind(simulate_rows(), cluster = rep(seq_len(60), 5))
  clustered <- mixture_regression(y ~ x, clustered_rows, 2, id = "cluster",
                                  n_starts = 1, seed = 1, vcov_type = "robust")
  expect_identical(clustered$inference$clustered_on, "cluster")
})

test_that("classes are ordered by size and the order is invariant to starts", {
  rows <- simulate_rows()
  one <- mixture_regression(y ~ x, rows, 2, n_starts = 2, seed = 1, vcov_type = "none")
  two <- mixture_regression(y ~ x, rows, 2, n_starts = 6, seed = 99, vcov_type = "none")
  shares <- get_results(one, "classes")$share
  expect_true(all(diff(shares) <= 0))
  expect_equal(one$params$beta, two$params$beta, tolerance = 1e-5)
})

test_that("row order does not change the fit", {
  rows <- simulate_rows()
  shuffled <- rows[rev(seq_len(nrow(rows))), ]
  original <- mixture_regression(y ~ x, rows, 2, n_starts = 1, seed = 1,
                                 vcov_type = "none")
  reversed <- mixture_regression(y ~ x, shuffled, 2, n_starts = 1, seed = 1,
                                 vcov_type = "none")
  expect_equal(original$log_likelihood, reversed$log_likelihood,
               tolerance = 1e-7)
})

test_that("posteriors and priors are probability distributions", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                            n_group_classes = 2, membership = ~ sleep, n_starts = 1,
                            seed = 1, vcov_type = "none")
  assignments <- get_results(fit, "assignments")
  probabilities <- assignments[c("probability_class_1", "probability_class_2")]
  expect_equal(unname(rowSums(probabilities)), rep(1, nrow(study_hours)),
               tolerance = 1e-12)
  classes <- predict(fit, type = "class_response")
  by_row <- tapply(classes$prior, classes$row, sum)
  expect_equal(unname(as.vector(by_row)), rep(1, nrow(study_hours)),
               tolerance = 1e-12)
})

test_that("a seed makes the fit reproducible and restores the caller's state", {
  set.seed(42)
  before <- stats::runif(1)
  set.seed(42)
  first <- mixture_regression(y ~ x, simulate_rows(), 2, n_starts = 3, seed = 7,
                              vcov_type = "none")
  after <- stats::runif(1)
  expect_identical(before, after)
  second <- mixture_regression(y ~ x, simulate_rows(), 2, n_starts = 3, seed = 7,
                               vcov_type = "none")
  expect_identical(first$starts, second$starts)
})

test_that("predict() on the fitted rows agrees with the fitted table", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, n_starts = 1, seed = 1)
  posterior <- predict(fit, study_hours, type = "posterior")
  expect_equal(posterior$probability_class_1,
               get_results(fit, "assignments")$probability_class_1,
               tolerance = 1e-12)
  response <- predict(fit, data.frame(hours = c(0, 10)))
  classes <- predict(fit, data.frame(hours = c(0, 10)), type = "class_response")
  expected <- tapply(classes$prior * classes$fitted, classes$row, sum)
  expect_equal(response$fitted, unname(as.vector(expected)))
})

test_that("simulate() reproduces the fitted class structure", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, n_starts = 1, seed = 1)
  draws <- simulate(fit, nsim = 2, seed = 3)
  expect_identical(dim(draws), c(nrow(study_hours), 2L))
  expect_identical(draws, simulate(fit, nsim = 2, seed = 3))
  refit_data <- study_hours
  refit_data$score <- draws$sim_1
  refit <- mixture_regression(score ~ hours, refit_data, 2, n_starts = 1, seed = 1,
                              vcov_type = "none")
  expect_equal(unname(refit$params$beta), unname(fit$params$beta),
               tolerance = 0.15)
})

test_that("shared coefficients are equal across classes and counted once", {
  fit <- mixture_regression(score ~ hours + sleep, study_hours, 2, common = ~ sleep,
                            n_starts = 1, seed = 1)
  table <- get_results(fit, "coefficients")
  expect_identical(sum(table$term == "sleep"), 1L)
  expect_identical(table$class[table$term == "sleep"], "common")
  expect_identical(fit$n_parameters, 2L * 2L + 1L + 2L + 1L)
})

test_that("enumeration returns one tidy row per model and marks the BIC choice", {
  skip_on_cran()
  classes <- enumerate_regressions(score ~ hours, study_hours, n_classes = 1:3,
                                   n_starts = 1, seed = 1)
  table <- as.data.frame(classes)
  expect_identical(table$n_classes, 1:3)
  expect_identical(sum(table$best_bic), 1L)
  expect_identical(table$n_classes[table$best_bic], 2L)
  chosen <- get_results(classes, "model", n_classes = 2)
  expect_s3_class(chosen, "latents_mixture_regression")
  expect_error(get_results(classes, "model", n_classes = 5),
               class = "latents_bad_argument")
})

test_that("the bootstrap test withholds a p-value when refits fail validation", {
  skip_on_cran()
  expect_warning(classes <- enumerate_regressions(score ~ hours, study_hours, n_classes = 1:2,
                                   n_starts = 1, seed = 1, bootstrap = 9,
                                   bootstrap_starts = 1),
                  class = "latents_failed_replicates")
  table <- as.data.frame(classes)
  expect_true(is.na(table$blrt_p_value[1]))
  expect_true(is.na(table$blrt_p_value[2]))
  expect_gt(table$blrt_flagged[2], 0L)
  expect_equal(table$blrt_replicates[2] + table$blrt_flagged[2], 9L)
})

test_that("every tidy table is a data frame", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                            n_group_classes = 2, n_starts = 1, seed = 1)
  tables <- get_results(fit, "all")
  expect_true(all(vapply(tables, is.data.frame, logical(1))))
  expect_identical(nrow(tables$groups), 150L)
  expect_identical(as.data.frame(fit), tables$coefficients)
})

test_that("print and summary are stable", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, n_starts = 2, seed = 1)
  expect_snapshot(print(fit))
  expect_snapshot(print(summary(fit), digits = 3))
})

test_that("plots draw for every view", {
  skip_on_cran()
  fit <- mixture_regression(score ~ hours, study_hours, 2, n_starts = 1, seed = 1)
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, what = "coefficients"))
  expect_invisible(plot(fit, what = "posteriors"))
})

test_that("misspecified requests are refused by class", {
  expect_error(mixture_regression(passed ~ hours, study_hours, 2, family = "binomial", n_starts = 1),
                                  class = "latents_not_identified")
  expect_error(mixture_regression(score ~ hours, study_hours, 2, class_level = "group", n_starts = 1),
                                  class = "latents_bad_argument")
  expect_error(mixture_regression(score ~ hours, study_hours, 2, n_group_classes = 2, n_starts = 1),
                                  class = "latents_bad_argument")
  expect_error(mixture_regression(score ~ hours, study_hours, 2, id = "student",
                                  class_level = "group", membership = ~ sleep, n_starts = 1),
                                  class = "latents_bad_data")
  expect_error(mixture_regression(score ~ hours, study_hours, 2, common = ~ sleep, n_starts = 1),
                                  class = "latents_bad_argument")
  expect_error(mixture_regression(score ~ hours, study_hours, 2, family = "poisson", n_starts = 1),
                                  class = "latents_bad_data")
  expect_error(mixture_regression(score ~ hours, study_hours, 2, family = "poisson",
                                  variance = "equal", n_starts = 1),
                                  class = "latents_bad_argument")
  gappy <- study_hours
  gappy$hours[1:3] <- NA
  expect_error(mixture_regression(score ~ hours, gappy, 2, n_starts = 1), class = "latents_missing_data")
  expect_warning(mixture_regression(score ~ hours, gappy, 2, missing = "omit",
                                    n_starts = 1, seed = 1, vcov_type = "none"),
                                    class = "latents_rows_dropped")
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
  expect_equal(fit$log_likelihood, direct_log_likelihood(fit),
               tolerance = 1e-10)
  expect_equal(sort(get_results(fit)$estimate[c(2, 4)]), c(-1, 2),
               tolerance = 0.25)
})

test_that("the recovery table cross-tabulates assignments against the truth", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, n_starts = 1, seed = 1,
                            vcov_type = "none")
  recovery <- get_results(fit, "recovery", data = study_hours,
                          truth = "strategy")
  expect_identical(names(recovery), c("assigned", "strategy", "n", "share"))
  expect_identical(sum(recovery$n), nrow(study_hours))
  expect_equal(as.vector(tapply(recovery$share, recovery$assigned, sum)),
               c(1, 1))
  two_level <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                                  n_group_classes = 2, n_starts = 1, seed = 1,
                                  vcov_type = "none")
  by_student <- get_results(two_level, "recovery", data = study_hours,
                            truth = "student_type", by = "group_class")
  expect_identical(sum(by_student$n), 150L)
  expect_error(get_results(two_level, "recovery", data = study_hours,
                           truth = "strategy", by = "group_class"),
               class = "latents_bad_data")
  expect_error(get_results(fit, "recovery"), class = "latents_bad_argument")
})

test_that("the classification table's rows are probability distributions", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, n_starts = 1, seed = 1,
                            vcov_type = "none")
  classification <- get_results(fit, "classification")
  expect_identical(nrow(classification), 4L)
  by_assigned <- tapply(classification$mean_posterior,
                        classification$assigned, sum)
  expect_equal(as.vector(by_assigned), c(1, 1), tolerance = 1e-12)
  # Its diagonal is the class table's mean posterior.
  diagonal <- subset(classification, assigned == class)
  expect_equal(diagonal$mean_posterior,
               get_results(fit, "classes")$mean_posterior, tolerance = 1e-12)
})
