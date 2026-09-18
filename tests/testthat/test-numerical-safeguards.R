test_that("failed starts are counted and do not hide a valid restart", {
  set.seed(381)
  synthetic <- data.frame(group = rep(seq_len(10), each = 8),
                          x = rnorm(80), y = rnorm(80))
  bad_start <- list(means = matrix(1e200, 1L, 2L),
                    variances = matrix(1, 1L, 2L),
                    profile_probabilities = matrix(1, 1L, 1L),
                    group_probabilities = 1)
  expect_warning(
    fit <- fit_multilpa(synthetic, c("x", "y"), "group", 1, 1,
                      n_starts = 2, start = bad_start, seed = 42),
    "1 of 2 starts failed")
  expect_equal(fit$n_failed_starts, 1L)
  expect_match(fit$starts$error[1L], "densities vanished|overflowed")
  expect_equal(fit$starts$log_likelihood[1L], -Inf)
  expect_equal(fit$best_start, 2L)
  expect_true(fit$converged)
  expect_error(fit_multilpa(synthetic, c("x", "y"), "group", 1, 1,
                          n_starts = 1, start = bad_start), "All 1 starts failed")
})

test_that("three profiles and group classes retain all matrix dimensions", {
  set.seed(809)
  group_index <- rep(seq_len(60L), each = 15L)
  true_groups <- rep(1:3, each = 20L)
  probabilities <- matrix(0.1, 3L, 3L)
  diag(probabilities) <- 0.8
  profiles <- vapply(group_index, function(group) {
    stopifnot(length(group) == 1L, group >= 1L)
    sample.int(3L, 1L, prob = probabilities[true_groups[group], ])
  }, integer(1))
  synthetic <- data.frame(group = group_index,
                          score = rnorm(length(profiles), c(-4, 0, 4)[profiles], 0.4))
  fit <- fit_multilpa(synthetic, "score", "group", 3L, 3L,
                    n_starts = 4L, seed = 998)
  expect_true(fit$converged)
  expect_equal(dim(fit$means), c(3L, 1L))
  expect_equal(dim(fit$variances), c(3L, 1L))
  expect_equal(dim(fit$profile_probabilities), c(3L, 3L))
  expect_equal(dim(fit$subject_posteriors), c(900L, 3L))
  expect_equal(dim(fit$group_posteriors), c(60L, 3L))
  expect_equal(fit$n_parameters, 14)
  expect_equal(sort(as.numeric(fit$means)), c(-4, 0, 4), tolerance = 0.1)
  expect_equal(unname(rowSums(fit$profile_probabilities)), rep(1, 3L), tolerance = 1e-10)
  expect_equal(sum(fit$effective_group_counts), 60, tolerance = 1e-10)
  expect_equal(sum(fit$effective_profile_counts), 900, tolerance = 1e-10)
  expect_true(all(diff(fit$log_likelihood_history) >= -1e-7))
})

test_that("an unbalanced-group update agrees with exhaustive assignments", {
  data <- data.frame(group = rep(c("large", "singleton", "middle", "pair"),
                                 times = c(4, 1, 3, 2)),
                     x = c(-1.2, 0.8, -0.6, 1.2, 0.3, -0.9, 0.9, 1.1, -0.5, 0.7))
  parameters <- list(means = matrix(c(-0.7, 0.8), 2L, 1L),
                      variances = matrix(c(0.7, 0.8), 2L, 1L),
                      profile_probabilities = rbind(c(0.6, 0.4), c(0.4, 0.6)),
                      group_probabilities = c(0.5, 0.5))
  reference <- enumerated_em_update(as.matrix(data["x"]), data$group, parameters)
  expect_warning(
    fit <- fit_multilpa(data, "x", "group", 2L, 2L, n_starts = 1L,
                      start = parameters, max_iter = 1L), "did not converge")
  expect_equal(unname(fit$group_probabilities), reference$group_probabilities,
               tolerance = 1e-10)
  expect_equal(unname(fit$profile_probabilities), reference$profile_probabilities,
               tolerance = 1e-10)
  expect_equal(unname(fit$means), unname(reference$means), tolerance = 1e-10)
  expect_equal(unname(fit$variances), unname(reference$variances), tolerance = 1e-10)
})
