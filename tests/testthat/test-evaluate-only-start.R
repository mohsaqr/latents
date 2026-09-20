# `max_iter = 0` with a supplied `start` is a request to score that parameter
# set. No update is performed, so no restart could improve on it; evaluating
# random initializations as well and keeping the highest would answer a
# different question with a different number. These tests fix the promise that
# the values handed in are the values scored, whatever `n_starts` says.

.evaluate_frame <- function(seed = 127L, n_groups = 30L, size = 8L) {
  set.seed(seed)
  n <- n_groups * size
  data.frame(school = rep(seq_len(n_groups), each = size),
             a = stats::rnorm(n), b = stats::rnorm(n, 5, 2))
}

# A deliberately implausible one-profile solution: far from anything EM would
# find, so a random start beating it would be unmistakable.
.evaluate_start <- function(mean_value = 20) {
  list(means = matrix(c(mean_value, mean_value), 1L),
       variances = matrix(c(1, 1), 1L),
       profile_probabilities = matrix(1, 1L, 1L),
       group_probabilities = 1)
}

test_that("max_iter = 0 scores the supplied start, and the closed form agrees", {
  frame <- .evaluate_frame()
  start <- .evaluate_start()
  evaluated <- multilpa(frame, c("a", "b"), "school", n_profiles = 1,
                        n_group_classes = 1, start = start, max_iter = 0,
                        n_starts = 1, seed = 1)
  # One profile, one group class, unit variances: the likelihood is just the
  # independent normal density, written from the definition rather than from
  # the package's own machinery.
  reference <- sum(stats::dnorm(frame$a, 20, 1, log = TRUE)) +
    sum(stats::dnorm(frame$b, 20, 1, log = TRUE))
  expect_equal(as.numeric(logLik(evaluated)), reference, tolerance = 1e-10)
  expect_identical(evaluated$iterations, 0L)
  expect_equal(evaluated$means, start$means, ignore_attr = TRUE)
})

test_that("evaluate-only ignores n_starts and never substitutes a random start", {
  frame <- .evaluate_frame()
  start <- .evaluate_start()
  single <- multilpa(frame, c("a", "b"), "school", n_profiles = 1,
                     n_group_classes = 1, start = start, max_iter = 0,
                     n_starts = 1, seed = 1)
  # The default n_starts is 10. Before this was fixed, nine random starts were
  # also evaluated and the best of them replaced the supplied values, returning
  # a much higher likelihood for parameters the caller never supplied.
  defaulted <- multilpa(frame, c("a", "b"), "school", n_profiles = 1,
                        n_group_classes = 1, start = start, max_iter = 0,
                        seed = 1)
  many <- multilpa(frame, c("a", "b"), "school", n_profiles = 1,
                   n_group_classes = 1, start = start, max_iter = 0,
                   n_starts = 25, seed = 1)
  expect_equal(defaulted$log_likelihood, single$log_likelihood, tolerance = 1e-12)
  expect_equal(many$log_likelihood, single$log_likelihood, tolerance = 1e-12)
  expect_equal(defaulted$means, start$means, ignore_attr = TRUE)
  expect_equal(many$means, start$means, ignore_attr = TRUE)
  # Exactly one start was run, and the result says so rather than hiding it.
  expect_identical(defaulted$best_start, 1L)
  expect_identical(nrow(as.data.frame(defaulted, what = "starts")), 1L)
  expect_identical(nrow(as.data.frame(many, what = "starts")), 1L)
})

test_that("a start supplied for a fit that does iterate still only seeds the first start", {
  frame <- .evaluate_frame()
  start <- .evaluate_start()
  # max_iter > 0 leaves the multiple-start search alone: the supplied values are
  # the first start and the remaining ones are random, as documented.
  searched <- multilpa(frame, c("a", "b"), "school", n_profiles = 1,
                       n_group_classes = 1, start = start, max_iter = 500,
                       n_starts = 3, seed = 1)
  expect_identical(nrow(as.data.frame(searched, what = "starts")), 3L)
  expect_gt(searched$log_likelihood, -76708.5)
})

test_that("max_iter = 0 without a start still scores every requested start", {
  frame <- .evaluate_frame()
  # There is nothing to evaluate at, so scoring the random initializations is
  # the only thing the call can mean, and it is left as it was.
  scored <- multilpa(frame, c("a", "b"), "school", n_profiles = 2,
                     n_group_classes = 2, max_iter = 0, n_starts = 4, seed = 2)
  expect_identical(nrow(as.data.frame(scored, what = "starts")), 4L)
  expect_identical(scored$iterations, 0L)
})

test_that("evaluate-only with a held measurement returns the held values", {
  frame <- .evaluate_frame()
  stage_one <- multilpa(frame, c("a", "b"), "school", n_profiles = 2,
                        n_group_classes = 1, n_starts = 3, seed = 1)
  start <- starting_values(stage_one, what = "measurement")
  evaluated <- multilpa(frame, c("a", "b"), "school", n_profiles = 2,
                        n_group_classes = 2, n_starts = 4, seed = 3,
                        max_iter = 0, start = start, fixed = "measurement")
  expect_equal(evaluated$means, stage_one$means)
  expect_equal(evaluated$variances, stage_one$variances)
  expect_identical(evaluated$iterations, 0L)
  expect_identical(evaluated$best_start, 1L)
  # The mixing values come from the start too, not from a random restart.
  expect_equal(evaluated$profile_probabilities,
               matrix(0.5, 2L, 2L), ignore_attr = TRUE)
})
