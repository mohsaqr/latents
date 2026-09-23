# A bootstrap likelihood-ratio test on fits made with `fixed` is only a test if
# the simulated statistics come from the same pair of constrained models as the
# observed one. Comparing p-values cannot detect a violation: a bootstrap that
# silently drops the constraint still returns a finite, plausible-looking
# number. These tests therefore record what the refits are actually given.

# Two nested constrained fits sharing one measurement solution, plus the data
# they were fitted to.
constrained_pair <- function(seed = 984) {
  set.seed(seed)
  group <- rep(seq_len(40), each = 12)
  high <- rbinom(length(group), 1, rep(rep(c(0.1, 0.9), 20), each = 12))
  data <- data.frame(g = group,
                     a = rnorm(length(group), high * 4),
                     b = rnorm(length(group), high * 3))
  measurement <- multilpa(data, c("a", "b"), "g", 2, 1, n_starts = 3, seed = 5,
                          tol = 1e-10)
  start <- starting_values(measurement, what = "measurement")
  list(data = data, start = start,
       null_model = multilpa(data, c("a", "b"), "g", 2, 1, n_starts = 1,
                             start = start, fixed = "measurement", seed = 5,
                             tol = 1e-10),
       alternative_model = multilpa(data, c("a", "b"), "g", 2, 2, n_starts = 3,
                                    start = start, fixed = "measurement",
                                    seed = 5, tol = 1e-10))
}

# Records every argument each bootstrap refit receives, then delegates to the
# real fitter so the bootstrap still completes. Comparing p-values cannot show
# that the constraint survived; this can.
refit_recorder <- function() {
  calls <- list()
  fit <- multilpa
  list(mock = function(...) {
         arguments <- list(...)
         calls[[length(calls) + 1L]] <<- arguments
         do.call(fit, arguments)
       },
       calls = function() calls)
}

test_that("every bootstrap refit carries the held blocks and their values", {
  skip_on_cran()
  models <- constrained_pair()
  recorder <- refit_recorder()
  testthat::local_mocked_bindings(multilpa = recorder$mock)
  result <- bootstrap_lrt(models$null_model, models$alternative_model,
                          models$data, iter = 2, n_starts = 2, seed = 45,
                          tol = 1e-8)
  calls <- recorder$calls()
  # Two refits per replicate, null and alternative.
  expect_identical(length(calls), 4L)
  held <- vapply(calls, function(arguments)
    paste(arguments$fixed, collapse = ","), character(1))
  expect_identical(held, rep("means,variances", 4L))
  # The values held, not merely the block names, must be the observed fit's.
  expect_true(all(vapply(calls, function(arguments)
    isTRUE(all.equal(arguments$start$means, models$start$means)) &&
      isTRUE(all.equal(arguments$start$variances, models$start$variances)),
    logical(1))))
  # And the refits must compare the same free-parameter counts as the observed
  # pair, which is what the dropped constraint changed.
  expect_identical(result$fixed, c("means", "variances"))
  expect_equal(models$null_model$n_parameters, 1)
  expect_equal(models$alternative_model$n_parameters, 3)
})

test_that("an unconstrained comparison is refitted without a constraint", {
  skip_on_cran()
  set.seed(53)
  data <- data.frame(g = rep(seq_len(30), each = 8), y = rnorm(240))
  small <- multilpa(data, "y", "g", 1, 1, variance_model = "equal",
                    n_starts = 2, seed = 8)
  large <- multilpa(data, "y", "g", 2, 1, variance_model = "equal",
                    n_starts = 2, seed = 8, max_iter = 3000, tol = 1e-7)
  recorder <- refit_recorder()
  testthat::local_mocked_bindings(multilpa = recorder$mock)
  result <- bootstrap_lrt(small, large, data, iter = 2, n_starts = 2,
                          max_iter = 3000, tol = 1e-7, seed = 42)
  calls <- recorder$calls()
  expect_identical(length(calls), 4L)
  expect_true(all(vapply(calls, function(arguments)
    length(arguments$fixed) == 0L && is.null(arguments$start), logical(1))))
  expect_identical(result$fixed, character(0))
})

test_that("bootstrap refits retain grand centering", {
  set.seed(91)
  data <- data.frame(g = rep(seq_len(20L), each = 6L),
                     y = c(stats::rnorm(60L, 28), stats::rnorm(60L, 32)))
  small <- multilpa(data, "y", "g", 1L, 1L, centering = "grand",
                    n_starts = 1L, seed = 1L)
  large <- multilpa(data, "y", "g", 2L, 1L, centering = "grand",
                    n_starts = 1L, seed = 1L)
  recorder <- refit_recorder()
  testthat::local_mocked_bindings(multilpa = recorder$mock)
  quietly(bootstrap_lrt(small, large, data, iter = 2L, n_starts = 1L,
                        seed = 2L))
  calls <- recorder$calls()
  expect_identical(length(calls), 4L)
  expect_true(all(vapply(calls, function(arguments)
    identical(arguments$centering, "grand"), logical(1))))
})

test_that("bootstrap refuses incomparable structures and person centering", {
  set.seed(92)
  data <- data.frame(g = rep(seq_len(12L), each = 6L),
                     a = stats::rnorm(72L), b = stats::rnorm(72L))
  fit <- function(n_group_classes, centering = "none", shape = "equal") {
    quietly(multilpa(data, c("a", "b"), "g", 2L, n_group_classes,
                      volume = "varying", shape = shape,
                      orientation = "axis", centering = centering,
                      n_starts = 1L, max_iter = 0L, seed = 1L))
  }
  expect_error(bootstrap_lrt(fit(1L), fit(2L, shape = "varying"),
                             data, iter = 2L),
               class = "multilpa_incomparable_models")
  expect_error(bootstrap_lrt(fit(1L), fit(2L, centering = "grand"),
                             data, iter = 2L),
               class = "multilpa_incomparable_models")
  expect_error(bootstrap_lrt(fit(1L, centering = "person"),
                             fit(2L, centering = "person"),
                             data, iter = 2L),
               class = "multilpa_unsupported_bootstrap")
})

test_that("the held constraint is named on the public surface", {
  skip_on_cran()
  models <- constrained_pair()
  result <- bootstrap_lrt(models$null_model, models$alternative_model,
                          models$data, iter = 2, n_starts = 2, seed = 45,
                          tol = 1e-8)
  test <- get_results(result, "test")
  expect_identical(test$fixed, "means, variances")
  expect_output(print(result), "conditional on it")
})

test_that("a constrained pair that is not nested is refused", {
  skip_on_cran()
  models <- constrained_pair()
  # Same blocks held, but at a measurement solution the null does not share.
  other <- multilpa(models$data, c("a", "b"), "g", 2, 1, n_starts = 1,
                    seed = 11, tol = 1e-10)
  shifted <- starting_values(other, what = "measurement")
  shifted$means <- shifted$means + 0.5
  mismatched <- multilpa(models$data, c("a", "b"), "g", 2, 2, n_starts = 2,
                         start = shifted, fixed = "measurement", seed = 5,
                         tol = 1e-10)
  expect_error(
    bootstrap_lrt(models$null_model, mismatched, models$data, iter = 2,
                  seed = 1),
    class = "multilpa_bad_nesting")
})

test_that("a held measurement is refused across different profile counts", {
  skip_on_cran()
  models <- constrained_pair()
  three <- multilpa(models$data, c("a", "b"), "g", 3, 1, n_starts = 1, seed = 5,
                    tol = 1e-10)
  larger <- multilpa(models$data, c("a", "b"), "g", 3, 1, n_starts = 1,
                     start = starting_values(three, what = "measurement"),
                     fixed = "measurement", seed = 5, tol = 1e-10)
  expect_error(
    bootstrap_lrt(models$null_model, larger, models$data, iter = 2, seed = 1),
    class = "multilpa_bad_nesting")
})

test_that("a constraint on only one of the two models is refused", {
  skip_on_cran()
  models <- constrained_pair()
  free_alternative <- multilpa(models$data, c("a", "b"), "g", 2, 2,
                               n_starts = 3, seed = 5, tol = 1e-10)
  expect_error(
    bootstrap_lrt(models$null_model, free_alternative, models$data, iter = 2,
                  seed = 1),
    class = "multilpa_bad_nesting")
  # And the other way round: an unconstrained null against a constrained
  # alternative is equally not a nested pair.
  free_null <- multilpa(models$data, c("a", "b"), "g", 2, 1, n_starts = 3,
                        seed = 5, tol = 1e-10)
  expect_error(
    bootstrap_lrt(free_null, models$alternative_model, models$data, iter = 2,
                  seed = 1),
    class = "multilpa_bad_nesting")
})

test_that("models holding different blocks are refused", {
  skip_on_cran()
  models <- constrained_pair()
  means_only <- multilpa(models$data, c("a", "b"), "g", 2, 2, n_starts = 3,
                         start = models$start, fixed = "means", seed = 5,
                         tol = 1e-10)
  expect_error(
    bootstrap_lrt(models$null_model, means_only, models$data, iter = 2,
                  seed = 1),
    class = "multilpa_bad_nesting")
})
