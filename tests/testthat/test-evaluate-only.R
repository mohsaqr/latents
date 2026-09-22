# Warm starts for categorical models, zero-iteration evaluation, and the
# information-carrying sample size used by the individual-level BIC.

categorical_frame <- function(n = 240L, seed = 11L) {
  set.seed(seed)
  profile <- rep(seq_len(2L), length.out = n)
  draw <- function(high) {
    ifelse(stats::runif(n) < ifelse(profile == 1L, high, 1 - high), "yes", "no")
  }
  data.frame(school = rep(seq_len(12L), each = n / 12L),
             a = draw(0.85), b = draw(0.8), c = draw(0.75),
             stringsAsFactors = FALSE)
}

test_that("starting_values round-trips a categorical fit exactly", {
  frame <- categorical_frame()
  items <- c("a", "b", "c")
  fit <- multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 2,
                  categorical = items, n_starts = 4, seed = 1)

  start <- starting_values(fit)
  expect_s3_class(start, "multilpa_start")
  expect_true("response_probabilities" %in% names(start))
  expect_length(start$response_probabilities, length(items))
  expect_equal(dim(start$response_probabilities[[1L]]), c(2L, 2L))

  evaluated <- multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 2,
                        categorical = items, n_starts = 1, start = start,
                        max_iter = 0)
  expect_identical(evaluated$iterations, 0L)
  expect_equal(as.numeric(logLik(evaluated)), as.numeric(logLik(fit)),
               tolerance = 1e-12)
})

test_that("max_iter = 0 scores externally supplied parameters without moving", {
  frame <- categorical_frame()
  items <- c("a", "b", "c")
  fit <- multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
                  categorical = items, n_starts = 4, seed = 1)
  start <- starting_values(fit)
  # Perturb one item's response probabilities: the evaluated likelihood must
  # move, and must fall, because the fit was at a maximum.
  start$response_probabilities[[1L]] <- matrix(c(0.5, 0.5, 0.5, 0.5), 2L, 2L)
  evaluated <- multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
                        categorical = items, n_starts = 1, start = start,
                        max_iter = 0)
  expect_lt(as.numeric(logLik(evaluated)), as.numeric(logLik(fit)))
  expect_identical(evaluated$iterations, 0L)
})

test_that("max_iter = 0 does not warn about convergence it was not asked to reach", {
  frame <- categorical_frame()
  items <- c("a", "b", "c")
  fit <- multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
                  categorical = items, n_starts = 4, seed = 1)
  expect_no_warning(
    multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
             categorical = items, n_starts = 1, start = starting_values(fit),
             max_iter = 0))
})

test_that("a malformed categorical start is refused", {
  frame <- categorical_frame()
  items <- c("a", "b", "c")
  fit <- multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
                  categorical = items, n_starts = 4, seed = 1)
  start <- starting_values(fit)

  dropped <- start
  dropped$response_probabilities <- NULL
  expect_error(
    multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
             categorical = items, n_starts = 1, start = dropped, max_iter = 0),
    "response_probabilities")

  unnormalized <- start
  unnormalized$response_probabilities[[2L]] <- matrix(c(0.9, 0.9, 0.9, 0.9), 2L, 2L)
  expect_error(
    multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
             categorical = items, n_starts = 1, start = unnormalized, max_iter = 0),
    "sum to one")

  wrong_shape <- start
  wrong_shape$response_probabilities[[3L]] <- matrix(0.25, 2L, 4L)
  expect_error(
    multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
             categorical = items, n_starts = 1, start = wrong_shape, max_iter = 0),
    "2 x 2")
})

test_that("max_iter must be a nonnegative integer", {
  frame <- categorical_frame()
  items <- c("a", "b", "c")
  expect_error(
    multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
             categorical = items, max_iter = -1),
    "nonnegative integer")
  expect_error(
    multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
             categorical = items, max_iter = 2.5),
    "nonnegative integer")
})

test_that("rows with no observed indicator are excluded from the BIC sample size", {
  frame <- categorical_frame()
  items <- c("a", "b", "c")
  fit <- multilpa(frame, items, "school", n_profiles = 2, n_group_classes = 1,
                  categorical = items, n_starts = 4, seed = 1)
  expect_identical(fit$n_informative, fit$n_observations)

  blanked <- frame
  blanked[c(5L, 60L, 120L), items] <- NA
  sparse <- multilpa(blanked, items, "school", n_profiles = 2, n_group_classes = 1,
                     categorical = items, n_starts = 4, seed = 1, missing = "fiml")
  # Reported, not warned about, so that resampling verbs do not emit one warning
  # per replicate.
  expect_no_warning(multilpa(blanked, items, "school", n_profiles = 2,
                             n_group_classes = 1, categorical = items,
                             n_starts = 1, seed = 1, missing = "fiml"))
  expect_output(print(sparse), "carry no observed indicator")
  expect_identical(get_results(summary(sparse), "model")$n_informative,
                   sparse$n_informative)
  expect_identical(sparse$n_observations, nrow(blanked))
  expect_identical(sparse$n_informative, nrow(blanked) - 3L)
  # An uninformative row changes neither the likelihood nor the penalty basis.
  expect_equal(sparse$bic_individual,
               -2 * as.numeric(logLik(sparse)) +
                 log(sparse$n_informative) * sparse$n_parameters)
})
