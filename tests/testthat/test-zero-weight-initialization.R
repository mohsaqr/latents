# Excluded independent units must not change starting states or the internal
# coordinate origin. Relative positive weights retain the established starts.

test_that("zero-weight Gaussian outliers cannot consume a starting profile", {
  set.seed(803)
  active <- data.frame(y = c(rnorm(50, -3, 0.5), rnorm(50, 3, 0.5)), w = 1)
  fit_at <- function(value) {
    quietly(lpa(rbind(active, data.frame(y = value, w = 0)), "y", 2,
                weights = "w", n_starts = 3, seed = 1))
  }
  reference <- fit_at(0)
  vapply(c(1e4, 1e20), function(value) {
    fitted <- fit_at(value)
    expect_equal(fitted$log_likelihood, reference$log_likelihood, tolerance = 1e-12)
    expect_identical(fitted$means, reference$means)
    expect_identical(fitted$variances, reference$variances)
    expect_equal(fitted$subject_posteriors[seq_len(nrow(active)), ],
                 reference$subject_posteriors[seq_len(nrow(active)), ])
    expect_equal(nrow(fitted$subject_posteriors), nrow(active) + 1L)
    expect_equal(unname(rowSums(fitted$subject_posteriors)), rep(1, nrow(active) + 1L))
    TRUE
  }, logical(1))
  expect_gt(diff(range(reference$means)), 5)
})

test_that("all start kinds use positive support for mixed and missing indicators", {
  set.seed(107)
  frame <- data.frame(a = rnorm(60), b = rnorm(60),
                      category = factor(rep(c("a", "b", "c"), 20)),
                      ordinal = ordered(rep(1:3, 20)), count = rpois(60, 3))
  frame$a[c(2, 9)] <- NA
  frame$b[c(3, 7)] <- NA
  frame$category[c(1, 4)] <- NA
  frame$ordinal[5] <- NA
  frame$count[6] <- NA
  x <- as.matrix(frame[c("a", "b")])
  codes <- matrix(as.integer(frame$category), ncol = 1)
  extra <- .latents_prepare_extra(frame, "ordinal", "count", "fiml")
  extended_x <- rbind(c(1e20, -1e20), x, c(NA, 1e12))
  extended_codes <- rbind(3L, codes, NA_integer_)
  extended_frame <- frame[c(1, seq_len(nrow(frame)), 2), ]
  extended_frame$count[c(1, nrow(extended_frame))] <- c(1e12, 0)
  extended_extra <- .latents_prepare_extra(extended_frame, "ordinal", "count", "fiml")
  groups <- rep(seq_len(15), each = 4)
  extended_groups <- c(1L, groups + 1L, 17L)
  vapply(c("ward", "kmeans", "partition", "centres"), function(kind) {
    args <- list(n_profiles = 2, n_types = 2, variance_model = "varying",
                 min_variance = 1e-6, start_index = 3,
                 covariance_model = "full", n_categories = 3,
                 hierarchical = identical(kind, "ward"),
                 kind = if (identical(kind, "ward")) "kmeans" else kind)
    set.seed(91)
    reference <- do.call(.multilpa_initialize,
                         c(list(x = x, group_index = groups, codes = codes,
                                extra = extra), args))
    reference_rng <- .Random.seed
    set.seed(91)
    fitted <- do.call(.multilpa_initialize,
                      c(list(x = extended_x, group_index = extended_groups,
                             codes = extended_codes, extra = extended_extra,
                             row_weights = c(0, rep(c(1, 2), 30), 0)), args))
    expect_identical(fitted, reference)
    expect_identical(.Random.seed, reference_rng)
    TRUE
  }, logical(1))
})

test_that("all-positive weights leave initialization and its RNG unchanged", {
  set.seed(118)
  x <- matrix(rnorm(120), 60, 2)
  initialize <- function(weights) {
    .multilpa_initialize(x, rep(seq_len(15), each = 4), 2, 2, "varying",
                         1e-6, 3, row_weights = weights)
  }
  set.seed(76)
  plain <- initialize(NULL)
  plain_rng <- .Random.seed
  set.seed(76)
  expect_identical(initialize(rep(c(0.5, 2), 30)), plain)
  expect_identical(.Random.seed, plain_rng)
  expect_identical(.multilpa_initialization_center(x), colMeans(x))
  expect_identical(.multilpa_initialization_center(x, rep(1, nrow(x))), colMeans(x))
})

test_that("excluded groups retain the full observation and group layout", {
  set.seed(223)
  active <- data.frame(id = rep(seq_len(20), each = 4),
                        a = rnorm(80, rep(c(-3, 3), 40), 0.6),
                        b = rnorm(80, rep(c(-2, 2), 40), 0.6), w = 1)
  fit_at <- function(value) {
    data <- rbind(data.frame(id = 0, a = value, b = -value, w = 0), active)
    quietly(multilpa(data, c("a", "b"), "id", 2, 2,
                     weights = "w", n_starts = 3, max_iter = 20, seed = 7))
  }
  reference <- fit_at(0)
  fitted <- fit_at(1e20)
  expect_identical(fitted$means, reference$means)
  expect_identical(fitted$variances, reference$variances)
  expect_equal(fitted$log_likelihood, reference$log_likelihood, tolerance = 1e-12)
  expect_equal(nrow(fitted$subject_posteriors), 81L)
  expect_equal(nrow(fitted$group_posteriors), 21L)
  expect_equal(unname(rowSums(fitted$group_posteriors)), rep(1, 21))
})

test_that("covariate and transition starts exclude zero-weight indicator outliers", {
  set.seed(304)
  active <- data.frame(id = rep(seq_len(30), each = 2), time = rep(1:2, 30),
                        y = rnorm(60, rep(c(-3, 3), each = 30), 0.6),
                        z = rnorm(60), w = 1)
  extend <- function(value) {
    rbind(active, data.frame(id = 31, time = 1:2,
                             y = c(value, -value), z = c(0, 1), w = 0))
  }
  cov_fit <- function(value) {
    quietly(lpa(extend(value), "y", 2, profile_covariates = "z", weights = "w",
                n_starts = 3, max_iter = 0, seed = 5))
  }
  reference <- cov_fit(0)
  fitted <- cov_fit(1e20)
  expect_identical(fitted$means, reference$means)
  expect_identical(fitted$variances, reference$variances)
  expect_equal(fitted$log_likelihood, reference$log_likelihood, tolerance = 1e-12)
  transition_fit <- function(value) {
    quietly(lta(extend(value), "y", "id", 2, time = "time",
                                n_group_classes = 1, weights = "w",
                                n_starts = 3, max_iter = 0, seed = 5))
  }
  reference <- transition_fit(0)
  fitted <- transition_fit(1e20)
  expect_identical(fitted$means, reference$means)
  expect_identical(fitted$variances, reference$variances)
  expect_equal(fitted$log_likelihood, reference$log_likelihood, tolerance = 1e-12)
  expect_equal(nrow(fitted$group_posteriors), 31L)
})

test_that("insufficient positive-weight indicator support fails clearly", {
  data <- data.frame(y = c(NA_real_, NA_real_, 2, 3), z = 0:3, w = c(1, 1, 0, 0))
  expect_error(quietly(lpa(data, c("y", "z"), 1, weights = "w", missing = "fiml")),
                "observed row with positive", class = "latents_bad_data")
  # The rows that carry weight hold one value, which is refused as the same
  # data without its zero-weight row is (ML-2, review 2026-10-26); it used to
  # reach the starts and fail there.
  data <- data.frame(y = c(1, 1, 4), w = c(1, 1, 0))
  expect_error(quietly(lpa(data, "y", 2, weights = "w", n_starts = 1)),
               class = "latents_bad_data")
})
