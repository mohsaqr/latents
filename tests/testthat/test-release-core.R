# Pre-release review 2026-10-27, cluster "core": CORE-2, CORE-3, CORE-4,
# CORE-5 and CORE-7.

# A weighted fit and the same fit with one appended row of weight zero.
release_core_ordinal_data <- function() {
  .latents_local_seed(4)
  n <- 120
  class <- stats::rbinom(n, 1, 0.5)
  data <- data.frame(y1 = stats::rnorm(n, class * 2),
                     o = pmin(pmax(round(stats::rnorm(n, 3 + class)), 1), 5),
                     q = stats::rpois(n, 2 + 3 * class), w = 1)
  list(data = data,
       padded = rbind(data, data.frame(y1 = 0, o = 9, q = 50, w = 0)))
}

test_that("an ordinal value seen only on a zero-weight row adds no category (CORE-2)", {
  sets <- release_core_ordinal_data()
  fit <- function(data) {
    lpa(data, c("y1", "o", "q"), 2, ordinal = "o", count = "q",
        n_starts = 2, seed = 1, weights = "w")
  }
  original <- fit(sets$data)
  # No `latents_unused_category` warning: the value is not a category.
  padded <- expect_no_warning(fit(sets$padded))
  expect_identical(padded$n_parameters, original$n_parameters)
  expect_equal(padded$bic, original$bic, tolerance = 1e-10)
  expect_equal(padded$log_likelihood, original$log_likelihood, tolerance = 1e-10)
  expect_identical(padded$extra_data$ordinal_levels, original$extra_data$ordinal_levels)
  expect_equal(padded$ordinal_intercepts, original$ordinal_intercepts, tolerance = 1e-6)
  # The zero-weight row's value has no fitted category.
  expect_true(is.na(padded$extra_data$ordinal[nrow(sets$padded), "o"]))
  # Refitting the rebuilt data gives the same model.
  rebuilt <- get_results(padded, "data")
  expect_identical(nrow(rebuilt), nrow(sets$padded))
  expect_equal(fit(rebuilt)$log_likelihood, original$log_likelihood, tolerance = 1e-10)
})

test_that("zero-weight rows cannot rescue a constant ordinal or all-zero count (CORE-2)", {
  .latents_local_seed(1)
  count_data <- data.frame(y1 = c(stats::rnorm(30, 3), stats::rnorm(30)), q = 0L, w = 1)
  count_padded <- rbind(count_data, data.frame(y1 = 0, q = 4L, w = 0))
  fit_count <- function(data) {
    lpa(data, c("y1", "q"), 2, count = "q", weights = "w", n_starts = 2, seed = 1)
  }
  expect_error(fit_count(count_data), class = "latents_bad_data")
  expect_error(fit_count(count_padded), class = "latents_bad_data")
  ordinal_data <- data.frame(y1 = c(stats::rnorm(30), stats::rnorm(30, 3)), o = 2L, w = 1)
  ordinal_padded <- rbind(ordinal_data, data.frame(y1 = 0, o = 1L, w = 0))
  fit_ordinal <- function(data) {
    lpa(data, c("y1", "o"), 2, ordinal = "o", weights = "w", n_starts = 2, seed = 1)
  }
  expect_error(fit_ordinal(ordinal_data), class = "latents_bad_data")
  expect_error(fit_ordinal(ordinal_padded), class = "latents_bad_data")
})

test_that("a missing ordinal or count value on a zero-weight row is absent from the fit (CORE-2)", {
  sets <- release_core_ordinal_data()
  data <- sets$data
  padded <- rbind(data, data.frame(y1 = 0, o = NA, q = NA, w = 0))
  fit <- function(d) {
    lpa(d, c("y1", "o", "q"), 2, ordinal = "o", count = "q",
        n_starts = 2, seed = 1, weights = "w")
  }
  original <- fit(data)
  # Under `missing = "error"`, as a missing categorical value on such a row is.
  actual <- expect_no_warning(fit(padded))
  expect_identical(actual$n_parameters, original$n_parameters)
  expect_equal(actual$log_likelihood, original$log_likelihood, tolerance = 1e-10)
  # A missing value on a row of positive weight is still refused, for either type.
  count_missing <- data
  count_missing$q[1] <- NA
  expect_error(fit(count_missing), class = "latents_bad_data")
  ordinal_missing <- data
  ordinal_missing$o[1] <- NA
  expect_error(fit(ordinal_missing), class = "latents_bad_data")
})

test_that("an unweighted ordinal fit keeps every observed category (CORE-2)", {
  data <- release_core_ordinal_data()$data
  data$o[1] <- 9
  fit <- lpa(data, c("y1", "o"), 2, ordinal = "o", n_starts = 1, seed = 1)
  expect_identical(fit$extra_data$ordinal_levels$o, c("1", "2", "3", "4", "5", "9"))
})

test_that("a pair item's zero-weight-only category does not enlarge the joint grid (CORE-3)", {
  items <- c("time_with_friends", "on_social_media", "sports", "walking")
  data <- student_esm[, items]
  data[] <- lapply(data, as.character)
  data$w <- 1
  padded <- rbind(data, data.frame(time_with_friends = "maybe", on_social_media = "no",
                                   sports = "no", walking = "no", w = 0))
  fit <- function(d) {
    lca(d, items, 2, local_dependence = "time_with_friends:on_social_media",
        weights = "w", n_starts = 1, seed = 1)
  }
  original <- fit(data)
  actual <- fit(padded)
  expect_identical(actual$n_parameters, original$n_parameters)
  expect_equal(actual$bic, original$bic, tolerance = 1e-10)
  expect_equal(actual$log_likelihood, original$log_likelihood, tolerance = 1e-10)
  expect_identical(actual$categorical_levels, original$categorical_levels)
  expect_length(actual$categorical_levels[["time_with_friends:on_social_media"]], 4L)
  # Every combination of the categories of positive weight stays in the grid,
  # observed or not, as before.
  thin <- data.frame(a = c("a", "a", "b", "c"), b = c("x", "y", "x", "x"),
                     w = c(1, 1, 1, 0))
  joint <- suppressWarnings(lca(thin, c("a", "b"), 1, local_dependence = "a:b",
                                weights = "w", n_starts = 1, seed = 1))
  expect_identical(joint$categorical_levels[["a:b"]], c("a:x", "a:y", "b:x", "b:y"))
})

test_that("residuals and diagnostics work after a zero-weight-only category (CORE-4)", {
  items <- c("time_with_friends", "on_social_media", "sports", "walking")
  data <- student_esm[, items]
  data[] <- lapply(data, as.character)
  data$w <- 1
  padded <- rbind(data, data.frame(time_with_friends = "no", on_social_media = "no",
                                   sports = "maybe", walking = "no", w = 0))
  fit <- function(d) lca(d, items, 2, weights = "w", n_starts = 1, seed = 1)
  original <- fit(data)
  actual <- fit(padded)
  # The excluded row carries weight zero, so every residual table is unchanged.
  residuals <- get_results(actual, "residuals")
  expect_s3_class(residuals, "data.frame")
  expect_equal(residuals, get_results(original, "residuals"), tolerance = 1e-8)
  expect_s3_class(diagnostics(actual), "multilpa_diagnostics")
  # A category the fit never saw on a row of positive weight is refused, by class.
  foreign <- get_results(actual, "data")
  foreign$sports[1] <- "never"
  expect_error(.multilpa_discrete_items(actual, foreign),
               class = "latents_bad_inference_data")
})

test_that("an all-categorical model with negative degrees of freedom warns (CORE-5)", {
  data <- student_esm[, c("sports", "walking", "reading")]
  # Two binary items, two classes: 5 parameters for 3 free cells.
  expect_warning(lca(data, c("sports", "walking"), 2, n_starts = 1, seed = 1),
                 class = "latents_unidentified")
  # Three binary items, two classes: 7 parameters for 7 free cells, df = 0.
  expect_no_warning(lca(data, c("sports", "walking", "reading"), 2,
                        n_starts = 1, seed = 1))
  expect_no_warning(lca(data, c("sports", "walking"), 1, n_starts = 1, seed = 1))
  # The estimates themselves are unchanged by the warning.
  fit <- suppressWarnings(lca(data, c("sports", "walking"), 2, n_starts = 1, seed = 1))
  expect_equal(fit$n_parameters, 5)
  # A model with a continuous indicator is not judged by counting cells.
  mixed <- data.frame(y = stats::rnorm(nrow(data)), sports = data$sports)
  fit_mixed <- withCallingHandlers(
    lpa(mixed, c("y", "sports"), 2, categorical = "sports", n_starts = 1, seed = 1),
    warning = function(w) {
      expect_false(inherits(w, "latents_unidentified"))
      invokeRestart("muffleWarning")
    })
  expect_s3_class(fit_mixed, "multilpa")
})

test_that("invalid enumerated arguments raise latents_bad_argument (CORE-7)", {
  data <- course_engagement
  vars <- c("browse", "lectures")
  fit <- function(...) lpa(data, vars, 2, n_starts = 1, seed = 1, ...)
  expect_error(fit(select_start = "x"), class = "latents_bad_argument")
  expect_error(fit(boundary = "x"), class = "latents_bad_argument")
  expect_error(fit(acceleration = "x"), class = "latents_bad_argument")
  expect_error(fit(centering = "x"), class = "latents_bad_argument")
  expect_error(fit(missing = "listwise"), class = "latents_bad_argument")
  expect_error(fit(missing = NA_character_), class = "latents_bad_argument")
  expect_error(fit(missing = c("error", "x")), class = "latents_bad_argument")
  expect_error(fit(variance_model = 1), class = "latents_bad_argument")
  # The message names the argument rather than `'arg'`.
  expect_error(fit(boundary = "x"), "`boundary`")
  # An unambiguous prefix is still accepted, as match.arg() accepts it.
  expect_identical(fit(missing = "fim")$missing, "fiml")
  expect_identical(fit(missing = c("error", "fiml"))$missing, "error")
})

test_that("verbs work on a fit whose zero-weight row holds an excluded or missing ordinal or count value (CORE-2)", {
  sets <- release_core_ordinal_data()
  fit <- function(d) {
    lpa(d, c("y1", "o", "q"), 2, ordinal = "o", count = "q",
        n_starts = 2, seed = 1, weights = "w")
  }
  original <- fit(sets$data)
  zero_rows <- list(excluded = data.frame(y1 = 0, o = 9, q = 3, w = 0),
                    missing = data.frame(y1 = 0, o = NA, q = NA, w = 0))
  lapply(zero_rows, function(row) {
    padded <- rbind(sets$data, row)
    actual <- fit(padded)
    predicted <- predict(actual)
    expect_s3_class(predicted, "data.frame")
    expect_identical(nrow(predicted), nrow(padded))
    expect_equal(unname(vcov(actual)), unname(vcov(original)), tolerance = 1e-6)
    expect_equal(unname(vcov(actual, data = padded)), unname(vcov(original)),
                 tolerance = 1e-6)
    expect_s3_class(get_results(actual, "residuals", data = padded), "data.frame")
  })
  # New data are still checked on every row.
  unseen <- sets$data
  unseen$o[1] <- 9
  expect_error(predict(original, newdata = unseen), class = "latents_bad_data")
})

test_that("a covariate model's internal base fit does not warn latents_unidentified (CORE-5)", {
  .latents_local_seed(2)
  n <- 300
  z <- stats::rnorm(n)
  class <- stats::rbinom(n, 1, stats::plogis(1.5 * z))
  data <- data.frame(a = stats::rbinom(n, 1, ifelse(class == 1, 0.85, 0.15)),
                     b = stats::rbinom(n, 1, ifelse(class == 1, 0.8, 0.2)), x = z)
  # The unconditional model is not identified by counting ...
  expect_warning(lca(data, c("a", "b"), 2, n_starts = 1, seed = 1),
                 class = "latents_unidentified")
  # ... but the covariate model is not that model, and is not judged by it.
  fit <- withCallingHandlers(
    lca(data, c("a", "b"), 2, profile_covariates = "x", n_starts = 1, seed = 1),
    warning = function(w) {
      expect_false(inherits(w, "latents_unidentified"))
      invokeRestart("muffleWarning")
    })
  expect_s3_class(fit, "multilpa_covariates")
})

test_that("covariate and transition fits follow the zero-weight rule for ordinal values (CORE-2, H4)", {
  .latents_local_seed(4)
  n <- 120
  class <- stats::rbinom(n, 1, 0.5)
  data <- data.frame(y1 = stats::rnorm(n, class * 2),
                     o = pmin(pmax(round(stats::rnorm(n, 3 + class)), 1), 5),
                     z = stats::rnorm(n), w = 1)
  covariate_fit <- function(d) {
    lpa(d, c("y1", "o"), 2, ordinal = "o", profile_covariates = "z",
        n_starts = 2, seed = 1, weights = "w")
  }
  original <- covariate_fit(data)
  zero_rows <- list(excluded = data.frame(y1 = 0, o = 9, z = 0, w = 0),
                    missing = data.frame(y1 = 0, o = NA, z = 0, w = 0))
  lapply(zero_rows, function(row) {
    padded <- rbind(data, row)
    actual <- covariate_fit(padded)
    expect_identical(actual$n_parameters, original$n_parameters)
    expect_identical(actual$extra_data$ordinal_levels, original$extra_data$ordinal_levels)
    expect_equal(actual$log_likelihood, original$log_likelihood, tolerance = 1e-10)
    expect_identical(dim(vcov(actual)), dim(vcov(original)))
  })

  .latents_local_seed(7)
  long <- data.frame(id = rep(seq_len(40), each = 3), time = rep(1:3, 40))
  state <- stats::rbinom(nrow(long), 1, 0.5)
  long$y <- stats::rnorm(nrow(long), state * 2)
  long$o <- pmin(pmax(round(stats::rnorm(nrow(long), 2 + state)), 1), 4)
  long$w <- 1
  transition_fit <- function(d) {
    lta(d, c("y", "o"), "id", 2, time = "time", ordinal = "o", weights = "w",
        n_starts = 2, seed = 1)
  }
  reference <- transition_fit(long)
  lapply(list(9, NA), function(value) {
    padded <- rbind(long, data.frame(id = 41, time = 1:3, y = 0, o = value, w = 0))
    actual <- transition_fit(padded)
    expect_identical(actual$n_parameters, reference$n_parameters)
    expect_identical(actual$extra_data$ordinal_levels, reference$extra_data$ordinal_levels)
    expect_equal(actual$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
    # The data check of inference re-encodes the rows as the fit did.
    expect_no_error(parameter_inference(actual, data = padded))
  })
})
