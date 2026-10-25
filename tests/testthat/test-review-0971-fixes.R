# Reporting fixes found while reading a 49-candidate multilevel LCA grid.

fixes_0971_grid <- function() {
  # Two distinct indicator rows, so three classes cannot be fitted.
  data <- data.frame(g = rep(seq_len(10L), each = 3L), y = rep(c(0L, 1L), 15L))
  suppressWarnings(enumerate_classes(data, "y", "g", n_profiles = 2:3,
                                     n_group_classes = 1L, categorical = "y",
                                     n_starts = 1L, seed = 1L))
}

test_that("a candidate that never fitted is reported as failed, not unconverged", {
  grid <- fixes_0971_grid()
  table <- as.data.frame(grid)
  expect_true(is.na(table$log_likelihood[2L]))
  expect_false(is.na(table$error[2L]))
  printed <- capture.output(print(grid))
  expect_true(any(grepl("1 candidate(s) failed to fit", printed, fixed = TRUE)))
  expect_false(any(grepl("did not converge", printed, fixed = TRUE)))
})

test_that("an LCA grid prints no empty covariance-structure column", {
  grid <- fixes_0971_grid()
  # The returned table keeps the column; only the printed view leaves it out.
  expect_true(all(is.na(as.data.frame(grid)$model)))
  header <- capture.output(print(grid))[2L]
  expect_false(grepl("\\bmodel\\b", header))
  summary_text <- capture.output(print(summary(grid)))
  expect_false(any(grepl("\\bmodel\\b", summary_text)))
  expect_true("model" %in% names(get_results(grid, "best")))
})

test_that("perfect classification has infinite odds, an undefined one NA", {
  odds <- .multilpa_odds_correct(c(1, 1 - 1e-12, 0.9, 0.9, NA),
                                 c(0.5, 0.5, 0.5, 1, 0.5))
  expect_identical(odds[1:2], c(Inf, Inf))
  expect_equal(odds[3L], 9)
  expect_true(all(is.na(odds[4:5])))
  # Invariant: below the rounding band the odds rise with the posterior.
  rising <- .multilpa_odds_correct(c(0.6, 0.8, 0.99, 0.999999), rep(0.3, 4))
  expect_true(all(diff(rising) > 0))
  probabilities <- cbind(c(1, 1, 0), c(0, 0, 1))
  object <- structure(list(subject_posteriors = probabilities), class = "multilpa")
  table <- get_results(object, "classification", level = "individuals")
  expect_identical(table$odds_correct_classification, c(Inf, Inf))
})

test_that("a k-means start that stops at its iteration cap is a recorded failed start", {
  set.seed(3)
  data <- data.frame(y = c(stats::rnorm(30), stats::rnorm(30, 4)))
  original <- stats::kmeans
  local_mocked_bindings(kmeans = function(...) {
    fit <- original(...)
    warning("did not converge in 100 iterations", call. = FALSE)
    fit$ifault <- 2L
    fit
  }, .package = "stats")
  warnings <- character()
  fit <- withCallingHandlers(
    multilpa(data, "y", NULL, 2L, n_starts = 3L, seed = 1L),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
  expect_s3_class(fit, "multilpa")
  # Only the summary of failed starts reaches the caller, never R's raw
  # k-means warning worded as though the model had not converged.
  expect_false(any(startsWith(warnings, "did not converge in")))
  expect_true(any(grepl("starts failed", warnings, fixed = TRUE)))
  starts <- get_results(fit, "starts")
  reasons <- starts$error[!is.na(starts$error)]
  expect_gt(length(reasons), 0L)
  expect_true(all(grepl("k-means clustering that draws this start", reasons,
                        fixed = TRUE)))
})

test_that("progress durations are given in words", {
  expect_identical(.latents_format_duration(30), "under a minute")
  expect_identical(.latents_format_duration(42 * 60), "42 min")
  expect_identical(.latents_format_duration(72 * 60), "1 h 12 min")
  expect_identical(.latents_count(1L, "worker"), "1 worker")
  expect_identical(.latents_count(12L, "worker"), "12 workers")
  expect_error(.latents_format_duration(-1))
})

test_that("a one-core grid reports no progress unless asked", {
  data <- subset(course_engagement, student <= 20)
  expect_no_message(enumerate_classes(data, c("browse", "lectures"), "student",
                                      n_profiles = 1:2, n_group_classes = 1L,
                                      model = "VVI", n_starts = 1, seed = 1))
  progress <- character()
  withCallingHandlers(
    enumerate_classes(data, c("browse", "lectures"), "student",
                      n_profiles = 1:2, n_group_classes = 1L, model = "VVI",
                      n_starts = 1, seed = 1, verbose = TRUE),
    latents_progress = function(condition) {
      progress <<- c(progress, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  expect_length(progress, 3L)
  expect_match(progress[1L], "^Preparing 2 candidates on 1 worker\\.")
  expect_match(progress[2L], "^Fitting 2 starts on 1 worker\\.")
  expect_match(progress[3L], "^Done in ")
  expect_error(enumerate_classes(data, c("browse", "lectures"), "student",
                                 n_profiles = 1:2, verbose = NA))
})
