# Numeric category summaries use fitted response probabilities, even with FIML.
.categorical_summary_fixture <- function(n_profiles = 2L, missing = FALSE,
                                         mixed = FALSE, covariates = FALSE) {
  set.seed(123)
  dat <- data.frame(
    group = rep(seq_len(12L), each = 10L),
    a = sample(1:5, 120L, replace = TRUE),
    b = sample(1:5, 120L, replace = TRUE),
    score = stats::rnorm(120L), z = stats::rnorm(120L))
  if (missing) dat$a[seq.int(1L, 120L, by = 7L)] <- NA_integer_
  multilpa(dat, if (mixed) c("a", "score", "b") else c("a", "b"),
    "group", n_profiles, 1L, categorical = c("a", "b"),
    profile_covariates = if (covariates) "z" else character(),
    missing = if (missing) "fiml" else "error", n_starts = 1L, seed = 1L)
}

test_that("categorical profile summaries draw probability-weighted scores", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture()
  fit$response_probabilities$a[,] <- rbind(
    c(.10, .40, .30, .15, .05), c(.05, .10, .15, .30, .40))
  fit$response_probabilities$b[,] <- rbind(
    c(.10, .15, .10, .25, .40), c(.25, .25, 0, .25, .25))
  expected <- list(mean = rbind(c(2.65, 3.70), c(3.90, 3.00)),
                   median = rbind(c(2, 4), c(4, 2)),
                   mode = rbind(c(2, 5), c(5, 1)))
  invisible(lapply(names(expected), function(statistic) {
    plotted <- expect_plot(plot(fit, what = "profiles", statistic = statistic,
                                intervals = FALSE))
    drawn <- plotted$data
    expect_equal(drawn$mean,
                 expected[[statistic]][cbind(drawn$profile, drawn$position)])
    expect_match(plotted$labels$y, statistic, ignore.case = TRUE)
    expect_true(all(drawn$mean >= 1 & drawn$mean <= 5))
    expect_false(anyNA(drawn$mean))
  }))
  expect_equal(plot(fit, intervals = FALSE)$data,
               plot(fit, statistic = "mean", intervals = FALSE)$data)
})

test_that("numeric ordering determines medians and modal ties", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture()
  scores <- c("5", "1", "3", "2", "4")
  fit$categorical_levels$a <- scores
  colnames(fit$response_probabilities$a) <- scores
  fit$response_probabilities$a[,] <- rbind(
    c(.20, .10, .30, .30, .10), c(.25, .25, 0, .25, .25))
  median <- expect_plot(plot(fit, statistic = "median", intervals = FALSE))$data
  mode <- expect_plot(plot(fit, statistic = "mode", intervals = FALSE))$data
  expect_equal(median$mean[median$position == 1L], c(3, 2))
  expect_equal(mode$mean[mode$position == 1L], c(2, 1))
})

test_that("categorical summary plots never trigger parameter inference", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture()
  local_mocked_bindings(parameter_inference = function(...) {
    stop("Plot unexpectedly requested parameter inference")
  })
  invisible(lapply(c("mean", "median", "mode"), function(statistic) {
    plotted <- expect_plot(plot(fit, statistic = statistic, intervals = TRUE))
    expect_true(all(is.na(plotted$data$lower)))
    expect_true(all(is.na(plotted$data$upper)))
    expect_match(plotted$labels$subtitle, "intervals.*unavailable",
                 ignore.case = TRUE)
  }))
})

test_that("one-class and FIML categorical fits retain their estimated summaries", {
  skip_if_not_installed("ggplot2")
  single <- .categorical_summary_fixture(n_profiles = 1L)
  expect_equal(nrow(expect_plot(plot(single, intervals = FALSE))$data), 2L)
  fit <- .categorical_summary_fixture(missing = TRUE)
  expect_true(anyNA(fit$categorical_data))
  expected <- vapply(fit$response_probabilities, function(probabilities) {
    drop(probabilities %*% as.numeric(colnames(probabilities)))
  }, numeric(fit$n_profiles))
  drawn <- expect_plot(plot(fit, intervals = FALSE))$data
  expect_equal(drawn$mean, unname(expected[cbind(drawn$profile, drawn$position)]))
  expect_true(all(is.finite(drawn$mean)))
  one_class_fiml <- .categorical_summary_fixture(n_profiles = 1L, missing = TRUE)
  empirical <- colMeans(one_class_fiml$categorical_data, na.rm = TRUE)
  drawn <- expect_plot(plot(one_class_fiml, intervals = FALSE))$data
  expect_equal(drawn$mean, unname(empirical[drawn$position]), tolerance = 1e-8)
})

test_that("categorical profiles reject ambiguous scores and standardized units", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture()
  expect_error(plot(fit, scale = "standardized", intervals = FALSE),
               class = "latents_bad_argument")
  invalid <- list(c("low", "medium", "high", "top", "bottom"),
                  c("1", "2", "3", "4", "Inf"),
                  c("1", "2", "3", "4", "NaN"),
                  c("1", "1.0", "3", "4", "5"))
  invisible(lapply(invalid, function(scores) {
    broken <- fit
    broken$categorical_levels$a <- scores
    colnames(broken$response_probabilities$a) <- scores
    expect_error(plot(broken, intervals = FALSE), "numeric|finite|distinct|unique",
                 class = "latents_bad_argument")
    expect_false("profiles" %in% .multilpa_supported_views(broken))
  }))
  expect_true("profiles" %in% .multilpa_supported_views(fit))
})

test_that("responses preserve their probabilities and refuse summary statistics", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture()
  drawn <- expect_plot(plot(fit, what = "responses", intervals = FALSE))$data
  expected <- vapply(fit$response_probabilities, function(probabilities) {
    probabilities[, "5"]
  }, numeric(fit$n_profiles))
  expect_equal(drawn$mean, unname(expected[cbind(drawn$profile, drawn$position)]))
  expect_error(plot(fit, what = "responses", statistic = "median"),
               class = "latents_bad_argument")
  expect_error(plot(fit, what = "sizes", statistic = "mode"),
               class = "latents_bad_argument")
})

test_that("mixed profiles preserve the continuous-only view for all statistics", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture(mixed = TRUE)
  # Nominal categorical indicators must not invalidate a continuous profile.
  fit$categorical_levels$a <- letters[1:5]
  colnames(fit$response_probabilities$a) <- letters[1:5]
  invisible(lapply(c("mean", "median", "mode"), function(statistic) {
    drawn <- expect_plot(plot(fit, statistic = statistic, intervals = FALSE))$data
    expect_equal(nrow(drawn), fit$n_profiles)
    expect_equal(drawn$mean, unname(fit$means[drawn$profile, "score"]))
    standardized <- expect_plot(plot(fit, statistic = statistic,
      scale = "standardized", intervals = FALSE))$data
    observed <- fit$indicator_data[, "score"]
    expect_equal(standardized$mean,
      (unname(fit$means[standardized$profile, "score"]) - mean(observed)) /
        stats::sd(observed))
  }))
})

test_that("membership covariate fits share the categorical summary API", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture(covariates = TRUE)
  expect_s3_class(fit, "multilpa_covariates")
  invisible(lapply(c("mean", "median", "mode"), function(statistic) {
    expect_plot(plot(fit, what = "profiles", statistic = statistic,
                     intervals = FALSE))
  }))
  expect_error(plot(fit, what = "responses", statistic = "median"),
               class = "latents_bad_argument")
})

test_that("all views apply a requested statistic only to the profiles view", {
  skip_if_not_installed("ggplot2")
  fit <- .categorical_summary_fixture()
  plots <- expect_plots(plot(fit, what = "all", statistic = "median",
                            intervals = FALSE))
  expect_true("profiles" %in% names(plots))
  expect_equal(plots$profiles$data,
               plot(fit, statistic = "median", intervals = FALSE)$data)
  expect_equal(plots$responses$data,
               plot(fit, what = "responses", intervals = FALSE)$data)
})

test_that("Gaussian profile statistics retain mean estimates and Wald intervals", {
  skip_if_not_installed("ggplot2")
  set.seed(56)
  dat <- data.frame(group = rep(1:12, each = 10L),
                    a = stats::rnorm(120L), b = stats::rnorm(120L))
  fit <- multilpa(dat, c("a", "b"), "group", 1L, 1L,
                  n_starts = 1L, seed = 1L)
  invisible(lapply(c("mean", "median", "mode"), function(statistic) {
    plotted <- expect_plot(plot(fit, statistic = statistic))
    drawn <- plotted$data
    expect_equal(drawn$mean, unname(colMeans(dat[c("a", "b")])[drawn$position]))
    expect_true(all(is.finite(drawn$lower)))
    expect_true(all(drawn$lower < drawn$mean & drawn$upper > drawn$mean))
    expect_equal(drawn, plot(fit, statistic = "mean")$data)
  }))
})
