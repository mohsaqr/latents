#' Capture thin-profile warnings from a fitted transition model
#' @param expr A transition fitting expression.
#' @return The fit and its thin-profile warning conditions.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

capture_lta_thin <- function(expr) {
  thin <- list()
  fit <- withCallingHandlers(expr, warning = function(condition) {
    if (inherits(condition, "latents_thin_profiles")) {
      thin[[length(thin) + 1L]] <<- condition
    }
    if (inherits(condition, .multilpa_expected_warnings)) invokeRestart("muffleWarning")
  })
  list(fit = fit, thin = thin)
}

#' Generate separated profiles at three irregular time values
#' @param n_rare Number of rows in the rare profile.
#' @return A data frame with grouping, time and two Gaussian indicators.
lta_thin_data <- function(n_rare = 3L) {
  stopifnot(n_rare > 0L, n_rare < 90L)
  set.seed(1)
  data.frame(g = rep(seq_len(30L), each = 3L),
             t = rep(c(2, 4, 9), 30L),
             a = c(rnorm(90L - n_rare), rnorm(n_rare, 20)),
             b = c(rnorm(90L - n_rare), rnorm(n_rare, 20)))
}

test_that("empty LTA rows cannot hide thin Gaussian profiles under either empty policy", {
  data <- lta_thin_data()
  # A partly observed measurement still contributes to the density.
  data$a[nrow(data)] <- NA_real_
  empty <- data.frame(g = rep(31:130, each = 3L), t = rep(c(2, 4, 9), 100L),
                      a = NA_real_, b = NA_real_)
  invisible(lapply(c("homogeneous", "occasion"), function(transitions) {
    invisible(lapply(c(TRUE, FALSE), function(drop_empty) {
      fit <- function(frame) capture_lta_thin(lta(
        frame, c("a", "b"), "g", 2, time = "t", missing = "fiml",
        transitions = transitions, drop_empty = drop_empty,
        n_starts = 1, seed = 1, max_iter = 500))
      original <- fit(data)
      appended <- fit(rbind(data, empty))
      expect_length(original$thin, 1L)
      expect_length(appended$thin, 1L)
      expect_match(paste(vapply(appended$thin, conditionMessage, character(1)),
                         collapse = " "), "3.0 Gaussian membership mass for 4")
      expect_equal(sum(appended$fit$informative), 90L)
      expect_gt(min(colSums(appended$fit$subject_posteriors)), 10)
      expect_equal(appended$fit$log_likelihood, original$fit$log_likelihood,
                   tolerance = 1e-7)
    }))
  }))
})

test_that("occasion Gaussian densities warn separately with occasion positions", {
  data <- lta_thin_data(6L)
  fit <- function(measurement) capture_lta_thin(lta(
    data, c("a", "b"), "g", 2, time = "t", measurement = measurement,
    transitions = "occasion", n_starts = 1, seed = 1, max_iter = 500))
  occasion <- fit("occasion")
  shared <- fit("invariant")
  counts <- rowsum(occasion$fit$subject_posteriors, data$t)
  expect_equal(unname(counts[, 2L]), rep(2, 3L), tolerance = 1e-8)
  expect_length(occasion$thin, 3L)
  messages <- vapply(occasion$thin, conditionMessage, character(1))
  expect_true(all(vapply(seq_len(3L), function(label) {
    any(grepl(sprintf("at occasion %s:", label), messages, fixed = TRUE))
  }, logical(1))))
  expect_length(shared$thin, 0L)
})

test_that("empty rows cannot inflate occasion-specific measurement support", {
  data <- lta_thin_data(6L)
  data$a[nrow(data)] <- NA_real_
  empty <- data.frame(g = rep(31:130, each = 3L), t = rep(c(2, 4, 9), 100L),
                      a = NA_real_, b = NA_real_)
  invisible(lapply(c(TRUE, FALSE), function(drop_empty) {
    result <- capture_lta_thin(lta(
      rbind(data, empty), c("a", "b"), "g", 2, time = "t", missing = "fiml",
      measurement = "occasion", drop_empty = drop_empty,
      n_starts = 1, seed = 1, max_iter = 500))
    expect_length(result$thin, 3L)
    expect_true(all(grepl("2.0 Gaussian membership mass for 4",
                         vapply(result$thin, conditionMessage, character(1)),
                         fixed = TRUE)))
    expect_equal(sum(result$fit$informative), 90L)
  }))
})

test_that("LTA thin measurement warnings use scaled sampling memberships", {
  invisible(lapply(c("invariant", "occasion"), function(measurement) {
    # Nine raw rare memberships support a density, but their small weights do not.
    data <- lta_thin_data(9L)
    data$w <- ifelse(data$g >= 28L, 0.2, 1)
    result <- capture_lta_thin(lta(
      data, c("a", "b"), "g", 2, time = "t", weights = "w",
      measurement = measurement, n_starts = 1, seed = 1, max_iter = 500))
    weighted <- result$fit$subject_posteriors *
      unname(result$fit$sampling_weights[as.character(data$g)])
    expect_equal(min(colSums(result$fit$subject_posteriors)), 9, tolerance = 1e-7)
    expect_lt(min(colSums(weighted)), 4)
    expect_length(result$thin, if (measurement == "occasion") 3L else 1L)
    # Conversely, a well-supported weighted profile must not warn merely
    # because its unweighted membership is small.
    data <- lta_thin_data()
    data$w <- ifelse(data$g == 30L, 10, 1)
    result <- capture_lta_thin(lta(
      data, c("a", "b"), "g", 2, time = "t", weights = "w",
      measurement = measurement, n_starts = 1, seed = 1, max_iter = 500))
    expect_equal(min(colSums(result$fit$subject_posteriors)), 3, tolerance = 1e-7)
    expect_length(result$thin, 0L)
  }))
})

test_that("categorical-only LTA rows do not support a Gaussian density", {
  data <- lta_thin_data(9L)
  data$category <- factor(ifelse(data$g >= 28L, "rare", "common"))
  data[data$g >= 29L, c("a", "b")] <- NA_real_
  invisible(lapply(c("homogeneous", "occasion"), function(transitions) {
    result <- capture_lta_thin(lta(
      data, c("a", "b", "category"), "g", 2, time = "t",
      categorical = "category", missing = "fiml", transitions = transitions,
      n_starts = 1, seed = 1, max_iter = 500))
    expect_gt(min(colSums(result$fit$subject_posteriors)), 8.9)
    expect_true(all(result$fit$informative))
    expect_length(result$thin, 1L)
    expect_match(paste(vapply(result$thin, conditionMessage, character(1)),
                       collapse = " "), "3.0 Gaussian membership mass for 4")
  }))
})
