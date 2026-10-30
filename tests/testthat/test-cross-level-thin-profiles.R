# The cross-level likelihood has two Gaussian measurement blocks; each is
# checked using the normalized sampling weights used by its own M-step.

#' Capture cross-level warnings without hiding unrelated condition classes
#' @param ... Arguments passed to multilpa().
#' @return The fitted model and its warning conditions.
skip_on_cran()

cross_level_capture_warnings <- function(...) {
  notices <- list()
  fit <- withCallingHandlers(multilpa(...), warning = function(notice) {
    notices[[length(notices) + 1L]] <<- notice
    invokeRestart("muffleWarning")
  })
  list(fit = fit, warnings = notices)
}

#' Separate a few groups at both Gaussian levels
#' @param n_small Number of separated groups.
#' @param group_sizes Number of ratings per group.
#' @return A data frame with a group ID and two Gaussian indicators.
cross_level_thin_data <- function(n_small = 3L, group_sizes = rep(2L, 40L)) {
  stopifnot(n_small > 0L, n_small < length(group_sizes), all(group_sizes > 0L))
  set.seed(7)
  data <- data.frame(g = rep(seq_along(group_sizes), group_sizes),
                     a = stats::rnorm(sum(group_sizes)),
                     b = stats::rnorm(sum(group_sizes)))
  data[data$g <= n_small, c("a", "b")] <-
    data[data$g <= n_small, c("a", "b")] + 20
  data$w <- ifelse(data$g <= n_small, 0.01, 1)
  data
}

#' Keep only thin-profile conditions for their separate assertions
#' @param captured Result of cross_level_capture_warnings().
#' @return The thin-profile warning conditions.
cross_level_thin_notices <- function(captured) {
  Filter(function(notice) inherits(notice, "latents_thin_profiles"), captured$warnings)
}

test_that("cross-level profile support restores normalized sampling weights", {
  data <- cross_level_thin_data()
  invisible(lapply(c("full_cross_level", "restricted_cross_level"), function(family) {
    captured <- cross_level_capture_warnings(
      data, c("a", "b"), "g", 2L, 1L, family = family, weights = "w",
      n_starts = 1L, seed = 1, max_iter = 300L)
    fit <- captured$fit
    notices <- cross_level_thin_notices(captured)
    expect_true(fit$converged)
    expect_false(fit$boundary)
    expect_length(notices, 1L)
    expect_match(conditionMessage(notices[[1L]]), "at the individual level")
    expect_gt(min(colSums(fit$subject_posteriors)), 4)
    weights <- unname(fit$sampling_weights[fit$sufficient_statistics$index])
    expect_lt(min(colSums(fit$subject_posteriors * weights)), 1)
    expect_equal(rowSums(fit$subject_posteriors), rep(1, nrow(data)))
    expect_equal(unname(fit$group_posteriors), matrix(1, 40L, 1L))
  }))
})

test_that("cross-level group Gaussian mixtures are checked exactly once", {
  data <- cross_level_thin_data(group_sizes = rep(3L, 40L))
  invisible(lapply(c("full_cross_level", "restricted_cross_level"), function(family) {
    captured <- cross_level_capture_warnings(
      data, c("a", "b"), "g", 1L, 2L, family = family,
      between_variance = "varying", n_starts = 1L, seed = 1, max_iter = 300L)
    fit <- captured$fit
    notices <- cross_level_thin_notices(captured)
    expect_true(fit$converged)
    expect_false(fit$boundary)
    expect_length(notices, 1L)
    expect_match(conditionMessage(notices[[1L]]), "at the group level \\(group classes\\)")
    expect_equal(sort(unname(colSums(fit$group_posteriors))), c(3, 37))
    expect_equal(unname(fit$subject_posteriors), matrix(1, nrow(data), 1L))
  }))
})

test_that("cross-level group support also restores sampling weights", {
  data <- cross_level_thin_data(n_small = 6L)
  captured <- cross_level_capture_warnings(
    data, c("a", "b"), "g", 1L, 2L, family = "full_cross_level", weights = "w",
    n_starts = 1L, seed = 1, max_iter = 300L)
  fit <- captured$fit
  notices <- cross_level_thin_notices(captured)
  expect_true(fit$converged)
  expect_length(notices, 1L)
  expect_match(conditionMessage(notices[[1L]]), "at the group level")
  expect_gt(min(colSums(fit$group_posteriors)), 4)
  expect_lt(min(colSums(fit$group_posteriors * fit$sampling_weights)), 1)
})

test_that("shared cross-level spreads do not raise a thin-profile notice", {
  data <- cross_level_thin_data(group_sizes = rep(3L, 40L))
  invisible(lapply(c("full_cross_level", "restricted_cross_level"), function(family) {
    captured <- cross_level_capture_warnings(
      data, c("a", "b"), "g", 2L, 2L, family = family,
      variance_model = "equal", between_variance = "equal",
      n_starts = 1L, seed = 1, max_iter = 300L)
    expect_true(captured$fit$converged)
    expect_length(cross_level_thin_notices(captured), 0L)
  }))
})

test_that("restricted support uses each factor's own weight normalization", {
  data <- cross_level_thin_data(group_sizes = c(rep(10L, 3L), rep(2L, 37L)))
  captured <- cross_level_capture_warnings(
    data, c("a", "b"), "g", 2L, 1L, family = "restricted_cross_level", weights = "w",
    n_starts = 1L, seed = 1, max_iter = 300L)
  fit <- captured$fit
  individual <- fit$gaussian_diagnostics[
    fit$gaussian_diagnostics$context == "individual profiles", , drop = FALSE]
  group <- fit$gaussian_diagnostics[
    fit$gaussian_diagnostics$context == "group classes", , drop = FALSE]
  raw_weights <- data$w
  normalized_weights <- raw_weights * (nrow(data) / sum(raw_weights))
  expect_equal(unname(individual$gaussian_mass),
               unname(colSums(fit$subject_posteriors * normalized_weights)))
  expect_equal(sum(individual$gaussian_mass), nrow(data))
  expect_equal(unname(group$gaussian_mass), 40)
  expect_length(cross_level_thin_notices(captured), 1L)
})

test_that("excluded groups do not determine full cross-level starts", {
  data <- cross_level_thin_data(n_small = 15L)
  data$w <- 1
  data$w[data$g == 40L] <- 0
  far <- data
  far[far$g == 40L, c("a", "b")] <- 1e20
  baseline <- cross_level_capture_warnings(
    data, c("a", "b"), "g", 1L, 2L, family = "full_cross_level", weights = "w",
    n_starts = 2L, seed = 1, max_iter = 300L)$fit
  rng <- .Random.seed
  moved <- cross_level_capture_warnings(
    far, c("a", "b"), "g", 1L, 2L, family = "full_cross_level", weights = "w",
    n_starts = 2L, seed = 1, max_iter = 300L)$fit
  expect_true(baseline$converged)
  expect_true(moved$converged)
  expect_identical(.Random.seed, rng)
  fields <- c("profile_means", "profile_variances", "composition", "group_means",
              "group_mean_variances", "group_probabilities", "log_likelihood", "starts")
  expect_identical(moved[fields], baseline[fields])
  expect_equal(unname(rowSums(moved$group_posteriors)), rep(1, 40L))
})

test_that("excluded groups cannot conceal constant or insufficient positive support", {
  data <- cross_level_thin_data(n_small = 15L)
  data$a <- 0
  data$w <- 1
  data$w[data$g == 40L] <- 0
  data$a[data$g == 40L] <- 1e4
  # `a` is constant in every group that carries weight; only the zero-weight
  # group varies. Since round 2 of the 2026-10-26 review (ADD2-3) that is
  # refused as the same data without the zero-weight group is, instead of
  # fitting at the variance floor with a boundary warning.
  expect_error(multilpa(data, c("a", "b"), "g", 2L, 1L, family = "full_cross_level",
                        weights = "w", n_starts = 1L, seed = 1, max_iter = 300L),
               class = "latents_bad_data")
  data <- data.frame(g = rep(1:3, each = 2L),
                     a = c(0, 0, 0, 0, 3, 4),
                     b = c(0, 0, 0, 0, 3, 4), w = c(1, 1, 1, 1, 0, 0))
  expect_error(multilpa(data, c("a", "b"), "g", 2L, 1L,
                        family = "full_cross_level", weights = "w",
                        n_starts = 1L, seed = 1),
               class = "latents_bad_data")
})
