# Review 2026-10-26, ML-1: a unit of sampling weight zero contributes nothing,
# so adding one changes no likelihood, no criterion and no sample size. Before,
# weights were normalised over every unit, and four zero-weight groups added to
# fifteen scaled the log likelihood by 19/15 and moved BIC from 622.99 to 784.74.

skip_on_cran()

zero_weight_data <- function() {
  withr::local_seed(1026)
  groups <- 16L
  data <- data.frame(g = rep(seq_len(groups), each = 6L), t = rep(0:5, groups))
  kind <- 1L + (data$g %% 2L)
  data$y <- c(-1.5, 1.5)[kind] + 0.3 * data$t + stats::rnorm(nrow(data), sd = 0.8)
  data$y2 <- c(1, -1)[kind] + stats::rnorm(nrow(data), sd = 0.8)
  data$w <- rep(c(1, 2), length.out = groups)[data$g]
  data
}

pad_with_zero_weights <- function(data) {
  added <- data[data$g <= 4L, ]
  added$g <- added$g + max(data$g)
  added$y <- added$y + 50
  added$y2 <- added$y2 - 50
  added$w <- 0
  rbind(data, added)
}

zero_weight_fitters <- list(
  multilpa = function(data) multilpa(data, c("y", "y2"), "g", n_profiles = 2,
                                     n_group_classes = 2, weights = "w", n_starts = 5,
                                     seed = 1, tol = 1e-10),
  lta = function(data) lta(data, c("y", "y2"), "g", n_profiles = 2, time = "t",
                           weights = "w", n_starts = 3, seed = 1, tol = 1e-10),
  additive = function(data) multilpa(data, "y", "g", n_group_classes = 1,
                                     family = "additive", weights = "w", n_starts = 1,
                                     seed = 1, tol = 1e-10),
  regression = function(data) mixture_regression(y ~ t, data, n_classes = 1, id = "g",
                                                 weights = "w", vcov_type = "none",
                                                 n_starts = 1, seed = 1, tol = 1e-10),
  growth = function(data) mixture_regression(y ~ t, data, n_classes = 1, id = "g",
                                             class_level = "group", random = ~ 1,
                                             weights = "w", vcov_type = "none",
                                             n_starts = 1, seed = 1, tol = 1e-10))

test_that("zero-weight units change no likelihood, criterion or sample size", {
  data <- zero_weight_data()
  padded <- pad_with_zero_weights(data)
  invisible(lapply(names(zero_weight_fitters), function(family) {
    fit <- zero_weight_fitters[[family]]
    original <- quietly(fit(data))
    expanded <- quietly(fit(padded))
    expect_equal(as.numeric(logLik(expanded)), as.numeric(logLik(original)),
                 tolerance = 1e-8, info = family)
    expect_equal(stats::BIC(expanded), stats::BIC(original), tolerance = 1e-8,
                 info = family)
    expect_identical(stats::nobs(expanded), stats::nobs(original), info = family)
  }))
})

test_that("zero-weight rows are not classified, but keep their posteriors", {
  data <- zero_weight_data()
  padded <- pad_with_zero_weights(data)
  original <- quietly(zero_weight_fitters$multilpa(data))
  expanded <- quietly(zero_weight_fitters$multilpa(padded))
  expect_equal(get_results(expanded, "entropy"), get_results(original, "entropy"),
               tolerance = 1e-8)
  expect_equal(get_results(expanded, "information_criteria"),
               get_results(original, "information_criteria"), tolerance = 1e-8)
  # The data still has its rows; every one of them keeps a posterior.
  expect_identical(expanded$n_observations, nrow(padded))
  expect_identical(nrow(get_results(expanded, "assignments")), nrow(padded))
})
