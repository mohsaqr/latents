.review_zero_weight_data <- function() {
  set.seed(1020)
  data <- data.frame(g = rep(seq_len(20), each = 5), t = rep(0:4, 20))
  data$y <- rep(stats::rnorm(20), each = 5) + 0.4 * data$t +
    stats::rnorm(100, sd = 0.7)
  data$weight <- rep(rep(c(1, 2), 10), each = 5)
  data
}

.review_zero_weight_fitters <- function() {
  list(
    additive = function(data) multilpa(data, "y", "g", n_group_classes = 1,
      family = "additive", weights = "weight", n_starts = 1, seed = 1, tol = 1e-10),
    regression = function(data) mixture_regression(y ~ t, data, n_classes = 1,
      id = "g", weights = "weight", vcov_type = "none", n_starts = 1,
      seed = 1, tol = 1e-10),
    growth = function(data) mixture_regression(y ~ t, data, n_classes = 1,
      id = "g", class_level = "group", random = ~ 1, weights = "weight",
      vcov_type = "none", n_starts = 1, seed = 1, tol = 1e-10))
}

test_that("zero-weight groups do not alter robust covariance across model families", {
  data <- .review_zero_weight_data()
  added <- data[data$g <= 5, ]
  added$g <- added$g + 20
  added$weight <- 0
  padded <- rbind(data, added)
  invisible(lapply(.review_zero_weight_fitters(), function(fit) {
    original <- fit(data)
    expanded <- fit(padded)
    expect_true(original$converged)
    expect_true(expanded$converged)
    expect_equal(unname(coef(expanded)), unname(coef(original)), tolerance = 2e-6)
    expect_equal(vcov(expanded, type = "robust"), vcov(original, type = "robust"),
                 tolerance = 2e-6)
  }))
})

test_that("zero-weight groups do not satisfy minimum independent unit requirements", {
  data <- .review_zero_weight_data()
  data <- data[data$g <= 3, ]
  data$y <- rep(c(-3, 0, 3), each = 5) + 0.4 * data$t +
    rep(c(-0.3, 0.4, -0.6, 0.8, -0.1), 3)
  added <- data
  added$g <- added$g + 3
  added$weight <- 0
  padded <- rbind(data, added)
  invisible(lapply(.review_zero_weight_fitters(), function(fit) {
    original <- fit(data)
    expanded <- fit(padded)
    expect_error(vcov(original, type = "robust"), class = "latents_too_few_groups")
    expect_error(vcov(expanded, type = "robust"), class = "latents_too_few_groups")
  }))
})
