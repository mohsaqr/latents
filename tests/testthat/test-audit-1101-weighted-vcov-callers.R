skip_on_cran()

# Audit 2026-11-01, follow-up: a weighted fit asked for an observed-information
# (or OPG) covariance is refused (`latents_unsupported_weights`). Where the fit's
# family has a bootstrap that runs on weighted fits (transition models, growth
# mixtures), the refusal names `method = "bootstrap"` beside robust inference;
# a plain mixture regression, which has no bootstrap, names robust inference
# only.

weighted_vcov_refusal <- function(expr) {
  condition <- tryCatch(expr, latents_unsupported_weights = function(e) e)
  expect_s3_class(condition, "latents_unsupported_weights")
  conditionMessage(condition)
}

weighted_lta_data <- function(n = 60L, occasions = 3L) {
  set.seed(41)
  state <- sample.int(2L, n * occasions, replace = TRUE)
  data <- data.frame(id = rep(seq_len(n), each = occasions),
                     time = rep(seq_len(occasions), n))
  data$y1 <- c(0, 3)[state] + stats::rnorm(n * occasions)
  data$y2 <- c(0.5, 3.5)[state] + stats::rnorm(n * occasions)
  data$w <- rep(sample(1:3, n, replace = TRUE), each = occasions)
  data
}

weighted_growth_data <- function(n = 60L, times = 0:4) {
  set.seed(21)
  class <- sample(1:2, n, replace = TRUE)
  data <- data.frame(id = rep(seq_len(n), each = length(times)),
                     time = rep(times, n))
  data$y <- c(2, 6)[class[data$id]] + c(0.5, -0.4)[class[data$id]] * data$time +
    rep(stats::rnorm(n, 0, 0.8), each = length(times)) +
    stats::rnorm(nrow(data), 0, 0.6)
  data$w <- rep(sample(1:3, n, replace = TRUE), each = length(times))
  data
}

test_that("weighted transition-model refusals name the bootstrap", {
  skip_on_cran()
  data <- weighted_lta_data()
  fit <- suppressWarnings(lta(data, c("y1", "y2"), "id", 2, time = "time",
                              weights = "w", n_starts = 2, seed = 1))
  # The bootstrap the message names does run on this weighted fit.
  boot <- parameter_inference(fit, data = data, method = "bootstrap", iter = 3,
                              n_starts = 1, seed = 1)
  expect_true(all(is.finite(boot$standard_error)))
  messages <- c(
    parameter_inference = weighted_vcov_refusal(
      parameter_inference(fit, vcov_type = "observed")),
    get_results = weighted_vcov_refusal(
      get_results(fit, "transitions", vcov_type = "opg")),
    vcov = weighted_vcov_refusal(vcov(fit, type = "observed")))
  expect_true(all(grepl("method = \"bootstrap\"", messages, fixed = TRUE)))
  expect_true(all(grepl("vcov_type = \"robust\"", messages, fixed = TRUE)))
})

test_that("weighted growth-mixture refusals name the bootstrap", {
  skip_on_cran()
  data <- weighted_growth_data()
  fit <- suppressWarnings(mixture_regression(
    y ~ time, data, n_classes = 2, id = "id", class_level = "group",
    random = ~ 1, weights = "w", n_starts = 2, seed = 1))
  expect_s3_class(fit, "latents_growth_mixture")
  boot <- parameter_inference(fit, method = "bootstrap", iter = 3, n_starts = 1,
                              seed = 1)
  expect_true(all(is.finite(boot$standard_error)))
  messages <- c(
    parameter_inference = weighted_vcov_refusal(
      parameter_inference(fit, vcov_type = "observed")),
    vcov = weighted_vcov_refusal(vcov(fit, type = "opg")),
    fit = weighted_vcov_refusal(mixture_regression(
      y ~ time, data, n_classes = 2, id = "id", class_level = "group",
      random = ~ 1, weights = "w", vcov_type = "opg", n_starts = 1, seed = 1)))
  expect_true(all(grepl("method = \"bootstrap\"", messages, fixed = TRUE)))
})

test_that("plain weighted mixture-regression refusals do not name a bootstrap", {
  skip_on_cran()
  set.seed(17)
  data <- data.frame(x = stats::rnorm(60))
  data$y <- ifelse(stats::rbinom(60, 1, 0.4) == 1, 1 + 2 * data$x, -1 - data$x) +
    stats::rnorm(60, 0, 0.5)
  data$w <- sample(1:3, 60, replace = TRUE)
  message <- weighted_vcov_refusal(
    mixture_regression(y ~ x, data, 2, weights = "w", vcov_type = "opg",
                       n_starts = 1, seed = 1))
  expect_false(grepl("bootstrap", message, fixed = TRUE))
  expect_true(grepl("vcov_type = \"robust\"", message, fixed = TRUE))
})
