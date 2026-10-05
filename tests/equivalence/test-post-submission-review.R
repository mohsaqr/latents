# Independent references and edge cases missed by the post-submission additions.
review_regression_data <- function() {
  withr::with_seed(943, {
    data <- expand.grid(time = 0:3, id = 1:32)
    data$school <- (data$id - 1L) %/% 4L + 1L
    data$off <- 4 + 2 * data$time
    data$y <- 3 + 0.6 * data$time + data$off +
      rep(stats::rnorm(32), each = 4) + stats::rnorm(nrow(data), sd = 0.7)
    data$binary <- stats::rbinom(nrow(data), 1, stats::plogis(-0.4 + 0.3 * data$time))
    data$weight <- rep(seq_len(32), each = 4)
    data
  })
}

test_that("one-class binary regression agrees with an ordinary logistic GLM", {
  data <- review_regression_data()
  reference <- stats::glm(binary ~ time, data, family = stats::binomial())
  fit <- mixture_regression(binary ~ time, data, n_classes = 1,
                            family = "binomial", n_starts = 0)
  expect_equal(fit$log_likelihood, as.numeric(logLik(reference)), tolerance = 1e-9)
  expect_equal(unname(fit$params$beta[, 1]), unname(stats::coef(reference)),
               tolerance = 1e-7)
  expect_equal(get_results(fit)$std_error, sqrt(unname(diag(stats::vcov(reference)))),
               tolerance = 1e-5)
})

test_that("trajectory estimates include offsets and errors are on the response scale", {
  data <- review_regression_data()
  fit <- mixture_regression(y ~ time + offset(off), data, 1, id = "id",
                            class_level = "group", n_starts = 0)
  rows <- data[1:4, ]
  inference <- latents:::.mixture_resolve_inference(fit, NULL)
  band <- latents:::.trajectory_band(fit, rows, inference, 0.95)
  expect_equal(band$estimate, predict(fit, rows)$fitted, tolerance = 1e-12)
  logistic <- mixture_regression(binary ~ time, data, 1, family = "binomial",
                                 id = "id", class_level = "group", n_starts = 0)
  inference <- latents:::.mixture_resolve_inference(logistic, NULL)
  band <- latents:::.trajectory_band(logistic, rows, inference, 0.95)
  reference <- stats::glm(binary ~ time, data, family = stats::binomial())
  predicted <- stats::predict(reference, rows, type = "link", se.fit = TRUE)
  expected <- predicted$se.fit * stats::plogis(predicted$fit) *
    stats::plogis(predicted$fit, lower.tail = FALSE)
  expect_equal(band$std_error, unname(expected), tolerance = 1e-5)
})
