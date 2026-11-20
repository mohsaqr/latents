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
  # stats::glm(family = binomial()): IRLS to its own deviance criterion.
  expect_reference(as.numeric(logLik(fit)), as.numeric(stats::logLik(reference)),
                   reference = "stats", covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one-class logistic regression",
                   tolerance = 1e-9, measure = "relative")
  # 1e-7 relative, the original tolerance: glm stops IRLS at a relative
  # deviance change of 1e-8, not at the exact optimum.
  expect_reference(coef(fit)[c("coefficient.class_1.(Intercept)", "coefficient.class_1.time")],
                   stats::coef(reference), reference = "stats",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "coefficients, one-class logistic regression",
                   tolerance = 1e-7, measure = "relative", key = c("(Intercept)", "time"))
  # 1e-5 relative, the original tolerance: glm's errors are the inverse
  # information at its IRLS stopping point, latents' at its own maximum.
  expect_reference(as.data.frame(get_results(fit))[c("term", "std_error")],
                   data.frame(term = c("(Intercept)", "time"),
                              std_error = sqrt(unname(diag(stats::vcov(reference))))),
                   reference = "stats",
                   covers = c("mixture_regression()", "get_results.latents_mixture_regression()",
                              "as.data.frame.latents_table()"),
                   quantity = "standard errors, one-class logistic regression",
                   tolerance = 1e-5, measure = "relative", key = "term")
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
  # The delta method on stats::predict.glm()'s link-scale errors; 1e-5
  # relative, the original tolerance (the two fits stop at their own criteria).
  expect_reference(band$std_error, unname(expected), reference = "stats",
                   covers = "mixture_regression()",
                   quantity = "response-scale standard errors of the trajectory, one-class logistic",
                   tolerance = 1e-5, measure = "relative",
                   key = sprintf("row %d", seq_len(nrow(rows))),
                   note = "Package side: the internal trajectory band (.trajectory_band) behind the trajectory tables and plot_trajectories() of a mixture_regression() fit; reference: stats::predict(type = \"link\", se.fit = TRUE) times the logistic density.")
})
