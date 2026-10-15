# Independent scalar Gaussian references for the growth engine's numerical
# fallback at large random-effect/residual variance ratios.
test_that("the growth quadratic preserves the random-intercept within-person scatter", {
  invisible(lapply(c(1e3, 1e5, 1e7, 1e8), function(amplitude) {
    y <- amplitude + c(-1, 0, 1)
    spec <- list(x = matrix(1, 3, 1), z = matrix(numeric(), 3, 0),
      random_design = matrix(1, 3, 1), y = y, offset = rep(0, 3), n = 3,
      group_index = rep(1L, 3))
    stats <- latents:::.growth_statistics(spec)
    result <- latents:::.growth_class_terms(stats, 0, matrix(amplitude, 1, 1), 1)
    quadratic <- sum((y - mean(y))^2) + 3 * mean(y)^2 / (1 + 3 * amplitude^2)
    reference <- -0.5 * (3 * log(2 * pi) + log(1 + 3 * amplitude^2) + quadratic)
    expect_equal(result$log_density, reference, tolerance = 1e-12)
  }))
})

test_that("extreme growth scores retain the random variance and residual derivatives", {
  amplitude <- 1e8
  spec <- list(x = matrix(1, 3, 1, dimnames = list(NULL, "(Intercept)")),
    z = matrix(numeric(), 3, 0), random_design = matrix(1, 3, 1,
      dimnames = list(NULL, "(Intercept)")), y = 2 * amplitude + c(-1, 0, 1),
    offset = rep(0, 3), n = 3, group_index = rep(1L, 3), n_groups = 1,
    n_classes = 1, variance = "varying", random_covariance = "equal",
    random_diagonal = TRUE, w = matrix(1, 1, 1, dimnames = list(NULL, "(Intercept)")),
    intercept_only = TRUE)
  params <- list(beta = matrix(0, 1, 1), common = numeric(), sigma2 = 1,
    random_covariance = list(matrix(amplitude^2, 1, 1)), random_scale = 1,
    gamma = matrix(0, 1, 1))
  scores <- latents:::.growth_scores(spec, latents:::.growth_statistics(spec), params)
  # The residual contrast variance is at its optimum (SS=2, df=2), while
  # mean(y)^2 / G - 1 = 4 - 1 is the log random-SD derivative.
  expect_equal(unname(scores[1, 1]), 3 * mean(spec$y) / (1 + 3 * amplitude^2), tolerance = 1e-10)
  expect_equal(unname(scores[1, 2]), 0, tolerance = 1e-10)
  expect_equal(unname(scores[1, 3]), 3, tolerance = 1e-10)
})
