# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Independent scalar Gaussian references for the growth engine's numerical
# fallback at large random-effect/residual variance ratios.
test_that("the growth quadratic preserves the random-intercept within-person scatter", {
  amplitudes <- c(1e3, 1e5, 1e7, 1e8)
  labels <- sprintf("amplitude %g", amplitudes)
  # The package's class density of three observations y = A + (-1, 0, 1)
  # with random-intercept variance A^2 and unit residual variance.
  engine <- vapply(amplitudes, function(amplitude) {
    y <- amplitude + c(-1, 0, 1)
    spec <- list(x = matrix(1, 3, 1), z = matrix(numeric(), 3, 0),
      random_design = matrix(1, 3, 1), y = y, offset = rep(0, 3), n = 3,
      group_index = rep(1L, 3))
    stats <- latents:::.growth_statistics(spec)
    latents:::.growth_class_terms(stats, 0, matrix(amplitude, 1, 1), 1)$log_density
  }, numeric(1))
  # The same density in closed form: the within-person scatter plus the mean
  # shrunk by 1 + 3 A^2, and log det V = log(1 + 3 A^2).
  closed_form <- vapply(amplitudes, function(amplitude) {
    y <- amplitude + c(-1, 0, 1)
    quadratic <- sum((y - mean(y))^2) + 3 * mean(y)^2 / (1 + 3 * amplitude^2)
    -0.5 * (3 * log(2 * pi) + log(1 + 3 * amplitude^2) + quadratic)
  }, numeric(1))
  expect_reference(stats::setNames(engine, labels), stats::setNames(closed_form, labels),
                   reference = "independent: closed-form random-intercept Gaussian log density",
                   covers = "mixture_regression()", quantity = "class log density",
                   tolerance = 1e-12, measure = "relative",
                   note = "Package side: the internal class density (.growth_class_terms) of the growth engine behind mixture_regression(random = ).")
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
  labels <- c("intercept", "log residual SD", "log random-intercept SD")
  # The residual contrast variance is at its optimum (SS=2, df=2), while
  # mean(y)^2 / G - 1 = 4 - 1 is the log random-SD derivative.
  # Scaled: the intercept score is about 2e-8 and the residual score is zero,
  # where a relative measure is undefined or judges rounding noise. 1e-10,
  # the original tolerance.
  expect_reference(stats::setNames(unname(scores[1, 1:3]), labels),
                   stats::setNames(c(3 * mean(spec$y) / (1 + 3 * amplitude^2), 0, 3), labels),
                   reference = "independent: closed-form scores of a three-observation random-intercept model",
                   covers = "mixture_regression()", quantity = "scores at a random-intercept SD of 1e8",
                   tolerance = 1e-10, measure = "scaled",
                   note = "Package side: the internal analytic scores (.growth_scores) behind vcov() of mixture_regression(random = ).")
})
