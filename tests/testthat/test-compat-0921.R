# Code written for latents 0.9.21, the release submitted to CRAN, keeps
# working. Renamed functions and arguments keep their former names.

test_that("candidate_fit() is get_candidate() under its 0.9.21 name", {
  set.seed(12)
  d <- data.frame(g = rep(1:30, each = 10), y = c(rnorm(150, -2), rnorm(150, 2)))
  grid <- enumerate_classes(d, "y", "g", n_profiles = 1:2, n_group_classes = 1:2,
                            model = "VVI", n_starts = 1, seed = 12)
  expect_identical(formals(candidate_fit), formals(get_candidate))
  expect_identical(candidate_fit(grid, 2, 1), get_candidate(grid, 2, 1))
  expect_identical(candidate_fit(grid, n_profiles = 2, n_group_classes = 2,
                                 model = "VVI"),
                   get_candidate(grid, 2, 2, "VVI"))
  expect_error(candidate_fit(grid, 1, 2), class = "latents_failed_candidate")
  expect_error(candidate_fit(grid, 5, 1), class = "latents_unknown_candidate")
})

test_that("a former argument name and its present name together are refused", {
  expect_identical(.latents_former_argument("aic", FALSE, NULL, c("what", "criterion")),
                   "aic")
  expect_identical(.latents_former_argument("aic", FALSE, "bic", c("what", "criterion")),
                   "bic")
  expect_error(.latents_former_argument("aic", TRUE, "bic", c("what", "criterion")),
               class = "latents_bad_argument")
})

test_that("an ignored argument warns by class and names itself", {
  expect_silent(.latents_warn_unused_arguments(list(), "plot()", ""))
  expect_warning(.latents_warn_unused_arguments(list(colour = "red", 1), "plot()",
                                                "Use ggplot2."),
                 "`colour`, an unnamed argument", class = "latents_unused_argument")
})
