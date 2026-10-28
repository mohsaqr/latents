library(latents)

capture_fit <- function(expr) {
  classes <- character()
  fit <- withCallingHandlers(expr, warning = function(condition) {
    classes <<- c(classes, class(condition)[1L])
    invokeRestart("muffleWarning")
  }, message = function(condition) invokeRestart("muffleMessage"))
  list(fit = fit, classes = classes)
}

testthat::test_that("patterned covariances use their actual parameter count", {
  set.seed(73)
  data <- as.data.frame(matrix(rnorm(80 * 5), 80, 5))
  names(data) <- letters[1:5]
  data[1:13, ] <- data[1:13, ] + 10
  result <- capture_fit(lpa(data, names(data), 2, local_dependence = "a:b",
                             n_starts = 1, seed = 1))
  testthat::expect_equal(result$fit$n_parameters, 23)
  testthat::expect_equal(unname(result$fit$effective_profile_counts), c(13, 67))
  testthat::expect_false("latents_thin_profiles" %in% result$classes)
})

testthat::test_that("fixed means and variances cannot collapse during a fit", {
  set.seed(73)
  data <- as.data.frame(matrix(rnorm(80 * 5), 80, 5))
  names(data) <- letters[1:5]
  data[1:13, ] <- data[1:13, ] + 10
  data <- rbind(data[1:3, ], data[14:80, ])
  start <- list(means = rbind(rep(10, 5), rep(0, 5)),
                variances = matrix(1, 2, 5),
                profile_probabilities = matrix(c(.05, .95), 1),
                group_probabilities = 1)
  result <- capture_fit(multilpa(data, names(data), id = NULL, n_profiles = 2,
                                  start = start, fixed = c("means", "variances"),
                                  n_starts = 1))
  testthat::expect_equal(result$fit$n_parameters, 1)
  testthat::expect_false("latents_thin_profiles" %in% result$classes)
})

testthat::test_that("empty subjects cannot erase a thin measurement profile", {
  set.seed(1)
  data <- data.frame(g = rep(1:30, each = 3), t = rep(1:3, 30),
                     a = c(rnorm(87), rnorm(3, 20)),
                     b = c(rnorm(87), rnorm(3, 20)))
  empty <- data.frame(g = rep(31:130, each = 3), t = rep(1:3, 100),
                      a = NA_real_, b = NA_real_)
  fit <- function(frame) capture_fit(lta(frame, c("a", "b"), "g", 2, time = "t",
                                         missing = "fiml", n_starts = 1,
                                         seed = 1, max_iter = 500))
  original <- fit(data)
  appended <- fit(rbind(data, empty))
  testthat::expect_equal(appended$fit$effective_profile_counts,
                         original$fit$effective_profile_counts, tolerance = 1e-8)
  testthat::expect_true("latents_thin_profiles" %in% original$classes)
  testthat::expect_true("latents_thin_profiles" %in% appended$classes)
})

testthat::test_that("occasion measurement assesses memberships per block", {
  set.seed(1)
  data <- data.frame(g = rep(1:30, each = 3), t = rep(1:3, 30),
                     a = c(rnorm(84), rnorm(6, 20)),
                     b = c(rnorm(84), rnorm(6, 20)))
  result <- capture_fit(lta(data, c("a", "b"), "g", 2, time = "t",
                             measurement = "occasion", n_starts = 1,
                             seed = 1, max_iter = 200))
  counts <- rowsum(result$fit$subject_posteriors, data$t)
  testthat::expect_equal(unname(counts[, 2L]), rep(2, 3), tolerance = 1e-8)
  testthat::expect_true("latents_thin_profiles" %in% result$classes)
})

testthat::test_that("staged bootstrap retains the first stage screening pool", {
  set.seed(6)
  data <- data.frame(g = rep(1:12, each = 4), a = rnorm(48), b = rnorm(48))
  original <- capture_fit(fit_staged(data, c("a", "b"), "g", 2, 2,
                                      n_starts = 1, n_screen = 8,
                                      search = "emEM", seed = 1, max_iter = 20))$fit
  arguments <- latents:::.multilpa_staged_refit_arguments(original)
  control <- latents:::.latents_refit_controls(original, n_starts = 1)
  # Reproduce the argument assembly in .multilpa_bootstrap_replicates().
  refit <- capture_fit(do.call(fit_staged, c(list(data = data), arguments,
    list(n_starts = 1, max_iter = 20, tol = 1e-8, n_screen = control$n_screen))))$fit
  testthat::expect_identical(original$stage_one$n_screen, 8L)
  testthat::expect_identical(refit$stage_one$n_screen, original$stage_one$n_screen)
})
