# Membership covariates became arguments to multilpa() in 0.11.2 rather than a
# separate verb. What matters is that the delegation forwards every argument
# unchanged, records the caller's own call, and refuses the three arguments the
# covariate path has no meaning for instead of accepting and ignoring them.

activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
# These tests check routing and arguments, not estimates, so a small subset
# and short runs are enough.
small <- engagement_small

test_that("naming a covariate returns a covariate fit, not a plain one", {
  fit <- quietly(multilpa(small, activity, "student", n_profiles = 2,
                  n_group_classes = 2, profile_covariates = "previous_grade",
                  n_starts = 1, max_iter = 20, seed = 1))
  expect_s3_class(fit, "multilpa_covariates")
  expect_false(inherits(fit, "multilpa"))
  expect_identical(fit$profile_covariates, "previous_grade")
  # The object reports the call the caller wrote, not the internal delegation.
  expect_identical(as.character(fit$call[[1L]]), "multilpa")
})

test_that("a group covariate alone also takes the covariate path", {
  # A group covariate must be constant within its group. The student's mean
  # prior grade is, and unlike `student_type` it does not separate the group
  # classes perfectly, so the logits stay finite and the fit converges.
  data <- course_engagement
  data$mean_grade <- stats::ave(data$previous_grade, data$student, FUN = mean)
  fit <- multilpa(data, activity, "student", n_profiles = 2,
                  n_group_classes = 2, group_covariates = "mean_grade",
                  n_starts = 2, max_iter = 4000, tol = 1e-8, seed = 1)
  expect_s3_class(fit, "multilpa_covariates")
  expect_identical(fit$group_covariates, "mean_grade")
  expect_length(fit$profile_covariates, 0L)
  expect_true(fit$converged)
  expect_true("group" %in% get_results(fit, "coefficients")$level)
})

test_that("no covariate leaves the covariate-free path untouched", {
  named <- quietly(multilpa(small, activity, "student", n_profiles = 2,
                    n_group_classes = 2, profile_covariates = character(),
                    group_covariates = character(), n_starts = 1,
                    max_iter = 20, seed = 1))
  plain <- quietly(multilpa(small, activity, "student", n_profiles = 2,
                    n_group_classes = 2, n_starts = 1, max_iter = 20,
                    seed = 1))
  expect_s3_class(named, "multilpa")
  expect_equal(named$log_likelihood, plain$log_likelihood)
  expect_equal(named$means, plain$means)
})

test_that("every forwarded argument reaches the estimator unchanged", {
  # A mis-wired argument in the delegation would change the fit, so the two
  # routes are compared on a specification that exercises the ones most easily
  # crossed: the covariance model, the variance floor and the tolerance.
  through_multilpa <- quietly(multilpa(
    small, activity, "student", n_profiles = 2,
    n_group_classes = 2, profile_covariates = "previous_grade",
    variance_model = "equal", covariance_model = "full", min_variance = 1e-4,
    tol = 1e-6, max_iter = 300L, n_starts = 1, seed = 1, time = "sequence"))
  direct <- quietly(multilpa:::.multilpa_fit_covariates(
    small, activity, "student", n_profiles = 2,
    n_group_classes = 2, profile_covariates = "previous_grade",
    variance_model = "equal", covariance_model = "full", min_variance = 1e-4,
    tol = 1e-6, max_iter = 300L, n_starts = 1, seed = 1, time = "sequence"))
  expect_equal(through_multilpa$log_likelihood, direct$log_likelihood)
  expect_equal(through_multilpa$means, direct$means)
  expect_equal(through_multilpa$profile_coefficients,
               direct$profile_coefficients)
  expect_equal(through_multilpa$covariances, direct$covariances)
  expect_identical(through_multilpa$variance_model, "equal")
  expect_identical(through_multilpa$covariance_model, "full")
  expect_identical(through_multilpa$time, "sequence")
})

test_that("the three arguments the covariate path cannot honour are refused", {
  call_with <- function(...) {
    multilpa(small, activity, "student", n_profiles = 2,
             n_group_classes = 2, profile_covariates = "previous_grade",
             n_starts = 1, max_iter = 20, seed = 1, ...)
  }
  plain <- quietly(multilpa(small, activity, "student", n_profiles = 2,
                    n_group_classes = 2, n_starts = 1, max_iter = 20,
                    seed = 1))
  expect_error(call_with(start = starting_values(plain)),
               class = "multilpa_bad_argument")
  expect_error(call_with(missing = "fiml"), class = "multilpa_bad_argument")
  expect_error(call_with(fixed = "means"), class = "multilpa_bad_argument")
  # The message names the argument, so the caller is told which one to drop.
  expect_error(call_with(missing = "fiml"), "`missing`")
  # All three at once are named together rather than one refusal at a time.
  expect_error(call_with(missing = "fiml", fixed = "means"), "`fixed`")
})

test_that("a non-character covariate name is refused before any fitting", {
  expect_error(multilpa(course_engagement, activity, "student", n_profiles = 2,
                        n_group_classes = 2, profile_covariates = 1),
               "must be a character vector")
  expect_error(multilpa(course_engagement, activity, "student", n_profiles = 2,
                        n_group_classes = 2, group_covariates = NA_character_),
               "must be a character vector")
})

test_that("a covariate fit offers the tables its family defines", {
  fit <- quietly(multilpa(small, activity, "student", n_profiles = 2,
                  n_group_classes = 2, profile_covariates = "previous_grade",
                  n_starts = 1, max_iter = 20, seed = 1))
  tables <- get_results(fit, "all")
  expect_true("coefficients" %in% names(tables))
  # A covariate model has no single profile prevalence: it varies with each
  # unit's covariates, so the table the covariate-free model has is absent.
  expect_false("profile_probabilities" %in% names(tables))
  expect_error(get_results(fit, "profile_probabilities"),
               class = "multilpa_bad_argument")
  expect_named(get_results(fit, "coefficients"),
               c("level", "outcome", "term", "parameter", "estimate"))
})
