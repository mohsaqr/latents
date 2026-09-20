.convenience_fit <- function(...) {
  multilpa(school_engagement,
           c("homework_hours", "participation", "interest"),
           id = "school", n_profiles = 2, n_group_classes = 2,
           n_starts = 4, seed = 1, ...)
}

test_that("a fit carries the data it was built from, exactly", {
  fit <- .convenience_fit()
  rebuilt <- as.data.frame(fit, what = "data")
  expect_identical(rebuilt,
                   school_engagement[c("school", "homework_hours",
                                       "participation", "interest")])
  # Columns the model never saw are not invented.
  expect_false("engaged" %in% names(rebuilt))
})

test_that("the round trip is exact for every model family", {
  set.seed(4)
  data <- data.frame(school = rep(1:20, each = 8), wave = rep(1:8, times = 20),
                     a = stats::rnorm(160), b = stats::rnorm(160),
                     passed = sample(c("yes", "no"), 160, TRUE),
                     z = stats::rnorm(160), stringsAsFactors = FALSE)
  exact <- function(fit, columns) {
    identical(.multilpa_model_frame(fit)[columns], data[columns])
  }
  expect_true(exact(multilpa(data, c("a", "b"), "school", 2, 2, n_starts = 2, seed = 1),
                    c("school", "a", "b")))
  # Categorical indicators are stored as codes and must come back as labels.
  expect_true(exact(multilpa(data, c("a", "b", "passed"), "school", 2, 2,
                             categorical = "passed", n_starts = 2, seed = 1),
                    c("school", "a", "b", "passed")))
  expect_true(exact(multilpa(data, c("a", "b"), "school", 2, 2, time = "wave",
                             n_starts = 2, seed = 1),
                    c("school", "wave", "a", "b")))
  expect_true(exact(fit_transitions(data, c("a", "b"), "school", n_profiles = 2,
                                    time = "wave", n_starts = 2, seed = 1),
                    c("school", "wave", "a", "b")))
  expect_true(exact(fit_random_intercept(data, c("a", "b"), "school", 2,
                                         n_starts = 2, seed = 1),
                    c("school", "a", "b")))
  expect_true(exact(fit_covariates(data, c("a", "b"), "school", 2, 1,
                                   profile_covariates = "z", n_starts = 2, seed = 1),
                    c("school", "a", "b")))
})

test_that("inference no longer has to be handed data it already holds", {
  fit <- .convenience_fit()
  expect_equal(suppressWarnings(parameter_inference(fit)),
               suppressWarnings(parameter_inference(fit, school_engagement)))
  expect_identical(bivariate_residuals(fit),
                   bivariate_residuals(fit, school_engagement))
  expect_equal(suppressWarnings(vcov(fit)),
               suppressWarnings(vcov(fit, school_engagement)))
  expect_equal(suppressWarnings(confint(fit)),
               suppressWarnings(confint(fit, data = school_engagement)))
  # A supplied frame is still validated against the fit, so the guard that
  # catches the wrong data has not been traded away for the convenience.
  wrong <- transform(school_engagement, homework_hours = homework_hours + 1)
  expect_error(suppressWarnings(parameter_inference(fit, wrong)),
               class = "multilpa_bad_inference_data")
})

test_that("descriptives describes what is there and says how it is nested", {
  described <- descriptives(school_engagement,
                            vars = c("homework_hours", "participation", "interest"),
                            id = "school")
  expect_s3_class(described, "data.frame")
  expect_identical(described$variable,
                   c("homework_hours", "participation", "interest"))
  expect_identical(names(described),
                   c("variable", "n", "n_missing", "mean", "sd", "min", "max",
                     "n_distinct", "n_groups", "icc"))
  expect_true(all(described$n == 720L))
  expect_true(all(described$n_groups == 60L))
  # Without an id there is nothing to report an icc against.
  plain <- descriptives(school_engagement, vars = "homework_hours")
  expect_false(any(c("n_groups", "icc") %in% names(plain)))
  expect_error(descriptives(school_engagement, vars = "not_a_column"))
})

test_that("the icc reports group structure where there is some, and none where there is not", {
  # Data with no group structure at all: the estimator should sit at zero, and
  # is allowed to go slightly below it rather than being truncated.
  set.seed(2)
  flat <- data.frame(g = rep(seq_len(30L), each = 10L), noise = stats::rnorm(300L))
  expect_lt(abs(descriptives(flat, vars = "noise", id = "g")$icc), 0.05)

  # Data that is nothing but group structure: every observation equals its
  # group's own value, so all of the variance is between groups.
  pure <- data.frame(g = rep(seq_len(20L), each = 10L),
                     y = rep(stats::rnorm(20L), each = 10L))
  expect_gt(descriptives(pure, vars = "y", id = "g")$icc, 0.99)

  # And the real data sit in between, which is why a two-level model is worth
  # fitting to them at all.
  expect_gt(descriptives(school_engagement, vars = "homework_hours",
                         id = "school")$icc, 0.05)
})

test_that("splitting by profile shows the grouping the model claims to explain", {
  fit <- .convenience_fit()
  overall <- descriptives(fit)
  by_profile <- descriptives(fit, by = "profile")
  expect_identical(names(by_profile)[1L], "profile")
  expect_identical(nrow(by_profile), 6L)
  expect_setequal(unique(by_profile$profile), c(1L, 2L))
  expect_identical(sum(by_profile$n), sum(overall$n))
  # The model's claim is that schools differ in their MIX, not in their
  # students. If that holds, conditioning on the profile leaves almost no
  # between-school variance behind.
  expect_gt(min(overall$icc), 0.05)
  # Every within-profile figure falls below the smallest overall one. Stated
  # that way rather than against a fixed threshold, because what is left inside
  # a profile is sampling noise around zero -- with 60 groups its standard error
  # is around 0.18 -- and its exact size is not the claim being made.
  expect_lt(max(abs(by_profile$icc)), min(overall$icc))
  expect_error(descriptives(fit, by = "nonsense"))
})

test_that("diagnostics gathers the verbs without changing what they return", {
  fit <- .convenience_fit()
  quality <- diagnostics(fit)
  expect_s3_class(quality, "multilpa_diagnostics")
  expect_identical(as.data.frame(quality, what = "entropy"), entropy_table(fit))
  expect_identical(as.data.frame(quality, what = "classification"),
                   classification_table(fit, level = "both"))
  expect_identical(as.data.frame(quality, what = "posteriors"),
                   average_posteriors(fit, level = "both"))
  expect_identical(as.data.frame(quality, what = "residuals"),
                   bivariate_residuals(fit))
  expect_error(as.data.frame(quality, what = "nonsense"))
  expect_output(print(quality), "Classification quality")
  expect_invisible(print(quality))
})

test_that("a model family without residuals says so instead of failing", {
  fit <- fit_random_intercept(school_engagement,
                              c("homework_hours", "participation", "interest"),
                              id = "school", n_profiles = 2, n_starts = 2, seed = 1)
  quality <- diagnostics(fit)
  expect_s3_class(as.data.frame(quality, what = "entropy"), "data.frame")
  expect_error(as.data.frame(quality, what = "residuals"),
               class = "multilpa_no_group_classes")
  expect_output(print(quality), "not available for this model family")
})

test_that("diagnostics draws only when asked", {
  fit <- .convenience_fit()
  draw({
    expect_silent(diagnostics(fit))
    expect_s3_class(diagnostics(fit, plots = TRUE), "multilpa_diagnostics")
    expect_invisible(plot(diagnostics(fit)))
  })
  expect_error(diagnostics(fit, plots = "yes"), "must be TRUE or FALSE")
})

test_that("report prints everything and returns the fit", {
  fit <- .convenience_fit()
  draw({
    expect_identical(suppressWarnings(report(fit, plots = FALSE)), fit)
    expect_output(report(fit, plots = FALSE), "Classification quality")
    expect_output(report(fit, plots = FALSE), "Profile means")
    # Drawing the bars now computes intervals from the fit's own data, so a
    # loosely converged fit is told so even through a convenience verb. The
    # warning is the package doing its job; it is asserted, not silenced.
    expect_warning(report(fit, plots = TRUE), "non-negligible score")
    tight <- multilpa(school_engagement,
                      c("homework_hours", "participation", "interest"),
                      id = "school", n_profiles = 2, n_group_classes = 2,
                      n_starts = 4, tol = 1e-10, seed = 1)
    # No warning, rather than no noise: report() prints by design, so
    # expect_silent() would fail on its own output.
    expect_warning(utils::capture.output(report(tight, plots = TRUE)),
                   regexp = NA)
  })
})

test_that("report draws only the views a given fit can supply", {
  gaussian <- .convenience_fit()
  # No categorical indicators and no time variable, so neither view is offered
  # and report() does not stop on the first refusal.
  expect_false("responses" %in% .multilpa_supported_views(gaussian))
  expect_false("sequences" %in% .multilpa_supported_views(gaussian))
  expect_true(all(c("profiles", "bars", "heatmap", "probabilities",
                    "entropy", "posteriors") %in%
                    .multilpa_supported_views(gaussian)))
  timed <- multilpa(engagement_panel, c("homework_hours", "participation"),
                    id = "student", n_profiles = 2, n_group_classes = 2,
                    time = "wave", n_starts = 2, seed = 1)
  expect_true("sequences" %in% .multilpa_supported_views(timed))
})

test_that("assignments put the estimate and the truth in the same row", {
  fit <- multilpa(school_engagement,
                  c("homework_hours", "participation", "interest"),
                  id = "school", n_profiles = 2, n_group_classes = 2,
                  n_starts = 4, tol = 1e-10, seed = 1)
  carried <- assignments(fit)
  expect_identical(nrow(carried), 720L)
  expect_true(all(c("profile", "group_class", "posterior_profile_1",
                    "posterior_profile_2") %in% names(carried)))
  expect_identical(carried$profile, fit$subject_profiles)
  # The rows the model saw come back untouched beside the assignments.
  expect_identical(carried[c("school", "homework_hours")],
                   school_engagement[c("school", "homework_hours")])

  # Supplying a frame keeps every column of it, including ones the model never
  # saw, which is the whole point: the comparison needs them in the same row.
  supplied <- assignments(fit, data = school_engagement)
  expect_true(all(names(school_engagement) %in% names(supplied)))
  expect_identical(supplied$engaged, school_engagement$engaged)
  expect_identical(supplied$profile, carried$profile)

  # The recovery check the vignette makes: one student of 720 misassigned.
  agreement <- xtabs(~ profile + engaged, data = supplied)
  expect_identical(sum(agreement), 720L)
  expect_lt(sum(agreement) - sum(diag(agreement[, c("FALSE", "TRUE")])), 5L)
})

test_that("assignments refuse a frame that cannot be aligned", {
  fit <- multilpa(school_engagement,
                  c("homework_hours", "participation", "interest"),
                  id = "school", n_profiles = 2, n_group_classes = 2,
                  n_starts = 4, tol = 1e-10, seed = 1)
  # Too few rows: silently recycling or truncating would be the alignment bug
  # this verb exists to prevent.
  expect_error(assignments(fit, data = school_engagement[seq_len(10L), ]),
               class = "multilpa_bad_nesting")
  # A column the assignments would overwrite is an error, not a replacement.
  expect_error(assignments(fit, data = transform(school_engagement, profile = 1L)),
               class = "multilpa_bad_data")
})
