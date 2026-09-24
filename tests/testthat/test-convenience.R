# Exhaustive property checks: they run on CI and locally, and are skipped on
# CRAN to keep the check within its time budget.
skip_on_cran()

.activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

.convenience_fit <- function(...) {
  multilpa(engagement_small, .activity,
           id = "student", n_profiles = 2, n_group_classes = 2,
           n_starts = 4, seed = 1, ...)
}

test_that("a fit carries the data it was built from, exactly", {
  fit <- .convenience_fit()
  rebuilt <- get_results(fit, "data")
  expect_identical(rebuilt, engagement_small[c("student", .activity)])
  # Columns the model never saw are not invented.
  expect_false(any(c("course", "previous_grade", "engagement", "student_type") %in%
                     names(rebuilt)))
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
  expect_true(exact(multilpa(data, c("a", "b"), "school", 2, 2, n_starts = 1, seed = 1),
                    c("school", "a", "b")))
  # Categorical indicators are stored as codes and must come back as labels.
  expect_true(exact(multilpa(data, c("a", "b", "passed"), "school", 2, 2,
                             categorical = "passed", n_starts = 1, seed = 1),
                    c("school", "a", "b", "passed")))
  expect_true(exact(multilpa(data, c("a", "b"), "school", 2, 2, time = "wave",
                             n_starts = 1, seed = 1),
                    c("school", "wave", "a", "b")))
  expect_true(exact(lta(data, c("a", "b"), "school", n_profiles = 2,
                                    time = "wave", n_starts = 1, seed = 1),
                    c("school", "wave", "a", "b")))
  expect_true(exact(multilpa(data, c("a", "b"), "school", 2, 1,
                                   profile_covariates = "z", n_starts = 1, seed = 1),
                    c("school", "a", "b")))
})

test_that("mixed fits retain numeric indicators and typed categorical values", {
  set.seed(42)
  data <- data.frame(
    school = rep(seq_len(16L), each = 6L),
    wave = rep(seq_len(6L), times = 16L),
    a = as.integer(round(5 * stats::rnorm(96L))),
    b = stats::rnorm(96L),
    passed = ordered(sample(c("yes", "no"), 96L, TRUE),
                     levels = c("yes", "no")),
    rating = sample(c(10L, 20L), 96L, TRUE),
    z = stats::rnorm(96L))
  indicators <- c("a", "b", "passed", "rating")
  categorical <- c("passed", "rating")
  fit <- quietly(multilpa(data, indicators, "school", 2L, 1L,
                           categorical = categorical, n_starts = 2L,
                           max_iter = 80L, seed = 1L))
  centered <- quietly(multilpa(data, indicators, "school", 2L, 1L,
                                categorical = categorical, centering = "person",
                                n_starts = 1L, max_iter = 40L, seed = 1L))
  transition <- quietly(lta(data, indicators, "school", n_profiles = 2L,
                            time = "wave", categorical = categorical,
                            n_starts = 2L, max_iter = 80L, seed = 1L))
  covariate <- quietly(multilpa(data, indicators, "school", 2L, 1L,
                                 categorical = categorical,
                                 profile_covariates = "z", n_starts = 2L,
                                 max_iter = 80L, seed = 1L))
  expected <- data[c("school", indicators)]
  expect_identical(get_results(fit, "data"), expected)
  restored <- get_results(centered, "data")
  expect_identical(restored[c("school", "a", "passed", "rating")],
                   expected[c("school", "a", "passed", "rating")])
  expect_equal(restored$b, expected$b, tolerance = 1e-14)
  expect_identical(get_results(transition, "data"),
                   data[c("school", "wave", indicators)])
  expect_identical(get_results(covariate, "data"), expected)
  expect_type(covariate$indicator_data, "double")
  expect_identical(colnames(covariate$indicator_data), c("a", "b"))
})

test_that("inference no longer has to be handed data it already holds", {
  fit <- .convenience_fit()
  expect_equal(quietly(parameter_inference(fit)),
               quietly(parameter_inference(fit, engagement_small)))
  expect_identical(get_results(fit, "residuals"),
                   get_results(fit, "residuals", data = engagement_small))
  expect_equal(quietly(vcov(fit)),
               quietly(vcov(fit, engagement_small)))
  expect_equal(quietly(confint(fit)),
               quietly(confint(fit, data = engagement_small)))
  # A supplied frame is still validated against the fit, so the guard that
  # catches the wrong data has not been traded away for the convenience.
  wrong <- transform(engagement_small, browse = browse + 1)
  expect_error(quietly(parameter_inference(fit, wrong)),
               class = "latents_bad_inference_data")
})

test_that("descriptives describes what is there and says how it is nested", {
  described <- descriptives(engagement_small, vars = .activity, id = "student")
  expect_s3_class(described, "data.frame")
  expect_identical(described$variable, .activity)
  expect_identical(names(described),
                   c("variable", "n", "n_missing", "mean", "sd", "min", "max",
                     "n_distinct", "n_groups", "icc"))
  expect_true(all(described$n == nrow(engagement_small)))
  expect_true(all(described$n_groups == length(unique(engagement_small$student))))
  # Without an id there is nothing to report an icc against.
  plain <- descriptives(engagement_small, vars = "browse")
  expect_false(any(c("n_groups", "icc") %in% names(plain)))
  expect_error(descriptives(engagement_small, vars = "not_a_column"))
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
  expect_gt(descriptives(engagement_small, vars = "browse",
                         id = "student")$icc, 0.05)
})

test_that("splitting by profile shows the grouping the model claims to explain", {
  fit <- .convenience_fit()
  overall <- descriptives(fit)
  by_profile <- descriptives(fit, by = "profile")
  expect_identical(names(by_profile)[1L], "profile")
  # One row per indicator per profile: five indicators, two profiles.
  expect_identical(nrow(by_profile), 10L)
  expect_setequal(unique(by_profile$profile), c(1L, 2L))
  expect_identical(sum(by_profile$n), sum(overall$n))
  # The model's claim is that students differ in their MIX, not in the courses
  # they take. If that holds, conditioning on the profile leaves almost no
  # between-student variance behind.
  expect_gt(min(overall$icc), 0.05)
  # Every within-profile figure falls below the smallest overall one. Stated
  # that way rather than against a fixed threshold, because what is left inside
  # a profile is sampling noise around zero -- it is free to come out slightly
  # negative -- and its exact size is not the claim being made.
  expect_lt(max(abs(by_profile$icc)), min(overall$icc))
  expect_error(descriptives(fit, by = "nonsense"))
})

test_that("diagnostics gathers the verbs without changing what they return", {
  fit <- .convenience_fit()
  quality <- diagnostics(fit)
  expect_s3_class(quality, "multilpa_diagnostics")
  expect_identical(get_results(quality, "entropy"), get_results(fit, "entropy"))
  expect_identical(get_results(quality, "classification"),
                   get_results(fit, "classification", level = "both"))
  expect_identical(get_results(quality, "average_posteriors"),
                   get_results(fit, "average_posteriors", level = "both"))
  expect_identical(get_results(quality, "residuals"),
                   get_results(fit, "residuals"))
  expect_error(get_results(quality, "nonsense"))
  expect_output(print(quality), "Classification quality")
  expect_invisible(print(quality))
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
    expect_identical(quietly(report(fit, plots = FALSE)), fit)
    expect_output(report(fit, plots = FALSE), "Classification quality")
    expect_output(report(fit, plots = FALSE), "-- profiles", fixed = TRUE)
    # Drawing the bars computes intervals from the fit's own data, so a loosely
    # converged fit is told so even through a convenience verb. The warning is
    # the package doing its job; it is asserted, not silenced.
    #
    # This previously used the default-tolerance fit above and so pinned a
    # miscalibration: the check cut a fixed 0.01 on the scaled score, which on
    # a sample this size fires at a displacement of a fraction of a percent of
    # a standard error. The criterion is now that displacement itself, so the
    # fixture has to be genuinely loose to trip it.
    loose <- .convenience_fit(tol = 1e-4)
    expect_warning(report(loose, plots = TRUE),
                   class = "latents_unconverged")
    tight <- .convenience_fit(tol = 1e-10)
    # No warning, rather than no noise: report() prints by design, so
    # expect_silent() would fail on its own output.
    expect_warning(utils::capture.output(report(tight, plots = TRUE)),
                   regexp = NA)
    # And the ordinary default-tolerance fit is quiet too, which is the point of
    # the recalibration: the happy path must not carry a qualification.
    expect_warning(utils::capture.output(report(fit, plots = TRUE)),
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
  timed <- multilpa(engagement_small, .activity,
                    id = "student", n_profiles = 2, n_group_classes = 2,
                    time = "sequence", n_starts = 2, seed = 1)
  expect_true("sequences" %in% .multilpa_supported_views(timed))
})

test_that("assignments put the estimate and the truth in the same row", {
  fit <- .convenience_fit(tol = 1e-10)
  carried <- get_results(fit, "assignments")
  expect_identical(nrow(carried), nrow(engagement_small))
  expect_true(all(c("profile", "group_class", "posterior_profile_1",
                    "posterior_profile_2") %in% names(carried)))
  expect_identical(carried$profile, fit$subject_profiles)
  # The rows the model saw come back untouched beside the assignments.
  expect_identical(carried[c("student", "browse")],
                   engagement_small[c("student", "browse")])

  # Supplying a frame keeps every column of it, including ones the model never
  # saw, which is the whole point: the comparison needs them in the same row.
  supplied <- get_results(fit, "assignments", data = engagement_small)
  expect_true(all(names(engagement_small) %in% names(supplied)))
  expect_identical(supplied$engagement, engagement_small$engagement)
  expect_identical(supplied$profile, carried$profile)

  # The recovery check the help page makes. Which fitted label lands on which
  # generated level is the optimiser's business, so the agreement is read off
  # both diagonals and the better one taken.
  agreement <- xtabs(~ profile + engagement, data = supplied)
  expect_identical(sum(agreement), nrow(engagement_small))
  matched <- max(sum(diag(agreement)), sum(agreement) - sum(diag(agreement)))
  expect_gt(matched / sum(agreement), 0.97)
})

test_that("assignments refuse a frame that cannot be aligned", {
  fit <- .convenience_fit(tol = 1e-10)
  # Too few rows: silently recycling or truncating would be the alignment bug
  # this verb exists to prevent.
  expect_error(get_results(fit, "assignments", data = head(engagement_small, 10L)),
               class = "latents_bad_inference_data")
  # A column the assignments would overwrite is an error, not a replacement.
  expect_error(get_results(fit, "assignments", data = transform(engagement_small, profile = 1L)),
               class = "latents_bad_data")
})
