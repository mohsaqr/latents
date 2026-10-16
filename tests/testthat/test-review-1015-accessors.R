# Regressions for the 2026-10-15 review of the accessors, methods, reporting
# verbs, descriptives, prediction and the condition catalogue.

review_activity <- c("browse", "lectures", "forum_read")

review_two_level <- function(...) {
  quietly(multilpa(engagement_small, vars = review_activity, id = "student",
                   n_profiles = 2, n_group_classes = 2, n_starts = 2, seed = 1,
                   ...))
}

review_single_level <- function(...) {
  quietly(multilpa(engagement_small, vars = review_activity, id = NULL,
                   n_profiles = 2, n_starts = 2, seed = 1, ...))
}

test_that("covariate-fit recovery compares a group-level truth with group classes", {
  fit <- review_two_level(profile_covariates = "previous_grade")
  # Before the fix the covariate fit's absent `group_sizes` made `max()` warn
  # and sent the group-level truth to the profile comparison.
  recovery <- expect_no_warning(get_results(
    fit, "assignments", data = engagement_small,
    truth = c("engagement", "student_type")))
  expect_identical(unique(recovery$assignment[recovery$truth == "student_type"]),
                   "group_class")
  expect_identical(unique(recovery$assignment[recovery$truth == "engagement"]),
                   "profile")
  # Each group is counted once at the group level.
  expect_identical(sum(recovery$n[recovery$truth == "student_type"]),
                   length(unique(engagement_small$student)))
  # The proportions are shares within each truth value.
  shares <- tapply(recovery$proportion, paste(recovery$truth, recovery$value), sum)
  expect_equal(as.vector(shares), rep(1, length(shares)))
})

test_that("descriptives of a single-level fit report no fabricated groups", {
  fit <- review_single_level()
  described <- descriptives(fit)
  expect_false(any(c("n_groups", "icc") %in% names(described)))
  expect_identical(described$variable, review_activity)
  by_profile <- descriptives(fit, by = "profile")
  expect_false(any(c("n_groups", "icc") %in% names(by_profile)))
  # A two-level fit keeps its grouping summaries.
  expect_true(all(c("n_groups", "icc") %in% names(descriptives(review_two_level()))))
})

test_that("descriptives of a fit refuse arguments it cannot use, by class", {
  fit <- review_two_level()
  expect_error(descriptives(fit, vars = "browse"), class = "latents_bad_argument")
  expect_error(descriptives(fit, by = "profile", id = "course"),
               class = "latents_bad_argument")
})

test_that("diagnostics() and report() work on a fit with a noise component", {
  fit <- review_single_level(noise = TRUE)
  quality <- expect_no_error(diagnostics(fit, plots = FALSE))
  expect_false("residuals" %in% names(get_results(quality, "all")))
  expect_error(get_results(quality, "residuals"), class = "latents_unsupported_noise")
  expect_output(print(quality), "not available")
  expect_output(suppressMessages(report(fit, plots = FALSE, rows = 1L)),
                "Classification quality")
})

test_that("summary print names single-level and class models as print() does", {
  single <- capture.output(print(summary(review_single_level()), rows = 1L))
  expect_match(single[1L], "^Latent profile analysis: 2 profiles$")
  expect_false(any(grepl("group class", single[1:3])))
  two <- capture.output(print(summary(review_two_level()), rows = 1L))
  expect_match(two[1L], "^Two-level latent profile analysis: 2 profiles, 2 group classes$")
})

test_that("summary records a valid missing-data mode for older fits", {
  fit <- review_two_level()
  fit$missing <- NULL
  expect_identical(summary(fit)$missing, "error")
})

test_that("a zero-row print depth prints no empty frames", {
  printed <- capture.output(print(summary(review_two_level()), rows = 0L))
  expect_false(any(grepl("<0 rows>", printed, fixed = TRUE)))
  expect_true(any(grepl("more rows", printed, fixed = TRUE)))
  expect_false(any(grepl("<0 rows>", capture.output(print(review_two_level(), rows = 0L)),
                         fixed = TRUE)))
})

test_that("empty measurement tables keep their standard-error columns", {
  continuous <- review_two_level()
  responses <- get_results(continuous, "responses")
  expect_identical(nrow(responses), 0L)
  expect_identical(names(responses),
                   c("profile", "indicator", "category", "probability", "threshold",
                     "probability_standard_error"))
  frame <- engagement_small
  frame$posted <- factor(ifelse(frame$forum_post > 0, "yes", "no"))
  categorical <- quietly(multilpa(frame, vars = "posted", id = "student",
                                  categorical = "posted", n_profiles = 2,
                                  n_group_classes = 1, n_starts = 2, seed = 1))
  profiles <- suppressMessages(get_results(categorical, "profiles"))
  expect_identical(nrow(profiles), 0L)
  expect_identical(names(profiles),
                   names(suppressMessages(get_results(continuous, "profiles"))))
})

test_that("starting_values() names its own argument in refusals", {
  expect_error(starting_values(list(a = 1)), "^`x`", class = "latents_bad_start")
  expect_error(starting_values(review_two_level(), covariance = "keep"), "^`x`",
               class = "latents_bad_start")
})

test_that("predict() refuses covariate and transition fits by class", {
  covariate <- review_two_level(profile_covariates = "previous_grade")
  expect_error(predict(covariate), class = "latents_unsupported_prediction")
  transitions <- quietly(lta(engagement_small, vars = review_activity, id = "student",
                             time = "sequence", n_profiles = 2, n_group_classes = 1,
                             n_starts = 2, seed = 1))
  expect_error(predict(transitions), class = "latents_unsupported_prediction")
})

test_that("every condition class the package raises is in the catalogue", {
  source_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(source_dir), "package sources are not available")
  sources <- unlist(lapply(list.files(source_dir, "[.]R$", full.names = TRUE),
                           readLines, warn = FALSE))
  raised <- unlist(regmatches(sources, gregexpr(
    "class = (c[(])?\"(latents|multilpa)_[a-z_]+\"(, *\"(latents|multilpa)_[a-z_]+\")*",
    sources)))
  raised <- unique(unlist(regmatches(raised, gregexpr("(latents|multilpa)_[a-z_]+", raised))))
  # Object classes share the call shape; they are not conditions.
  objects <- c("latents_additive_data", "latents_family_enumeration",
               "latents_growth_mixture", "latents_mixture_regression",
               "latents_plots", "latents_pooled", "latents_prior",
               "latents_regression_enumeration", "latents_transition_enumeration",
               "multilpa_covariates", "multilpa_enumeration", "multilpa_transitions")
  catalogue <- readLines(file.path(source_dir, "conditions.R"), warn = FALSE)
  documented <- unique(unlist(regmatches(catalogue, gregexpr(
    "`(latents|multilpa)_[a-z_]+`", catalogue))))
  documented <- gsub("`", "", documented, fixed = TRUE)
  expect_identical(sort(setdiff(raised, c(objects, documented))), character())
})

test_that("descriptives of a data frame name an argument they do not use", {
  expect_warning(described <- descriptives(engagement_small, vars = "browse",
                                           group = "student"),
                 class = "latents_unused_argument")
  expect_identical(described, descriptives(engagement_small, vars = "browse"))
})
