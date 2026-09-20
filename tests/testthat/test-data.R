test_that("the bundled data have the shape they are documented with", {
  expect_s3_class(school_engagement, "data.frame")
  expect_identical(dim(school_engagement), c(720L, 6L))
  expect_identical(names(school_engagement),
                   c("school", "term", "engaged", "homework_hours",
                     "participation", "interest"))
  expect_false(anyNA(school_engagement))
  expect_type(school_engagement$engaged, "logical")
  # One row per student, so the school-term pair is a key, not a repeated measure.
  expect_equal(anyDuplicated(school_engagement[c("school", "term")]), 0L)
  expect_length(unique(school_engagement$school), 60L)

  expect_s3_class(engagement_panel, "data.frame")
  expect_identical(dim(engagement_panel), c(480L, 5L))
  expect_identical(names(engagement_panel),
                   c("student", "wave", "state", "homework_hours", "participation"))
  expect_false(anyNA(engagement_panel))
  expect_setequal(unique(engagement_panel$state), c(1L, 2L))
  expect_equal(anyDuplicated(engagement_panel[c("student", "wave")]), 0L)
  expect_length(unique(engagement_panel$student), 120L)
})

test_that("the schools differ in their mix, which is what the data are for", {
  share <- as.vector(tapply(school_engagement$engaged,
                            school_engagement$school, mean))
  # The documented hook: a pooled number near a half that describes no school.
  expect_equal(round(mean(school_engagement$engaged), 2), 0.48)
  expect_lt(min(share), 0.2)
  expect_gt(max(share), 0.9)
  # Two kinds of school, not a continuum: the shares are bimodal, so a
  # two-class group model has something to find.
  expect_gt(sum(share > 0.5), 20L)
  expect_gt(sum(share < 0.5), 20L)
})

test_that("a fit recovers the measurement model the data were generated from", {
  fit <- multilpa(school_engagement,
                  vars = c("homework_hours", "participation", "interest"),
                  id = "school", n_profiles = 2, n_group_classes = 2,
                  n_starts = 6, seed = 1)
  profiles <- as.data.frame(fit)
  # Profile labels are arbitrary, so the comparison is against the sorted
  # generating means rather than against a label.
  generated <- c(4.2, 8.0, 4.0, 7.2, 4.3, 7.4)
  fitted_means <- unlist(lapply(
    split(profiles$mean, profiles$indicator), sort), use.names = FALSE)
  expect_equal(sort(fitted_means), sort(generated), tolerance = 0.1)
  expect_equal(fit$n_profiles, 2L)
  expect_true(fit$converged)
})

test_that("a fit recovers the two-level structure, not just the profiles", {
  fit <- multilpa(school_engagement,
                  vars = c("homework_hours", "participation", "interest"),
                  id = "school", n_profiles = 2, n_group_classes = 2,
                  n_starts = 6, seed = 1)
  prevalence <- as.data.frame(fit, what = "profile_probabilities")
  # The whole claim of the package: the two group classes have opposite mixes,
  # so the lines cross. A pooled model would report one number near the middle.
  by_class <- tapply(prevalence$probability,
                     list(prevalence$group_class, prevalence$profile), identity)
  expect_gt(max(by_class[1L, ]), 0.65)
  expect_gt(max(by_class[2L, ]), 0.65)
  expect_false(identical(which.max(by_class[1L, ]), which.max(by_class[2L, ])))
})

test_that("the panel recovers the persistence it was generated with", {
  moves <- fit_transitions(engagement_panel,
                           vars = c("homework_hours", "participation"),
                           id = "student", time = "wave", n_profiles = 2,
                           n_starts = 6, seed = 1)
  stable <- transitions(moves, stable = TRUE)
  # Generated with persistence 0.85. Four waves of 120 students recovers it to
  # within a few points; asserting more than that would be asserting noise.
  staying <- stable$probability[stable$stable]
  expect_length(staying, 2L)
  expect_equal(mean(staying), 0.85, tolerance = 0.08)
})
