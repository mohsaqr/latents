# The bundled dataset is `course_engagement`, and these tests are the contract
# `?course_engagement` states. Every documented claim is asserted here, so the
# help page cannot drift away from the data it describes.

.activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

# Persistence of the shipped `engagement` from one enrolment to the next,
# within a student. Brackets stay inside this helper rather than in the tests.
.observed_persistence <- function(data) {
  last <- nrow(data)
  within_student <- data$student[-1L] == data$student[-last]
  from <- data$engagement[-last][within_student]
  to <- data$engagement[-1L][within_student]
  as.vector(diag(prop.table(table(from, to), margin = 1L)))
}

# Share of the sample assigned to each fitted profile that agrees with a truth
# column, maximised over the arbitrary labelling of the fitted classes.
.label_free_agreement <- function(estimate, truth) {
  agreement <- table(estimate, truth)
  matched <- max(sum(diag(agreement)), sum(agreement) - sum(diag(agreement)))
  matched / sum(agreement)
}

.course_fit <- function(...) {
  multilpa(course_engagement, .activity, id = "student", n_profiles = 2,
           n_group_classes = 2, n_starts = 4, seed = 1, ...)
}

test_that("the bundled data have the shape they are documented with", {
  expect_s3_class(course_engagement, "data.frame")
  expect_identical(dim(course_engagement), c(1422L, 11L))
  expect_identical(names(course_engagement),
                   c("student", "course", "sequence", "browse", "lectures",
                     "forum_read", "forum_post", "attendance",
                     "previous_grade", "engagement", "student_type"))
  expect_false(anyNA(course_engagement))

  # The identifiers are integer, the measurements and the covariate numeric,
  # and the two generating truths factors.
  expect_type(course_engagement$student, "integer")
  expect_type(course_engagement$course, "integer")
  expect_type(course_engagement$sequence, "integer")
  expect_type(course_engagement$previous_grade, "double")
  expect_true(all(vapply(course_engagement[.activity], is.numeric, logical(1L))))
  # Both truth columns are factors, with the levels the help page lists and in
  # the order it lists them.
  expect_s3_class(course_engagement$engagement, "factor")
  expect_identical(levels(course_engagement$engagement),
                   c("disengaged", "engaged"))
  expect_s3_class(course_engagement$student_type, "factor")
  expect_identical(levels(course_engagement$student_type),
                   c("committed", "wavering"))
  # The truth column is deliberately NOT called `profile`: `assignments()` adds
  # a column of that name, and a dataset carrying it already would make the
  # verb refuse to join its own output onto the truth. Pinned here because the
  # help page states it and a workflow depends on it.
  expect_false("profile" %in% names(course_engagement))

  # 106 students and 32 courses, both numbered from one with no gaps.
  expect_identical(sort(unique(course_engagement$student)), seq_len(106L))
  expect_identical(sort(unique(course_engagement$course)), seq_len(32L))

  # One row per enrolment: a student takes a given course once, and occupies a
  # given position in their own order once.
  expect_equal(anyDuplicated(course_engagement[c("student", "course")]), 0L)
  expect_equal(anyDuplicated(course_engagement[c("student", "sequence")]), 0L)
  # And a student is one kind throughout, so `student_type` is a property of
  # the nesting unit rather than of the row.
  kinds <- vapply(split(as.character(course_engagement$student_type),
                        course_engagement$student),
                  function(kind) length(unique(kind)), integer(1L))
  expect_true(all(kinds == 1L))
})

test_that("the sequence runs 1..n inside every student, and the panel is ragged", {
  per_student <- split(course_engagement$sequence, course_engagement$student)
  clean <- vapply(per_student,
                  function(s) identical(s, seq_along(s)), logical(1L))
  expect_true(all(clean))

  taken <- vapply(per_student, length, integer(1L))
  expect_identical(sum(taken), 1422L)
  # Documented as "up to fifteen": the longest sequence reaches the cap and the
  # shortest falls well below it, which is what makes the panel ragged. A test
  # that only checked the maximum would pass on a balanced panel too.
  expect_identical(max(taken), 15L)
  expect_lt(min(taken), 15L)
  expect_gt(length(unique(taken)), 1L)
  expect_false(all(taken == taken[[1L]]))
})

test_that("the indicators are log counts and the covariate is standardised", {
  finite <- vapply(course_engagement[.activity],
                   function(v) all(is.finite(v)), logical(1L))
  non_negative <- vapply(course_engagement[.activity],
                         function(v) all(v >= 0), logical(1L))
  expect_true(all(finite))
  # A log count cannot be negative. The generator asserts this rather than
  # clamping, because a clamp would build the floor effect that looked like an
  # extra class; the shipped data must therefore satisfy it on its own.
  expect_true(all(non_negative))
  # The documented back-transform: `expm1()` reads an indicator as events, and
  # a `browse` of 3.9 is about 48 page views.
  expect_equal(round(expm1(3.9)), 48)

  # `previous_grade` is documented as standardised across the cohort.
  expect_equal(mean(course_engagement$previous_grade), 0, tolerance = 1e-3)
  expect_equal(stats::sd(course_engagement$previous_grade), 1, tolerance = 0.01)
})

test_that("the shipped truth columns describe the data they came with", {
  # `"engaged"` is documented as the engaged pattern, so every indicator must
  # be higher there. Asserting the direction, not just a difference, is what
  # keeps the labels in the help page attached to the right rows.
  separation <- vapply(.activity, function(indicator) {
    means <- tapply(course_engagement[[indicator]],
                    course_engagement$engagement, mean)
    as.vector(diff(means))
  }, numeric(1L))
  expect_true(all(separation > 0.5))

  engaged <- course_engagement$engagement == "engaged"
  share <- tapply(engaged, course_engagement$student_type, mean)
  # "A committed student is engaged in about four enrolments in five, a
  # wavering student in about one in four."
  expect_equal(share[["committed"]], 0.80, tolerance = 0.05)
  expect_equal(share[["wavering"]], 0.25, tolerance = 0.10)
  expect_gt(share[["committed"]], share[["wavering"]])
  # Both kinds contain both patterns: the difference is the mix, not the
  # presence of a pattern, which is the two-level quantity the package fits.
  expect_true(all(share > 0) && all(share < 1))

  # "Pooled across the sample the engaged share describes neither kind."
  pooled <- mean(engaged)
  expect_gt(pooled - share[["wavering"]], 0.2)
  expect_gt(share[["committed"]] - pooled, 0.2)

  # Student by student, the same claim again: the pooled figure describes few
  # individual students, and the shares are two clusters rather than a
  # continuum, so a two-class group model has something to find.
  by_student <- tapply(engaged, course_engagement$student, mean)
  expect_lt(min(by_student), 0.2)
  expect_gt(max(by_student), 0.9)
  expect_gt(sum(by_student > 0.5), 20L)
  expect_gt(sum(by_student < 0.5), 20L)
})

test_that("a fit recovers the measurement model the data were generated from", {
  fit <- .course_fit()
  profiles <- as.data.frame(fit)
  # Profile labels are arbitrary, so the comparison is against the sorted
  # generating means rather than against a label. From data-raw the click
  # measures are drawn at disengaged (2.90, 2.70, 2.40, 2.10) and engaged
  # (3.90, 3.60, 3.90, 3.20). `attendance` is excluded because it is not drawn
  # from a mean of its own: it is built from the click measures and then
  # blurred, which shifts both of its profile means away from the constant the
  # generator starts it at. Its own contract is the local dependence asserted
  # below.
  clicks <- subset(profiles, indicator != "attendance")
  generated <- c(2.90, 2.70, 2.40, 2.10, 3.90, 3.60, 3.90, 3.20)
  expect_equal(sort(clicks$mean), sort(generated), tolerance = 0.05)
  expect_identical(fit$n_profiles, 2L)
  expect_true(fit$converged)
  # And the enrolments land where they were generated from. `engagement` is
  # the truth the model never saw; 98% agreement is what this fixture reaches.
  expect_gt(.label_free_agreement(fit$subject_profiles,
                                  course_engagement$engagement), 0.97)
})

test_that("attendance is locally dependent on the clicks, deliberately", {
  # `attendance` is the count of days the student was active, so it is not an
  # independent measurement: a student is recorded present BECAUSE they
  # clicked. It is built from the click measures and then blurred, which leaves
  # a correlation that the shared profile alone does not account for. That is
  # the whole reason the dataset can show `bivariate_residuals()` finding a
  # real violation, so it is pinned rather than assumed.
  clicks <- c("browse", "lectures", "forum_read", "forum_post")
  within_profile <- function(level) {
    rows <- subset(course_engagement, engagement == level)
    vapply(rows[clicks], function(v) stats::cor(rows$attendance, v), numeric(1L))
  }
  residual_correlation <- c(within_profile("disengaged"),
                            within_profile("engaged"))
  # Conditioning on the generating profile does not remove it: every click
  # measure still moves with attendance inside each profile.
  expect_true(all(residual_correlation > 0.1))
  expect_gt(max(residual_correlation), 0.2)
  # And pooled, where the shared profile contributes as well, it is stronger
  # still -- every pooled correlation sits above every within-profile one.
  pooled <- vapply(course_engagement[clicks],
                   function(v) stats::cor(course_engagement$attendance, v),
                   numeric(1L))
  expect_true(all(pooled > 0.5))
  expect_gt(min(pooled), max(residual_correlation))
})

test_that("the dependence is exactly what bivariate_residuals reports", {
  fit <- .course_fit()
  residuals <- bivariate_residuals(fit)
  with_attendance <- subset(residuals, indicator_1 == "attendance" |
                              indicator_2 == "attendance")
  without <- subset(residuals, indicator_1 != "attendance" &
                      indicator_2 != "attendance")
  expect_identical(nrow(with_attendance), 8L)
  expect_identical(nrow(without), 12L)
  # Every pair that involves the derived indicator has a larger residual than
  # every pair that does not. A local-independence violation the generator put
  # there on purpose is a violation the verb has to be able to find.
  expect_gt(min(abs(with_attendance$residual)), max(abs(without$residual)))
  expect_lt(min(with_attendance$p_value), 1e-6)
})

test_that("the truth joins onto the assignments without a name clash", {
  # The documented workflow of the help page's own example. It only works
  # because no column of the data is called `profile`; a dataset that carried
  # one would make `assignments()` refuse, and the reader would have to drop
  # the column by hand before comparing the estimate against the truth.
  fit <- .course_fit()
  joined <- assignments(fit, data = course_engagement)
  expect_identical(nrow(joined), 1422L)
  expect_true(all(c("profile", "group_class", "engagement", "student_type") %in%
                    names(joined)))
  # Every column of the supplied frame comes back untouched beside the fit's.
  expect_identical(joined$engagement, course_engagement$engagement)
  expect_identical(joined$student_type, course_engagement$student_type)
  expect_identical(joined$profile, fit$subject_profiles)
  recovery <- xtabs(~ profile + engagement, data = joined)
  expect_identical(sum(recovery), 1422L)
})

test_that("a fit recovers the two-level structure, not just the profiles", {
  fit <- .course_fit()
  prevalence <- as.data.frame(fit, what = "profile_probabilities")
  # The whole claim of the package: the two group classes have opposite mixes,
  # so the lines cross. A pooled model would report one number near the middle.
  by_class <- tapply(prevalence$probability,
                     list(prevalence$group_class, prevalence$profile), identity)
  expect_gt(max(by_class[1L, ]), 0.65)
  expect_gt(max(by_class[2L, ]), 0.65)
  expect_false(identical(which.max(by_class[1L, ]), which.max(by_class[2L, ])))
  # `student_type` is the truth behind those group classes, and the students
  # sort into them the way they were generated. Compared row by row, because
  # the fit's own group order is its business and not something to reproduce
  # here.
  assigned <- assignments(fit, data = course_engagement)
  expect_gt(.label_free_agreement(assigned$group_class, assigned$student_type),
            0.8)
})

test_that("the sequence recovers the persistence it was generated with", {
  skip_on_cran()
  moves <- fit_transitions(course_engagement, vars = .activity, id = "student",
                           time = "sequence", n_profiles = 2, n_starts = 4,
                           seed = 1)
  staying <- transitions(moves, stable = TRUE)$probability
  expect_length(staying, 2L)
  # Persistence was generated per student kind -- committed (0.60, 0.90),
  # wavering (0.90, 0.70) -- so there is no single generating number to compare
  # against. The pooled persistence of the shipped `engagement` column is that
  # mixture, and it is what a one-class transition model can recover. Both are
  # sorted because the fitted labels are arbitrary.
  expect_equal(sort(staying), sort(.observed_persistence(course_engagement)),
               tolerance = 0.05)
})

test_that("previous_grade predicts engagement, as the data are documented to", {
  skip_on_cran()
  fit <- .course_fit()
  # Which fitted profile is the engaged one is a property of this seed, and the
  # sign of the coefficient below is only readable once it is pinned.
  profiles <- as.data.frame(fit)
  average <- tapply(profiles$mean, profiles$profile, mean)
  expect_gt(average[["1"]], average[["2"]])

  slope <- subset(as.data.frame(r3step(fit, course_engagement,
                                       covariates = "previous_grade")),
                  term == "previous_grade")
  expect_identical(nrow(slope), 1L)
  expect_identical(slope$outcome, "class_1")
  # A student who did well in the previous course is more likely to be engaged
  # in this one. Reported with its interval, not its p-value alone.
  expect_gt(slope$estimate, 0)
  expect_gt(slope$conf_low, 0)
  expect_lt(slope$p_value, 1e-6)
})
