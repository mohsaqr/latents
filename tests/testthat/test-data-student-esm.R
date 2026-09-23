# `student_esm` is the bundled categorical dataset, and these tests are the
# contract `?student_esm` states. Every documented claim is asserted here.

.leisure <- c("time_with_friends", "on_social_media", "tv_video_games",
              "listened_music", "sports", "walking", "reading",
              "part_time_job")
.affect <- c("happy", "relaxed", "worried", "exhausted")

test_that("student_esm has the documented shape and types", {
  expect_s3_class(student_esm, "data.frame")
  expect_identical(dim(student_esm), c(2582L, 15L))
  expect_identical(names(student_esm),
                   c("student", "day", "beep", .leisure, .affect))
  expect_false(anyNA(student_esm))

  expect_type(student_esm$student, "integer")
  expect_type(student_esm$day, "integer")
  expect_type(student_esm$beep, "integer")
  expect_true(all(vapply(student_esm[.leisure], \(v)
    is.factor(v) && identical(levels(v), c("no", "yes")), logical(1L))))
  expect_true(all(vapply(student_esm[.affect], \(v)
    is.integer(v) && all(v >= 1L & v <= 7L), logical(1L))))
})

test_that("student_esm nests 15 to 51 prompts in each of 100 students", {
  per_student <- as.vector(table(student_esm$student))
  expect_identical(sort(unique(student_esm$student)), seq_len(100L))
  expect_identical(range(per_student), c(15L, 51L))
  expect_identical(stats::median(per_student), 24)
  expect_identical(range(student_esm$day), c(0L, 13L))
  expect_identical(range(student_esm$beep), c(1L, 5L))
  expect_identical(anyDuplicated(student_esm[c("student", "day", "beep")]), 0L)
})

test_that("the documented floor on `worried` is in the data", {
  expect_equal(mean(student_esm$worried == 1L), 0.497, tolerance = 5e-4)
})

test_that("the first week supports a two-level latent class model", {
  skip_on_cran()
  first_week <- subset(student_esm, day <= 6)
  expect_identical(nrow(first_week), 1422L)
  lca <- multilpa(first_week, vars = .leisure, id = "student",
                  n_profiles = 2, n_group_classes = 2,
                  categorical = .leisure, n_starts = 3, seed = 1)
  model <- get_results(lca, "model")
  expect_true(model$converged)
  expect_false(model$boundary)
  responses <- get_results(lca, "responses")
  expect_true(all(responses$probability > 1e-6 &
                    responses$probability < 1 - 1e-6))
  # Response probabilities of each item sum to one within a profile.
  sums <- stats::aggregate(probability ~ profile + indicator, responses, sum)
  expect_true(all(abs(sums$probability - 1) < sqrt(.Machine$double.eps)))
})
