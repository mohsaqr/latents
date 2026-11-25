# The contract `?engagement_momentum` states: its shape, columns, scale and
# the structure it was generated from, which a multilevel latent profile model
# must recover.

.momentum_vars <- c("browse", "lectures", "forum_read", "forum_post", "active_days")

test_that("engagement_momentum has the documented shape and columns", {
  expect_s3_class(engagement_momentum, "data.frame")
  expect_identical(dim(engagement_momentum), c(1310L, 12L))
  expect_identical(names(engagement_momentum),
                   c("student", "course", "sequence", .momentum_vars, "grade",
                     "previous_grade", "state", "trajectory"))
  expect_false(anyNA(engagement_momentum))
  expect_identical(length(unique(engagement_momentum$student)), 106L)
  expect_identical(levels(engagement_momentum$state),
                   c("active", "average", "disengaged"))
  expect_identical(levels(engagement_momentum$trajectory),
                   c("engaged", "fluctuating", "mostly_disengaged"))
})

test_that("every student has 10 to 15 courses in sequence and one trajectory", {
  by_student <- split(engagement_momentum, engagement_momentum$student)
  lengths_ok <- vapply(by_student, function(d) {
    nrow(d) >= 10L && nrow(d) <= 15L && identical(d$sequence, seq_len(nrow(d)))
  }, logical(1))
  expect_true(all(lengths_ok))
  one_trajectory <- vapply(by_student, function(d) {
    length(unique(d$trajectory)) == 1L
  }, logical(1))
  expect_true(all(one_trajectory))
})

test_that("each indicator is standardized within every course", {
  standardized <- vapply(.momentum_vars, function(v) {
    by_course <- split(engagement_momentum[[v]], engagement_momentum$course)
    all(vapply(by_course, function(z) {
      abs(mean(z)) < 0.02 && (length(z) < 2L || abs(stats::sd(z) - 1) < 0.02)
    }, logical(1)))
  }, logical(1))
  expect_true(all(standardized))
})

test_that("a three-profile, three-class model recovers the states and trajectories", {
  fit <- multilpa(engagement_momentum, .momentum_vars, id = "student",
                  n_profiles = 3, n_group_classes = 3, n_starts = 3, seed = 1)
  recovery <- get_results(fit, "assignments", data = engagement_momentum,
                          truth = c("state", "trajectory"))
  agreement <- function(level) {
    rows <- subset(recovery, assignment == level)
    sum(tapply(rows$n, rows$class, max)) / sum(rows$n)
  }
  expect_gt(agreement("profile"), 0.85)
  expect_gt(agreement("group_class"), 0.85)
  # The generating model has conditionally independent indicators, so no
  # residual association should survive the multiplicity adjustment.
  residuals <- get_results(fit, "residuals")
  expect_true(all(residuals$p_adjusted > 0.05))
  # Higher grades go with more active states.
  outcome <- three_step(fit, data = engagement_momentum, outcome = "grade")
  expect_gt(diff(range(as.data.frame(outcome)$estimate)), 1)
})
