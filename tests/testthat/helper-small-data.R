# Tests check mechanics, not estimates, so none of them fits the full bundled
# data: 20 students (274 enrolments) are enough, and keep the suite fast.
# test-data.R is the exception, because it tests the dataset itself.
engagement_small <- subset(course_engagement, student <= 20)
row.names(engagement_small) <- NULL

# The emEM starts table: the screened pool (five times `n_starts` by default)
# with `n_starts` of them run to the end.
expect_emem_starts <- function(starts, n_starts) {
  expect_identical(nrow(starts), 5L * as.integer(n_starts))
  expect_identical(sum(starts$stage == "final"), as.integer(n_starts))
}
