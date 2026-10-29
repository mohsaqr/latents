# Comparisons with independent references moved from tests/testthat/test-tidy-diagnostics.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("Cramer's V matches an independent computation", {
  # Calibration against stats::chisq.test, not against our own arithmetic.
  counts <- matrix(c(30, 10, 5, 55), 2L, 2L)
  reference <- quietly(stats::chisq.test(counts, correct = FALSE))
  expect_equal(.multilpa_cramers_v(counts),
               sqrt(unname(reference$statistic) / sum(counts)))
  # For a two-by-two table Cramer's V is the absolute indicator correlation.
  long <- data.frame(a = rep(c(0, 0, 1, 1), c(30, 5, 10, 55)),
                     b = rep(c(0, 1, 0, 1), c(30, 5, 10, 55)))
  expect_equal(.multilpa_cramers_v(counts), abs(stats::cor(long$a, long$b)))
  three <- matrix(c(20, 5, 3, 4, 18, 25), 3L, 2L)
  reference_three <- quietly(stats::chisq.test(three, correct = FALSE))
  expect_equal(.multilpa_cramers_v(three),
               sqrt(unname(reference_three$statistic) / (sum(three) * 1L)))
  # A product of its own margins carries no association at all.
  expect_equal(.multilpa_cramers_v(outer(c(30, 70), c(40, 60)) / 100), 0)
  # Undefined rather than zero where there is no table to speak of.
  expect_true(is.na(.multilpa_cramers_v(matrix(0, 2L, 2L))))
  expect_true(is.na(.multilpa_cramers_v(matrix(c(3, 4), 2L, 1L))))
})
