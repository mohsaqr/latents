# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-tidy-diagnostics.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("Cramer's V matches an independent computation", {
  note <- "Package side: the internal Cramer's V (.multilpa_cramers_v) behind the association columns of get_residuals()."
  # Calibration against stats::chisq.test, not against our own arithmetic.
  # testthat's default tolerance (sqrt of the machine epsilon, relative) of the
  # original test, kept.
  counts <- matrix(c(30, 10, 5, 55), 2L, 2L)
  reference <- quietly(stats::chisq.test(counts, correct = FALSE))
  expect_reference(latents:::.multilpa_cramers_v(counts),
                   sqrt(unname(reference$statistic) / sum(counts)),
                   reference = "stats", covers = "get_residuals()",
                   quantity = "Cramer's V from chisq.test(), two-by-two table",
                   tolerance = sqrt(.Machine$double.eps), measure = "relative", note = note)
  # For a two-by-two table Cramer's V is the absolute indicator correlation.
  long <- data.frame(a = rep(c(0, 0, 1, 1), c(30, 5, 10, 55)),
                     b = rep(c(0, 1, 0, 1), c(30, 5, 10, 55)))
  expect_reference(latents:::.multilpa_cramers_v(counts), abs(stats::cor(long$a, long$b)),
                   reference = "stats", covers = "get_residuals()",
                   quantity = "Cramer's V as the absolute correlation from cor(), two-by-two table",
                   tolerance = sqrt(.Machine$double.eps), measure = "relative", note = note)
  three <- matrix(c(20, 5, 3, 4, 18, 25), 3L, 2L)
  reference_three <- quietly(stats::chisq.test(three, correct = FALSE))
  expect_reference(latents:::.multilpa_cramers_v(three),
                   sqrt(unname(reference_three$statistic) / (sum(three) * 1L)),
                   reference = "stats", covers = "get_residuals()",
                   quantity = "Cramer's V from chisq.test(), three-by-two table",
                   tolerance = sqrt(.Machine$double.eps), measure = "relative", note = note)
  # A product of its own margins carries no association at all.
  expect_equal(latents:::.multilpa_cramers_v(outer(c(30, 70), c(40, 60)) / 100), 0)
  # Undefined rather than zero where there is no table to speak of.
  expect_true(is.na(latents:::.multilpa_cramers_v(matrix(0, 2L, 2L))))
  expect_true(is.na(latents:::.multilpa_cramers_v(matrix(c(3, 4), 2L, 1L))))
})
