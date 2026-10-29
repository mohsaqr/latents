# mixture_regression() against mixtools, MoEClust and flexmix and against
# three published analyses (regression-packages/): SAS PROC FMM Example 40.3,
# Grun and Leisch (2008) and Benaglia et al. (2009). Tolerances, fixed before
# the comparison: software, the log likelihood to 1e-6 plus 1e-10 of its size
# and the rest to 1e-4 of max(1, |value|); published values, half a unit of
# the last printed digit. Disagreements are allowed only where
# regression-packages/known-differences.csv explains them.

skip_if_not_installed("mixtools")
skip_if_not_installed("MoEClust")
skip_if_not_installed("flexmix")
skip_if_not_installed("withr")

regression_packages_dir <- file.path("regression-packages")
source(file.path(regression_packages_dir, "compare-functions.R"), local = TRUE)

test_that("mixture_regression() reproduces three packages and three publications", {
  known <- utils::read.csv(file.path(regression_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_regression()
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    any(known$case == comparison$case[[i]] & known$reference == comparison$reference[[i]] &
          vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
  }, logical(1))
  unexplained <- subset(comparison, !agrees & !explained)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  published <- subset(comparison, kind == "published")
  expect_identical(nrow(published), 20L)
  expect_true(all(published$agrees))
  expect_gte(sum(comparison$agrees), 71L)
})
