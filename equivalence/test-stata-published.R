# lpa() and lca() against genuine Stata 19 output printed in its manual
# (Stata [SEM] Examples 52g and 50g) and Masyn (2013) Table 25.12, including
# its adjusted LMR p-values (stata-published/). Tolerances are fixed in
# stata-published/compare-functions.R. Disagreements are allowed only where
# stata-published/known-differences.csv explains them.

stata_dir <- file.path("stata-published")
source(file.path(stata_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("lpa() and lca() reproduce Stata's printed output and Masyn's table", {
  known <- utils::read.csv(file.path(stata_dir, "known-differences.csv"), stringsAsFactors = FALSE)
  comparison <- compare_all_stata(stata_dir)
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 154L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(registered$stale), 0L)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
})
