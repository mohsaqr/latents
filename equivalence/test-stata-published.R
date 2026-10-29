# lpa() and lca() against genuine Stata 19 output printed in its manual
# (Stata [SEM] Examples 52g and 50g) and Masyn (2013) Table 25.12, including
# its adjusted LMR p-values (stata-published/). Tolerances are fixed in
# stata-published/compare-functions.R. Disagreements are allowed only where
# stata-published/known-differences.csv explains them.

stata_dir <- file.path("stata-published")
source(file.path(stata_dir, "compare-functions.R"), local = TRUE)

test_that("lpa() and lca() reproduce Stata's printed output and Masyn's table", {
  known <- utils::read.csv(file.path(stata_dir, "known-differences.csv"), stringsAsFactors = FALSE)
  comparison <- compare_all_stata(stata_dir)
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    any(known$case == comparison$case[[i]] & known$reference == comparison$reference[[i]] &
          vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
  }, logical(1))
  unexplained <- subset(comparison, !agrees & !explained)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
})
