# pool_imputations() against mitml and mitools (pool-packages/). Tolerances are fixed in pool-packages/compare-functions.R.
# Disagreements are allowed only where pool-packages/known-differences.csv explains
# them, and every entry there must still match one.

invisible(lapply(c("mitml", "mitools", "mice", "withr"), function(p) skip_if_not_installed(p)))

pool_packages_dir <- file.path("pool-packages")
source(file.path(pool_packages_dir, "compare-functions.R"), local = TRUE)

test_that("pool-packages agrees with its references", {
  known <- utils::read.csv(file.path(pool_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_pool()
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    any(known$case == comparison$case[[i]] & known$reference == comparison$reference[[i]] &
          vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
  }, logical(1))
  unexplained <- subset(comparison, !agrees & !explained)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  expect_true(all(vapply(seq_len(nrow(known)), function(r) {
    any(!comparison$agrees & comparison$case == known$case[[r]] &
          comparison$reference == known$reference[[r]] &
          grepl(known$quantity[[r]], comparison$quantity))
  }, logical(1))))
})
