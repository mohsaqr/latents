# pool_imputations() against mitml and mitools (pool-packages/). Tolerances are fixed in pool-packages/compare-functions.R.
# Disagreements are allowed only where pool-packages/known-differences.csv explains
# them, and every entry there must still match one.

invisible(lapply(c("mitml", "mitools", "mice", "withr"), function(p) skip_if_not_installed(p)))

pool_packages_dir <- file.path("pool-packages")
source(file.path(pool_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("pool-packages agrees with its references", {
  known <- utils::read.csv(file.path(pool_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_pool()
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 36L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  expect_identical(nrow(registered$stale), 0L)
})
