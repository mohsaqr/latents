# enumerate_lpa() against mclust::mclustBIC(), enumerate_lca() against poLCA,
# and bootstrap_lrt() against mclust::mclustBootstrapLRT() and mclust's
# vignette (enumeration-packages/). Tolerances are fixed in
# enumeration-packages/compare-functions.R. Disagreements are allowed only
# where enumeration-packages/known-differences.csv explains them.

invisible(lapply(c("mclust", "poLCA", "withr"), function(p) skip_if_not_installed(p)))

enumeration_packages_dir <- file.path("enumeration-packages")
source(file.path(enumeration_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("enumeration grids and the bootstrap LRT reproduce mclust and poLCA", {
  known <- utils::read.csv(file.path(enumeration_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_enumeration()
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 168L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  # Every register entry still explains a disagreement.
  expect_identical(nrow(registered$stale), 0L)
  expect_gte(sum(comparison$agrees), 164L)
})
