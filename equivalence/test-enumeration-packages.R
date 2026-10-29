# enumerate_lpa() against mclust::mclustBIC(), enumerate_lca() against poLCA,
# and bootstrap_lrt() against mclust::mclustBootstrapLRT() and mclust's
# vignette (enumeration-packages/). Tolerances are fixed in
# enumeration-packages/compare-functions.R. Disagreements are allowed only
# where enumeration-packages/known-differences.csv explains them.

invisible(lapply(c("mclust", "poLCA", "withr"), function(p) skip_if_not_installed(p)))

enumeration_packages_dir <- file.path("enumeration-packages")
source(file.path(enumeration_packages_dir, "compare-functions.R"), local = TRUE)

test_that("enumeration grids and the bootstrap LRT reproduce mclust and poLCA", {
  known <- utils::read.csv(file.path(enumeration_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_enumeration()
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    any(known$case == comparison$case[[i]] & known$reference == comparison$reference[[i]] &
          vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
  }, logical(1))
  unexplained <- subset(comparison, !agrees & !explained)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  # Every register entry still explains a disagreement.
  expect_true(all(vapply(seq_len(nrow(known)), function(r) {
    any(!comparison$agrees & comparison$case == known$case[[r]] &
          comparison$reference == known$reference[[r]] &
          grepl(known$quantity[[r]], comparison$quantity))
  }, logical(1))))
  expect_gte(sum(comparison$agrees), 164L)
})
