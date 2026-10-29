# lta() against LMest, seqHMM, Bartolucci et al. (2017) on RLMSdat and
# Pennoni et al. (2024) on data_market_sim (lta-packages/). Tolerances are
# fixed in lta-packages/compare-functions.R.
# Disagreements are allowed only where lta-packages/known-differences.csv explains
# them, and every entry there must still match one.

invisible(lapply(c("LMest", "seqHMM", "TraMineR", "withr"), function(p) skip_if_not_installed(p)))

lta_packages_dir <- file.path("lta-packages")
source(file.path(lta_packages_dir, "compare-functions.R"), local = TRUE)

test_that("lta-packages agrees with its references", {
  known <- utils::read.csv(file.path(lta_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_lta_packages()
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
  expect_gte(sum(comparison$agrees), 48L)
})
