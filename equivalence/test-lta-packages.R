# lta() against LMest, seqHMM, Bartolucci et al. (2017) on RLMSdat and
# Pennoni et al. (2024) on data_market_sim (lta-packages/). Tolerances are
# fixed in lta-packages/compare-functions.R.
# Disagreements are allowed only where lta-packages/known-differences.csv explains
# them, and every entry there must still match one.

invisible(lapply(c("LMest", "seqHMM", "TraMineR", "withr"), function(p) skip_if_not_installed(p)))

lta_packages_dir <- file.path("lta-packages")
source(file.path(lta_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("lta-packages agrees with its references", {
  known <- utils::read.csv(file.path(lta_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_lta_packages()
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 52L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  expect_identical(nrow(registered$stale), 0L)
  expect_gte(sum(comparison$agrees), 48L)
})
