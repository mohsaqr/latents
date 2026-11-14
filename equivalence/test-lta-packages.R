# lta() against LMest, seqHMM, Bartolucci et al. (2017) on RLMSdat and
# Pennoni et al. (2024) on data_market_sim (lta-packages/). Tolerances are
# fixed in lta-packages/compare-functions.R.
# Disagreements are allowed only where lta-packages/known-differences.csv explains
# them, and every entry there must still match one.

invisible(lapply(c("LMest", "seqHMM", "TraMineR", "withr"), function(p) skip_if_not_installed(p)))

lta_packages_dir <- file.path("lta-packages")
source(file.path(lta_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "lta-packages", compare_all_lta_packages, n_rows = 52L,
  groups = c("R packages", "published analyses"),
  coverage = function(result) testthat::expect_gte(sum(result$comparison$agrees), 48L))
