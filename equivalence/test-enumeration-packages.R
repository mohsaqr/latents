# enumerate_lpa() against mclust::mclustBIC(), enumerate_lca() against poLCA,
# and bootstrap_lrt() against mclust::mclustBootstrapLRT() and mclust's
# vignette (enumeration-packages/). Tolerances are fixed in
# enumeration-packages/compare-functions.R. Disagreements are allowed only
# where enumeration-packages/known-differences.csv explains them.

invisible(lapply(c("mclust", "poLCA", "withr"), function(p) skip_if_not_installed(p)))

enumeration_packages_dir <- file.path("enumeration-packages")
source(file.path(enumeration_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "enumeration-packages", compare_all_enumeration, n_rows = 168L, groups = "R packages",
  coverage = function(result) testthat::expect_gte(sum(result$comparison$agrees), 164L))
