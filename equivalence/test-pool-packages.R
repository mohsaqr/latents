# pool_imputations() against mitml and mitools (pool-packages/). Tolerances are fixed in pool-packages/compare-functions.R.
# Disagreements are allowed only where pool-packages/known-differences.csv explains
# them, and every entry there must still match one.

invisible(lapply(c("mitml", "mitools", "mice", "withr"), function(p) skip_if_not_installed(p)))

pool_packages_dir <- file.path("pool-packages")
source(file.path(pool_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "pool-packages", compare_all_pool, n_rows = 36L, groups = "R packages", kind = "R package")
