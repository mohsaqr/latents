# mixture_regression() against mixtools, MoEClust and flexmix and against
# three published analyses (regression-packages/): SAS PROC FMM Example 40.3,
# Grun and Leisch (2008) and Benaglia et al. (2009). Tolerances, fixed before
# the comparison: software, the log likelihood to 1e-6 plus 1e-10 of its size
# and the rest to 1e-4 of max(1, |value|); published values, half a unit of
# the last printed digit. Disagreements are allowed only where
# regression-packages/known-differences.csv explains them.

skip_if_not_installed("mixtools")
skip_if_not_installed("MoEClust")
skip_if_not_installed("flexmix")
skip_if_not_installed("withr")

regression_packages_dir <- file.path("regression-packages")
source(file.path(regression_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "regression-packages", compare_all_regression, n_rows = 76L,
  groups = c("R packages", "published analyses"),
  coverage = function(result) {
    published <- subset(result$comparison, kind == "published")
    testthat::expect_identical(nrow(published), 20L)
    testthat::expect_true(all(published$agrees))
    # Rows explained by a "sometimes" entry (a platform-dependent reference)
    # count with the agreeing rows; the floor still bounds what the register waives.
    optional <- result$explained %in% which(result$known$occurs == "sometimes")
    testthat::expect_gte(sum(result$comparison$agrees | optional), 71L)
  })
