# lpa() and lca() against genuine Stata 19 output printed in its manual
# (Stata [SEM] Examples 52g and 50g) and Masyn (2013) Table 25.12, including
# its adjusted LMR p-values (stata-published/). Tolerances are fixed in
# stata-published/compare-functions.R. Disagreements are allowed only where
# stata-published/known-differences.csv explains them.

stata_dir <- file.path("stata-published")
source(file.path(stata_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "stata-published", function() compare_all_stata(stata_dir), n_rows = 154L,
  groups = c("Stata", "published analyses"))
