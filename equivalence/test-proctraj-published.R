# Latent class growth models (mixture_regression() with one class per person
# and a shared residual standard deviation) against SAS PROC TRAJ output
# printed in Andruff et al. (2009), Tutorials in Quantitative Methods for
# Psychology 5(1), 11-24 (proctraj-published/). Tolerances and the standard-
# error factor are fixed in proctraj-published/compare-functions.R.
# Disagreements are allowed only where proctraj-published/known-differences.csv
# explains them.

proctraj_dir <- file.path("proctraj-published")
source(file.path(proctraj_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "proctraj-published", function() compare_proctraj(proctraj_dir), n_rows = 44L, groups = "SAS")
