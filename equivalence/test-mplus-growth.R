# Growth mixtures in mixture_regression() against genuine Mplus 9 Demo runs
# (mplus-growth/, copied from the JStats fixtures): latent class growth at one
# to three classes on two datasets and growth mixtures with class-specific
# random intercept and slope covariances. Mplus prints three decimals; the
# tolerance is half a unit plus 1e-6 of the value. Disagreements are allowed
# only where mplus-growth/known-differences.csv explains them.

mplus_growth_dir <- file.path("mplus-growth")
source(file.path(mplus_growth_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "mplus-growth", function() compare_all_mplus_growth(mplus_growth_dir), n_rows = 212L,
  groups = "Mplus")
