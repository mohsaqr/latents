# Growth mixtures in mixture_regression() against genuine Mplus 9 Demo runs
# (mplus-growth/, copied from the JStats fixtures): latent class growth at one
# to three classes on two datasets and growth mixtures with class-specific
# random intercept and slope covariances. Mplus prints three decimals; the
# tolerance is half a unit plus 1e-6 of the value. Disagreements are allowed
# only where mplus-growth/known-differences.csv explains them.

mplus_growth_dir <- file.path("mplus-growth")
source(file.path(mplus_growth_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("growth mixtures reproduce genuine Mplus runs", {
  known <- utils::read.csv(file.path(mplus_growth_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_mplus_growth(mplus_growth_dir)
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 212L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(registered$stale), 0L)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$quantity, collapse = "; "))
})
