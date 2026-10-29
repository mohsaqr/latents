# Growth mixtures in mixture_regression() against genuine Mplus 9 Demo runs
# (mplus-growth/, copied from the JStats fixtures): latent class growth at one
# to three classes on two datasets and growth mixtures with class-specific
# random intercept and slope covariances. Mplus prints three decimals; the
# tolerance is half a unit plus 1e-6 of the value. Disagreements are allowed
# only where mplus-growth/known-differences.csv explains them.

mplus_growth_dir <- file.path("mplus-growth")
source(file.path(mplus_growth_dir, "compare-functions.R"), local = TRUE)

test_that("growth mixtures reproduce genuine Mplus runs", {
  known <- utils::read.csv(file.path(mplus_growth_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_mplus_growth(mplus_growth_dir)
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    any(known$case == comparison$case[[i]] & known$reference == comparison$reference[[i]] &
          vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
  }, logical(1))
  unexplained <- subset(comparison, !agrees & !explained)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$quantity, collapse = "; "))
})
