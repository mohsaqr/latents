test_that("the Lo-Mendell-Rubin adjustment reproduces genuine Mplus TECH11", {
  reference <- readRDS(equivalence_fixture("mplus", "lmr-tech11.rds"))
  expect_identical(nrow(reference), 2L)
  adjusted <- reference$statistic / (1 + 1 / (reference$df * log(reference$n)))
  # Mplus prints TECH11 to three decimals; that is the comparison precision.
  expect_lt(max(abs(adjusted - reference$mplus_adjusted)), 5e-3)
})

