# Ordinal and count indicators against Latent GOLD 6.1 (latentgold-ordinal/):
# seven simulated cases -- ordinal (adjacent-category logit with equally
# spaced scores), Poisson and negative-binomial counts with class-varying or
# shared dispersion, mixed with continuous indicators, one two-level -- fitted
# by multilpa() and compared with Latent GOLD's retained listings.
#
# Tolerances, as the verification matrix reads this kit: the log likelihood
# to 6e-5 and every class size, category probability and mean to 1e-4
# (Latent GOLD prints four decimals).

ordinal_dir <- file.path("latentgold-ordinal")
source(file.path(ordinal_dir, "compare-functions.R"), local = TRUE)
source("package-kit-runner.R", local = TRUE)

test_that("ordinal and count indicators reproduce Latent GOLD", {
  snapshot <- matrix_capture_source("..")
  comparison <- compare_ordinal_kit(ordinal_dir)
  kit_save_comparison(ordinal_dir, comparison, snapshot)
  expect_identical(nrow(comparison), 7L)
  expect_true(all(abs(comparison$ll_difference) <= 6e-5),
              info = paste(comparison$case, signif(comparison$ll_difference, 3), collapse = "; "))
  expect_true(all(comparison$max_parameter_difference <= 1e-4),
              info = paste(comparison$case, signif(comparison$max_parameter_difference, 3),
                           collapse = "; "))
})
