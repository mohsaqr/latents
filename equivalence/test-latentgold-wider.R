# Mixture regression and sampling weights against Latent GOLD 6.1
# (latentgold-wider/): Gaussian, Poisson and logistic mixture regressions with
# a membership covariate, their model-based and robust standard errors, and
# sampling weights in a weighted mixture regression, lpa(), lta() and a
# two-level multilpa(). Latent GOLD's listings are retained in returned/, so
# this runs offline. Every estimate and standard error Latent GOLD prints is
# compared, after expressing the latents fit in Latent GOLD's coding
# (reference-coded logits, mean contrasts, variances) with delta-method
# standard errors.
#
# Tolerances, fixed before the comparison: half a unit of Latent GOLD's
# fourth printed decimal (5e-5) plus 2e-6 for the two programs' convergence;
# the log likelihood to 6e-5; parameter counts exactly.

wider_dir <- file.path("latentgold-wider")
source(file.path(wider_dir, "compare-functions.R"), local = TRUE)

test_that("the kit is the one Latent GOLD was run on", {
  recorded <- utils::read.table(file.path(wider_dir, "returned", "INPUTS.md5"),
                                col.names = c("md5", "file"))
  files <- file.path(wider_dir, "kit", recorded$file)
  expect_identical(unname(tools::md5sum(files)), recorded$md5)
})

test_that("mixture regressions and weighted fits reproduce Latent GOLD", {
  cases <- readRDS(file.path(wider_dir, "kit", "cases.rds"))
  expect_identical(names(cases), c("r01", "r02", "r03", "r04", "r05", "w01", "w02", "w03"))
  comparison <- compare_all(cases, file.path(wider_dir, "returned"))
  expect_identical(nrow(comparison), 172L)
  disagreements <- subset(comparison, !agrees)
  expect_identical(nrow(disagreements), 0L,
                   info = paste(disagreements$case, disagreements$quantity, collapse = "; "))
})
