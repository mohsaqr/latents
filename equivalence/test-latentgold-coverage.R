# Growth mixtures, lca(), lpa(), lta() and class enumeration against Latent
# GOLD 6.1 (latentgold-coverage/): latent class growth and random-effect growth
# mixtures (exact, as factor models), binary, nominal and covariate latent
# class models, VVI/EEI/VVV profile models, Gaussian, binary and
# covariate-driven transition models, and lca()/lpa() enumeration at one to
# four classes. Latent GOLD's listings are retained in returned/, so this runs
# offline. Every estimate and standard error Latent GOLD prints is compared
# after expressing the latents fit in Latent GOLD's coding, with delta-method
# standard errors, plus every row's posteriors.
#
# Tolerances, fixed before the comparison: half a unit of Latent GOLD's fourth
# printed decimal (5e-5) plus 2e-6; the log likelihood to 6e-5, BIC to 1.2e-4;
# parameter counts exactly; posteriors to 2e-5. Disagreements are allowed only
# where known-differences.csv explains them, and every entry there must still
# match one.

coverage_dir <- file.path("latentgold-coverage")
source(file.path(coverage_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("the coverage kit is the one Latent GOLD was run on", {
  inputs <- utils::read.table(file.path(coverage_dir, "returned", "INPUTS.md5"),
                              col.names = c("md5", "file"))
  outputs <- utils::read.table(file.path(coverage_dir, "returned", "OUTPUTS.md5"),
                               col.names = c("md5", "file"))
  expect_identical(unname(tools::md5sum(file.path(coverage_dir, "kit", inputs$file))), inputs$md5)
  expect_identical(unname(tools::md5sum(file.path(coverage_dir, "returned", outputs$file))),
                   outputs$md5)
})

test_that("growth, latent class, profile, transition and enumeration fits reproduce Latent GOLD", {
  cases <- readRDS(file.path(coverage_dir, "kit", "cases.rds"))
  expect_length(cases, 22L)
  known <- utils::read.csv(file.path(coverage_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all(cases, file.path(coverage_dir, "returned"))
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 440L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$quantity, collapse = "; "))
  expect_identical(nrow(registered$stale), 0L)
  expect_gte(sum(comparison$agrees), 434L)
})
