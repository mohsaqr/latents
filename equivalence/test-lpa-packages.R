# lpa() against Rmixmod (VVI, EEI, VVV, EEE), EMCluster (VVV) and mixtools
# (VVV, EEE) on six datasets (lpa-packages/), compared at each reference's
# solution: likelihood, shares, means, (co)variances and every posterior.
# Tolerances, fixed before the comparison: the log likelihood to 1e-6 plus
# 1e-10 of its size; parameters to 1e-4 of max(1, |value|); posteriors to
# 1e-4. Disagreements are allowed only where lpa-packages/known-differences.csv
# explains them, and every entry there must still match one.

skip_if_not_installed("Rmixmod")
skip_if_not_installed("EMCluster")
skip_if_not_installed("mixtools")
skip_if_not_installed("mclust")
skip_if_not_installed("MASS")
skip_if_not_installed("withr")

lpa_packages_dir <- file.path("lpa-packages")
source(file.path(lpa_packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("lpa() reproduces independent Gaussian mixture implementations", {
  known <- utils::read.csv(file.path(lpa_packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_lpa()
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 156L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  expect_identical(nrow(registered$stale), 0L)
  expect_gte(sum(comparison$agrees), 155L)
})
