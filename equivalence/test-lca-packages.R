# lca() against glca, randomLCA, BayesLCA (EM), multilevLCA and depmixS4 on
# eight datasets (lca-packages/): maximized log likelihood, class shares,
# item response probabilities and, for glca and depmixS4, every posterior.
# Tolerances, fixed before the comparison: the log likelihood to 1e-6, the
# rest to 1e-4. Disagreements are allowed only where
# lca-packages/known-differences.csv explains them (multilevLCA's probability
# bounds and single start), and every entry there must still match one.

skip_if_not_installed("glca")
skip_if_not_installed("randomLCA")
skip_if_not_installed("BayesLCA")
skip_if_not_installed("multilevLCA")
skip_if_not_installed("depmixS4")
skip_if_not_installed("poLCA")
skip_if_not_installed("withr")

packages_dir <- file.path("lca-packages")
source(file.path(packages_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "lca-packages", compare_all_lca, n_rows = 211L,
  groups = c("R packages", "published analyses"),
  coverage = function(result) {
    comparison <- result$comparison
    testthat::expect_setequal(unique(comparison$reference[comparison$kind == "R package"]),
                              c("glca", "randomLCA", "BayesLCA", "multilevLCA", "depmixS4"))
    testthat::expect_setequal(unique(comparison$reference[comparison$kind == "published"]),
                              "Kleinman and Horton (2011), SAS and R blog, Example 8.21")
    # Rows explained by a "sometimes" entry (a platform-dependent reference)
    # count with the agreeing rows; the floor still bounds what the register waives.
    optional <- result$explained %in% which(result$known$occurs == "sometimes")
    testthat::expect_gte(sum(comparison$agrees | optional), 203L)
  })
