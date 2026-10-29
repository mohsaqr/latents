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

test_that("lca() reproduces five independent latent class implementations", {
  known <- utils::read.csv(file.path(packages_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_lca()
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    any(known$case == comparison$case[[i]] & known$reference == comparison$reference[[i]] &
          vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
  }, logical(1))
  unexplained <- subset(comparison, !agrees & !explained)
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  expect_true(all(vapply(seq_len(nrow(known)), function(r) {
    any(!comparison$agrees & comparison$case == known$case[[r]] &
          comparison$reference == known$reference[[r]] &
          grepl(known$quantity[[r]], comparison$quantity))
  }, logical(1))))
  expect_setequal(unique(comparison$reference),
                  c("glca", "randomLCA", "BayesLCA", "multilevLCA", "depmixS4"))
  expect_gte(sum(comparison$agrees), 203L)
})
