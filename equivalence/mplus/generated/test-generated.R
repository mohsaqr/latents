# Regression checks against independently fitted, retained Mplus9 outputs.
artifact_dir <- if (file.exists("comparison.rds")) "." else
  file.path("equivalence", "mplus", "generated")
comparisons <- readRDS(file.path(artifact_dir, "comparison.rds"))
testthat::test_that("both synthetic multilevel models reproduce genuine Mplus estimates", {
  testthat::expect_identical(names(comparisons), c("varying", "equal"))
  invisible(lapply(comparisons, function(comparison) {
    stopifnot(is.list(comparison), is.numeric(comparison$differences))
    differences <- comparison$differences
    testthat::expect_true(comparison$fit$converged)
    testthat::expect_true(all(differences[c("means", "variances",
      "profile_probabilities", "group_probabilities")] < 1e-6))
    testthat::expect_true(all(differences[c("subject_posteriors", "group_posteriors")] < 1e-6))
    testthat::expect_true(all(differences[c("log_likelihood", "aic", "bic_individual")] < 1e-4))
    output <- readLines(file.path(artifact_dir, paste0(comparison$variance_model, ".out")))
    testthat::expect_true(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", output)))
    testthat::expect_true(any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", output)))
    testthat::expect_false(any(grepl("*** ERROR", output, fixed = TRUE)))
  }))
})
