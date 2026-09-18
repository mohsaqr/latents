# Mplus9 independently fitted these synthetic models with 100/20 starts.
# Every reference estimate and probability is parsed from actual Mplus files.
test_that("native multilevel LPA reproduces Mplus9 on two synthetic variance models", {
  invisible(lapply(c("varying", "equal"), function(variance_model) {
    stopifnot(is.character(variance_model), length(variance_model) == 1L)
    reference <- readRDS(test_path("..", "fixtures", "mplus",
      paste0("twolevel-synthetic-", variance_model, ".rds")))
    fit <- multilpa(reference$data, c("y1", "y2"), "clus", 2L, 2L,
      variance_model = variance_model, n_starts = 10L, tol = 1e-12,
      max_iter = 5000L, seed = 5739)
    # Align profile and group labels independently; neither label order is fixed.
    profile_order <- order(fit$means[, "y1"])[rank(reference$means[, 1])]
    aligned_profile <- fit$profile_probabilities[, profile_order, drop = FALSE]
    group_order <- order(aligned_profile[, 1])[rank(reference$profile_probabilities[, 1])]
    expect_true(fit$converged)
    expect_false(fit$boundary)
    expect_equal(fit$n_parameters, reference$n_parameters)
    expect_equal(unname(fit$means[profile_order, ]), reference$means, tolerance = 1e-6)
    expect_equal(unname(fit$variances[profile_order, ]), reference$variances, tolerance = 1e-6)
    expect_equal(unname(aligned_profile[group_order, ]), unname(reference$profile_probabilities),
      tolerance = 1e-6)
    expect_equal(unname(fit$group_probabilities[group_order]), reference$group_probabilities,
      tolerance = 1e-6)
    expect_equal(unname(fit$subject_posteriors[, profile_order]), reference$subject_posteriors,
      tolerance = 1e-6)
    expect_equal(unname(fit$group_posteriors[, group_order]), reference$group_posteriors,
      tolerance = 1e-6)
    expect_lt(abs(fit$log_likelihood - reference$log_likelihood), 1e-4)
    expect_lt(abs(fit$aic - reference$aic), 1e-4)
    expect_lt(abs(fit$bic_individual - reference$bic_individual), 1e-4)
    expect_match(reference$provenance$version, "Mplus VERSION 9 DEMO")
    expect_length(reference$provenance$md5, 5L)
  }))
})
