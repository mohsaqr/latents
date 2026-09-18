test_that("full covariance with missing indicators reproduces genuine two-level Mplus", {
  invisible(lapply(c("varying", "equal"), function(variance_model) {
    stopifnot(is.character(variance_model), length(variance_model) == 1L)
    reference <- readRDS(test_path("..", "fixtures", "mplus",
                                   paste0("twolevel-missing-full-", variance_model, ".rds")))
    fit <- fit_ml_lpa(reference$data, c("y1", "y2"), "clus", 2, 2,
                      covariance_model = "full", missing = "fiml",
                      variance_model = variance_model, n_starts = 4, seed = 20260917,
                      tol = 1e-13, max_iter = 10000)
    profile_order <- order(fit$means[, "y1"])[rank(reference$means[, 1])]
    group_order <- order(fit$profile_probabilities[, profile_order[1L]])[
      rank(reference$profile_probabilities[, 1])]
    expect_true(fit$converged)
    expect_false(fit$boundary)
    expect_equal(fit$n_parameters, reference$n_parameters)
    expect_lt(abs(fit$log_likelihood - reference$log_likelihood), 5e-5)
    expect_lt(abs(fit$aic - reference$aic), 5e-4)
    expect_lt(abs(fit$bic_individual - reference$bic_individual), 5e-4)
    expect_lt(max(abs(fit$means[profile_order, ] - reference$means)), 1e-5)
    expect_lt(max(abs(fit$covariances[, , profile_order] - reference$covariances)), 1e-5)
    expect_lt(max(abs(fit$profile_probabilities[group_order, profile_order] -
                       reference$profile_probabilities)), 1e-5)
    expect_lt(max(abs(fit$group_probabilities[group_order] - reference$group_probabilities)), 1e-5)
    expect_lt(max(abs(fit$subject_posteriors[, profile_order] - reference$subject_posteriors)), 1e-5)
    expect_lt(max(abs(fit$group_posteriors[, group_order] - reference$group_posteriors)), 1e-5)
    expect_equal(fit$indicator_data, as.matrix(reference$data[c("y1", "y2")]))
    expect_match(reference$provenance$oracle, "Genuine Mplus")
  }))
})
