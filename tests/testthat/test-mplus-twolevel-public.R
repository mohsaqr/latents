test_that("matched multilevel fits reproduce Mplus on public continuous data", {
  invisible(lapply(c("varying", "equal"), function(variance_model) {
    stopifnot(is.character(variance_model), length(variance_model) == 1L)
    reference <- readRDS(test_path("..", "fixtures", "mplus",
                                   paste0("twolevel-public-", variance_model, ".rds")))
    expected <- reference$expected
    fit <- multilpa(reference$data, reference$indicators, "clus", 2L, 2L,
                      variance_model = variance_model, n_starts = 10L,
                      max_iter = 10000L, tol = 1e-14, seed = 5739)
    profile_order <- order(fit$means[, "y1"])[rank(expected$means[, 1L])]
    group_order <- order(fit$profile_probabilities[, profile_order[1L]])[
      rank(expected$profile_probabilities[, 1L])]
    expect_equal(fit$n_groups, 110L)
    expect_equal(fit$n_observations, 1000L)
    expect_true(fit$converged)
    expect_false(fit$boundary)
    expect_equal(fit$n_parameters, unname(expected$metrics["n_parameters"]))
    # Explicit ABSOLUTE tolerances: no scaling by 7,000-point log likelihoods.
    expect_lt(abs(fit$log_likelihood - expected$metrics["log_likelihood"]), 0.00005)
    expect_lt(abs(fit$aic - expected$metrics["aic"]), 0.0005)
    expect_lt(abs(fit$bic_individual - expected$metrics["bic_individual"]), 0.0005)
    # Iterative optimizers have different stopping rules; this predeclared 1e-4
    # absolute tolerance includes estimation differences, not just saved rounding.
    expect_lt(max(abs(fit$means[profile_order, ] - expected$means)), 1e-4)
    expect_lt(max(abs(fit$variances[profile_order, ] - expected$variances)), 1e-4)
    expect_lt(max(abs(fit$profile_probabilities[group_order, profile_order] -
                       expected$profile_probabilities)), 1e-4)
    expect_lt(max(abs(fit$group_probabilities[group_order] - expected$group_probabilities)), 1e-4)
    expect_lt(max(abs(fit$subject_posteriors[, profile_order] - expected$subject_posteriors)), 1e-4)
    expect_lt(max(abs(fit$group_posteriors[, group_order] - expected$group_posteriors)), 1e-4)
    expect_identical(fit$subject_profiles |> match(profile_order),
                     max.col(expected$subject_posteriors, ties.method = "first"))
    expect_identical(fit$group_classes |> match(group_order),
                     max.col(expected$group_posteriors, ties.method = "first"))
    expect_match(reference$provenance$model, "NOT original")
    expect_match(reference$provenance$version, "Mplus VERSION 9 DEMO")
  }))
})
