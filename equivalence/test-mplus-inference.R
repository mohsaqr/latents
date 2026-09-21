test_that("ordinary ML Hessian standard errors reproduce genuine Mplus9 results", {
  # These are the q standard errors following q estimates in the retained
  # validation/mplus/generated/{varying,equal}-results.dat files. TECH1 fixes
  # their ordering; the accompanying RDS fixtures retain file MD5 provenance.
  mplus_standard_errors <- list(
    varying = c(0.025488887, 0.033533281, 0.022282673, 0.038562637,
      0.024286238, 0.030894151, 0.020636523, 0.033392006,
      0.25820121, 0.11595498, 0.16936181),
    equal = c(0.02512813, 0.032521195, 0.015157493, 0.025388461,
      0.024630471, 0.03187695, 0.25820121, 0.1234565, 0.16936146))
  invisible(lapply(names(mplus_standard_errors), function(variance_model) {
    stopifnot(is.character(variance_model))
    reference <- readRDS(equivalence_fixture("mplus",
      paste0("twolevel-synthetic-", variance_model, ".rds")))
    parameters <- reference[c("means", "variances", "profile_probabilities", "group_probabilities")]
    fit <- multilpa(reference$data, c("y1", "y2"), "clus", 2, 2,
      variance_model = variance_model, start = parameters, n_starts = 1,
      max_iter = 5000, tol = 1e-13)
    information <- parameter_inference(fit, reference$data)
    q <- fit$n_parameters
    n_measurement <- q - 3L
    # Natural measurement estimates first; Mplus uses group logit, a CW
    # intercept at CB=2, and a CW logit contrast CB=1 minus CB=2.
    jacobian <- matrix(0, q, q)
    ordering <- if (variance_model == "varying") c(1, 2, 5, 6, 3, 4, 7, 8) else c(1, 2, 5, 6, 3, 4)
    measurement_derivatives <- c(rep(1, 4), if (variance_model == "varying")
      as.vector(t(fit$variances)) else fit$variances[1L, ])
    jacobian[cbind(seq_len(n_measurement), ordering)] <- measurement_derivatives[ordering]
    jacobian[q - 2L, q] <- 1
    jacobian[q - 1L, q - 1L] <- 1
    jacobian[q, c(q - 2L, q - 1L)] <- c(1, -1)
    mplus_scale_se <- sqrt(diag(jacobian %*% attr(information, "covariance_unconstrained") %*% t(jacobian)))
    # Mplus and finite-difference score derivatives differ by < 4e-8 here.
    expect_lt(max(abs(mplus_scale_se - mplus_standard_errors[[variance_model]])), 1e-6)
    expect_lt(attr(information, "scaled_score"), 1e-5)
    expect_equal(information$standard_error,
      parameter_inference(fit, reference$data, step = 5e-5)$standard_error, tolerance = 1e-6)
  }))
})

test_that("full-covariance FIML Hessian standard errors reproduce genuine Mplus9", {
  invisible(lapply(c("varying", "equal"), function(variance_model) {
    stopifnot(is.character(variance_model))
    reference <- readRDS(equivalence_fixture("mplus",
      paste0("twolevel-missing-full-", variance_model, ".rds")))
    start <- reference[c("means", "covariances", "profile_probabilities", "group_probabilities")]
    fit <- multilpa(reference$data, c("y1", "y2"), "clus", 2, 2,
      variance_model = variance_model, covariance_model = "full", missing = "fiml",
      start = start, n_starts = 1, tol = 1e-13, max_iter = 5000)
    information <- parameter_inference(fit, reference$data)
    q <- fit$n_parameters
    n_measurement <- q - 3L
    ordering <- if (variance_model == "varying") c(1, 2, 5, 6, 7, 3, 4, 8, 9, 10) else
      c(1, 2, 5, 6, 7, 3, 4)
    jacobian <- matrix(0, q, q)
    jacobian[seq_len(n_measurement), ] <- .multilpa_inference_jacobian(fit)[ordering, ]
    jacobian[q - 2L, q] <- 1
    jacobian[q - 1L, q - 1L] <- 1
    jacobian[q, c(q - 2L, q - 1L)] <- c(1, -1)
    observed_standard_errors <- sqrt(diag(jacobian %*% attr(information, "covariance_unconstrained") %*% t(jacobian)))
    expect_lt(max(abs(observed_standard_errors - reference$mplus_standard_errors)), 1e-6)
    expect_lt(attr(information, "scaled_score"), 1e-4)
  }))
})

