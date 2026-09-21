test_that("robust standard errors reproduce genuine Mplus MLR results", {
  reference <- readRDS(equivalence_fixture("mplus", "twolevel-robust.rds"))
  dat <- reference$data
  fit <- multilpa(dat, c("y1", "y2"), "clus", 2, 2, variance_model = "varying",
                    start = starting_values(reference), n_starts = 1,
                    max_iter = 5000, tol = 1e-13)
  information <- parameter_inference(fit, dat, vcov_type = "robust")
  expect_identical(attr(information, "vcov_type"), "robust")
  q <- fit$n_parameters
  n_measurement <- q - 3L
  # Mplus reports a group logit, a within-class intercept at cb = 2, and the
  # cb = 1 minus cb = 2 contrast, with variances rather than log variances.
  ordering <- c(1, 2, 5, 6, 3, 4, 7, 8)
  derivatives <- c(rep(1, 4), as.vector(t(fit$variances)))
  jacobian <- matrix(0, q, q)
  jacobian[cbind(seq_len(n_measurement), ordering)] <- derivatives[ordering]
  jacobian[q - 2L, q] <- 1
  jacobian[q - 1L, q - 1L] <- 1
  jacobian[q, c(q - 2L, q - 1L)] <- c(1, -1)
  standard_errors <- sqrt(diag(jacobian %*% attr(information, "covariance_unconstrained") %*%
                                 t(jacobian)))
  expect_lt(max(abs(standard_errors - reference$mplus_standard_errors)), 1e-6)
  expect_lt(abs(attr(information, "scaling_correction") - reference$mplus_scaling), 1e-6)
  indices <- get_data(fit, "information_criteria", format = "long")
  # A criterion that has no sample-size convention -- AIC counts parameters, not
  # units -- carries NA there rather than a sentinel level, so it is selected
  # by NA and not by a magic string.
  select <- function(criterion, convention) {
    wanted <- if (is.na(convention)) is.na(indices$convention) else
      !is.na(indices$convention) & indices$convention == convention
    indices$value[indices$criterion == criterion & wanted]
  }
  expect_lt(abs(select("aic", NA_character_) - reference$mplus_aic), 1e-3)
  expect_lt(abs(select("bic", "individuals") - reference$mplus_bic), 1e-3)
  expect_lt(abs(select("sabic", "individuals") - reference$mplus_sabic), 1e-3)
})

