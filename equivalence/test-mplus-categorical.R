.categorical_indicators <- function() c("u1", "u2", "u3", "u4", "u5")

.categorical_fit <- function(data, ...) {
  multilpa(data, .categorical_indicators(), "clus", n_profiles = 2,
           n_group_classes = 2, categorical = .categorical_indicators(),
           n_starts = 8, seed = 4, ...)
}

test_that("thresholds and their standard errors reproduce a genuine Mplus run", {
  fixture <- readRDS(equivalence_fixture("mplus", "twolevel-categorical.rds"))
  fit <- .categorical_fit(fixture$data)
  inference <- parameter_inference(fit, fixture$data)
  theta <- multilpa:::.multilpa_coefficients(fit, "unconstrained")
  # The tidy decomposition is the canonical one, so the block is selected by
  # its `parameter` column rather than by parsing a name back apart.
  labels <- multilpa:::.multilpa_coefficient_labels(fit, "unconstrained")
  logits <- which(labels$parameter == "response_logit")
  errors <- sqrt(diag(attr(inference, "covariance_unconstrained")))[logits]

  # For a binary indicator the unconstrained coordinate log(p0 / p1) is exactly
  # the threshold Mplus reports, so the two are directly comparable.
  order_by_first <- order(fit$response_probabilities[[1L]][, 1L], decreasing = TRUE)
  mplus_threshold <- c(1.451, 1.052, 1.392, -1.471, -1.727,
                       -1.944, -1.508, -1.153, 1.129, 1.792)
  mplus_error <- c(0.118, 0.101, 0.112, 0.116, 0.131,
                   0.145, 0.117, 0.105, 0.104, 0.139)
  profile <- as.integer(sub("^profile_", "", labels$outcome[logits]))
  indicator <- sub(":.*$", "", labels$term[logits])
  key <- paste(match(profile, order_by_first), indicator)
  expected_key <- paste(rep(c(1L, 2L), each = 5), .categorical_indicators())

  expect_equal(unname(theta[logits])[match(expected_key, key)], mplus_threshold,
               tolerance = 1e-3)
  expect_equal(unname(errors)[match(expected_key, key)], mplus_error,
               tolerance = 2e-3)
  expect_equal(fit$log_likelihood, fixture$mplus_log_likelihood, tolerance = 1e-6)
})

test_that("categorical fits reproduce genuine Mplus two-level results", {
  reference <- readRDS(equivalence_fixture("mplus", "twolevel-categorical.rds"))
  vars <- paste0("u", 1:5)
  fit <- multilpa(reference$data, vars, "clus", n_profiles = 2,
                    n_group_classes = 2, categorical = vars,
                    n_starts = 40, seed = 20260918, tol = 1e-13, max_iter = 20000)
  expect_true(fit$converged)
  expect_equal(fit$n_parameters, as.numeric(reference$n_parameters))
  thresholds <- vapply(fit$response_probabilities, function(block) {
    as.vector(.multilpa_categorical_thresholds(block))
  }, numeric(2))
  profile_order <- order(thresholds[, 1L])[rank(reference$mplus_thresholds[, 1L])]
  group_order <- order(fit$profile_probabilities[, profile_order[1L]])[
    rank(reference$mplus_profile_probabilities[, 1L])]
  expect_lt(max(abs(thresholds[profile_order, ] - reference$mplus_thresholds)), 1e-4)
  expect_lt(max(abs(fit$profile_probabilities[group_order, profile_order] -
                      reference$mplus_profile_probabilities)), 1e-5)
  expect_lt(max(abs(fit$group_probabilities[group_order] -
                      reference$mplus_group_probabilities)), 1e-5)
  expect_lt(abs(fit$log_likelihood - reference$mplus_log_likelihood), 1e-3)
})

