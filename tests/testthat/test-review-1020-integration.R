test_that("paired mixed measurement retains direct-effect inference through data adapters", {
  set.seed(12020)
  n <- 400
  data <- data.frame(g = rep(seq_len(40), each = 10), z = stats::rnorm(n),
                     a = factor(stats::rbinom(n, 1, 0.4)),
                     b = factor(stats::rbinom(n, 1, 0.6)),
                     ordinal = ordered(sample(0:2, n, TRUE)),
                     count = stats::rpois(n, 3))
  data$y <- 1 + 0.7 * data$z + stats::rnorm(n)
  data$weight <- rep(rep(c(1, 2), 20), each = 10)
  vars <- c("a", "b", "ordinal", "count", "y")
  fit <- function(frame) quietly(multilpa(frame, vars, "g", 1, 1,
    categorical = c("a", "b"), ordinal = "ordinal", count = "count",
    local_dependence = "a:b", direct_effects = "y ~ z", weights = "weight",
    missing = "fiml", n_starts = 1, seed = 1, tol = 1e-12))
  original <- fit(data)
  shifted_data <- transform(data, z = z + 2026)
  shifted <- fit(shifted_data)
  original_inference <- parameter_inference(original)
  shifted_inference <- parameter_inference(shifted)
  direct_rows <- original_inference$parameter == "direct_effect"
  expect_equal(shifted_inference$estimate[direct_rows],
               original_inference$estimate[direct_rows], tolerance = 1e-8)
  expect_equal(shifted_inference$standard_error[direct_rows],
               original_inference$standard_error[direct_rows], tolerance = 1e-7)
  expect_equal(predict(shifted, shifted_data, type = "density"),
               predict(original, data, type = "density"), tolerance = 1e-8)
  # The Gaussian block is a one-profile weighted normal regression even
  # when independent categorical, ordinal and count blocks are present.
  reference <- stats::lm(y ~ z, data = shifted_data, weights = weight)
  means <- shifted_inference$parameter == "mean"
  expect_equal(shifted_inference$estimate[c(which(means), which(direct_rows))],
               unname(stats::coef(reference)), tolerance = 1e-8)
  frame <- .multilpa_resolve_data(shifted, NULL)
  expect_equal(predict(shifted, frame, type = "density"),
               predict(shifted, shifted_data, type = "density"), tolerance = 1e-10)
  refitted <- do.call(multilpa, c(list(data = frame, n_starts = 1, seed = 1,
    tol = 1e-12), .multilpa_refit_arguments(shifted)))
  expect_equal(refitted$log_likelihood, shifted$log_likelihood, tolerance = 1e-9)
  expect_equal(refitted$n_parameters, shifted$n_parameters)
  expect_equal(attr(parameter_inference(refitted), "covariance"),
               attr(shifted_inference, "covariance"), tolerance = 1e-7)
  # Independent empty groups with extreme covariates supply neither
  # information nor regression geometry, including in mixed measurement.
  extra <- shifted_data[seq_len(20), ]
  extra$g <- rep(41:42, each = 10)
  extra[vars] <- NA
  extra$z <- 1e5
  extra$weight <- 0
  expanded <- fit(rbind(shifted_data, extra))
  expanded_inference <- parameter_inference(expanded)
  expect_true(expanded$converged)
  expect_equal(expanded_inference$estimate, shifted_inference$estimate, tolerance = 1e-7)
  expect_equal(attr(expanded_inference, "covariance"),
               attr(shifted_inference, "covariance"), tolerance = 1e-6)
})
