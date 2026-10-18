# Sampling weights in a mixture against Mplus (mplus/mixture-weights/): a
# three-profile VVI model of a Monte Carlo dataset with log-normal row
# weights, fitted by genuine Mplus 9 Demo with WEIGHT (MLR: the weighted
# pseudo log likelihood and its sandwich standard errors). Mplus scales the
# weights to sum to n, as this package does. Mplus's results file keeps the
# estimates and their standard errors to 8 significant digits, per class in
# TECH1 order: the 4 means, the 4 variances, then the class logits.

test_that("a weighted mixture reproduces Mplus's pseudo likelihood, posteriors and robust SEs", {
  directory <- file.path("mplus", "mixture-weights")
  lines <- read_mplus_lines(file.path(directory, "weights.out"))
  expect_true(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)))
  expect_true(any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)))
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(directory, "weights-results.dat"), quiet = TRUE)
  estimates <- values[seq_len(npar)]
  errors <- values[npar + seq_len(npar)]
  vars <- paste0("y", 1:4)
  data <- utils::read.table(file.path(directory, "weights.dat"),
                            col.names = c("id", vars, "w"))
  fit <- suppressWarnings(lpa(data, vars, n_profiles = 3, model = "VVI", weights = "w",
                              n_starts = 30, seed = 1, tol = 1e-10, max_iter = 5000))
  expect_identical(as.integer(fit$n_parameters), npar)
  expect_lt(abs(fit$log_likelihood - values[2L * npar + 2L]), 5e-4)
  saved <- utils::read.table(file.path(directory, "weights-posteriors.dat"))
  matched <- montecarlo_relabel(fit$subject_posteriors, as.matrix(saved[, 5:7]))
  expect_lt(matched$gap, 1e-3)
  # Mplus's estimates and SEs by class (rows: Mplus classes; columns: y1-y4).
  block <- function(from, offset) t(vapply(0:2, function(class)
    from[class * 8L + offset + 1:4], numeric(4)))
  inference <- as.data.frame(parameter_inference(fit))
  ours <- function(which, column) {
    rows <- inference[inference$parameter == which, , drop = FALSE]
    t(vapply(sprintf("profile_%d", matched$order), function(profile)
      vapply(vars, function(v) rows[[column]][rows$outcome == profile & rows$term == v],
             numeric(1)), numeric(4)))
  }
  expect_lt(max(abs(ours("mean", "estimate") - block(estimates, 0L))), 1e-3)
  expect_lt(max(abs(ours("variance", "estimate") - block(estimates, 4L))), 1e-3)
  # Robust (sandwich) standard errors. This package applies the CR1 small-sample
  # factor sqrt(n / (n - 1)) (0.9.31); Mplus's MLR does not. Measured on
  # 2026-10-18 the ratio was 1.00099 to 1.00105 for every parameter, against
  # sqrt(500 / 499) = 1.001002; with the factor removed the two agree to 1e-4.
  cr1 <- sqrt(nrow(data) / (nrow(data) - 1))
  expect_lt(max(abs(ours("mean", "standard_error") / cr1 / block(errors, 0L) - 1)), 1e-4)
  expect_lt(max(abs(ours("variance", "standard_error") / cr1 / block(errors, 4L) - 1)), 1e-4)
})
