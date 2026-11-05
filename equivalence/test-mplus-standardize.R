# `scaling = "grand"` against Mplus's DEFINE: STANDARDIZE (genuine Mplus 9
# Demo, mplus/standardize/). Mplus divides by the standard deviation with
# denominator n; grand scaling uses the same population standard deviation.
# Both hand standardization and the package's transform therefore reproduce
# Mplus directly, including its likelihood, without a sample-SD correction.

test_that("grand scaling matches Mplus STANDARDIZE with the population SD", {
  directory <- file.path("mplus", "standardize")
  lines <- read_mplus_lines(file.path(directory, "standardize.out"))
  expect_true(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)))
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(directory, "standardize-results.dat"), quiet = TRUE)
  mplus_ll <- values[2L * npar + 2L]
  vars <- paste0("y", 1:4)
  data <- utils::read.table(file.path(directory, "standardize.dat"),
                            col.names = c("id", vars))
  fit <- function(frame, ...) lpa(frame, vars, n_profiles = 3, model = "VVI",
                                  n_starts = 20, seed = 1, tol = 1e-10, ...)
  # Denominator n by hand: Mplus exactly.
  by_n <- data
  by_n[vars] <- lapply(data[vars], \(v) (v - mean(v)) / sqrt(mean((v - mean(v))^2)))
  hand <- fit(by_n)
  expect_identical(as.integer(hand$n_parameters), npar)
  expect_lt(abs(hand$log_likelihood - mplus_ll), 5e-4)
  # The built-in transform uses the same denominator, with no likelihood shift.
  scaled <- fit(data, centering = "grand", scaling = "grand")
  expect_lt(abs(scaled$log_likelihood - mplus_ll), 5e-4)
  expect_equal(scaled$scaling_factors[1L, ],
                vapply(data[vars], function(v) sqrt(mean((v - mean(v))^2)), numeric(1)),
                tolerance = 1e-12)
  expect_identical(max.col(scaled$subject_posteriors), max.col(hand$subject_posteriors))
  # And the posteriors are Mplus's, row by row.
  saved <- utils::read.table(file.path(directory, "standardize-posteriors.dat"))
  posteriors <- as.matrix(saved[, 5:7])
  matched <- montecarlo_relabel(scaled$subject_posteriors, posteriors)
  expect_lt(matched$gap, 1e-3)
})
