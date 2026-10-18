# `scaling = "grand"` against Mplus's DEFINE: STANDARDIZE (genuine Mplus 9
# Demo, mplus/standardize/). Mplus divides by the standard deviation with
# denominator n; this package, like sd() and scale(), uses n - 1. The two
# therefore differ by a known constant and nothing else: standardized with
# denominator n by hand, the package reproduces Mplus exactly, and with its
# own convention its log likelihood is Mplus's plus
# n * sum(log(sqrt(n / (n - 1)))) over the indicators, with the same classes.

test_that("grand scaling matches Mplus STANDARDIZE up to the n versus n - 1 constant", {
  directory <- file.path("mplus", "standardize")
  lines <- read_mplus_lines(file.path(directory, "standardize.out"))
  expect_true(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)))
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(directory, "standardize-results.dat"), quiet = TRUE)
  mplus_ll <- values[2L * npar + 2L]
  vars <- paste0("y", 1:4)
  data <- utils::read.table(file.path(directory, "standardize.dat"),
                            col.names = c("id", vars))
  n <- nrow(data)
  fit <- function(frame, ...) lpa(frame, vars, n_profiles = 3, model = "VVI",
                                  n_starts = 20, seed = 1, tol = 1e-10, ...)
  # Denominator n by hand: Mplus exactly.
  by_n <- data
  by_n[vars] <- lapply(data[vars], \(v) (v - mean(v)) / sqrt(mean((v - mean(v))^2)))
  hand <- fit(by_n)
  expect_identical(as.integer(hand$n_parameters), npar)
  expect_lt(abs(hand$log_likelihood - mplus_ll), 5e-4)
  # The package's own convention: the same fit, shifted by the known constant.
  scaled <- fit(data, centering = "grand", scaling = "grand")
  constant <- n * length(vars) * log(sqrt(n / (n - 1)))
  expect_lt(abs(scaled$log_likelihood - (mplus_ll + constant)), 5e-4)
  expect_identical(max.col(scaled$subject_posteriors), max.col(hand$subject_posteriors))
  # And the posteriors are Mplus's, row by row.
  saved <- utils::read.table(file.path(directory, "standardize-posteriors.dat"))
  posteriors <- as.matrix(saved[, 5:7])
  matched <- montecarlo_relabel(scaled$subject_posteriors, posteriors)
  expect_lt(matched$gap, 1e-3)
})
