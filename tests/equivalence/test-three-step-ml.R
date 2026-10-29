# Comparisons with independent references moved from tests/testthat/test-three-step-ml.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("both variance models converge and agree with Latent GOLD's step three", {
  # Latent GOLD 6.1, `step3 modal ml` on this package's case c09 with the
  # class sizes fixed at the average posteriors (equivalence/latentgold/):
  # log likelihoods, means and variances as printed there.
  skip_if_not(file.exists(test_path("..", "..", "equivalence", "latentgold", "cases.R")),
              "the Latent GOLD kit is not available")
  source(test_path("..", "..", "equivalence", "latentgold", "cases.R"), local = TRUE)
  case <- .lg_step_case(.lg_step_data(seed = 109L))
  fit <- case$fit()
  pieces <- .multilpa_level_assignments(fit, "individuals")
  varying <- .multilpa_step_ml_fit(pieces, case$data$distal, "varying")
  equal <- .multilpa_step_ml_fit(pieces, case$data$distal, "equal")
  expect_true(varying$converged)
  expect_true(equal$converged)
  # Within half a unit of Latent GOLD's last printed digit (1e-4 for the
  # likelihood of the equal-variance model, which Latent GOLD itself stops at
  # its own tolerance).
  within <- function(actual, printed, half_unit = 5e-5) {
    expect_true(all(abs(actual - printed) <= half_unit),
                info = paste(format(actual, digits = 8), collapse = " "))
  }
  within(varying$log_likelihood, -1858.5077)
  within(equal$log_likelihood, -1858.5415, 1e-4)
  # Latent GOLD numbers its classes as this package does on this case.
  within(varying$means, c(0.0029, 9.9668))
  within(varying$variances, c(4.0837, 3.9502))
  within(equal$means, c(0.0004, 9.9641))
  within(equal$variances[1L], 4.0231)
})
