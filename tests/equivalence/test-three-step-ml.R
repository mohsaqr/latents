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
  pieces <- latents:::.multilpa_level_assignments(fit, "individuals")
  varying <- latents:::.multilpa_step_ml_fit(pieces, case$data$distal, "varying")
  equal <- latents:::.multilpa_step_ml_fit(pieces, case$data$distal, "equal")
  expect_true(varying$converged)
  expect_true(equal$converged)
  note <- function(listing) {
    sprintf("Typed from Latent GOLD 6.1's listing equivalence/latentgold/returned/%s. Package side: the internal step-three ML fit (.multilpa_step_ml_fit) behind three_step(method = 'ml').",
            listing)
  }
  classes <- c("class_1", "class_2")
  # Within half a unit of Latent GOLD's last printed digit (1e-4 for the
  # likelihood of the equal-variance model, which Latent GOLD itself stops at
  # its own tolerance).
  expect_reference(varying$log_likelihood, -1858.5077, reference = "Latent GOLD",
                   covers = "three_step()", quantity = "step-three log likelihood, class-specific variances",
                   tolerance = 5e-5, measure = "absolute", note = note("c09_distal_ml_varying.lst"))
  expect_reference(equal$log_likelihood, -1858.5415, reference = "Latent GOLD",
                   covers = "three_step()", quantity = "step-three log likelihood, equal variances",
                   tolerance = 1e-4, measure = "absolute", note = note("c09_distal_ml_equal.lst"))
  # Latent GOLD numbers its classes as this package does on this case.
  expect_reference(varying$means, c(0.0029, 9.9668), reference = "Latent GOLD",
                   covers = "three_step()", quantity = "step-three means, class-specific variances",
                   tolerance = 5e-5, measure = "absolute", key = classes,
                   note = note("c09_distal_ml_varying.lst"))
  expect_reference(varying$variances, c(4.0837, 3.9502), reference = "Latent GOLD",
                   covers = "three_step()", quantity = "step-three variances, class-specific variances",
                   tolerance = 5e-5, measure = "absolute", key = classes,
                   note = note("c09_distal_ml_varying.lst"))
  expect_reference(equal$means, c(0.0004, 9.9641), reference = "Latent GOLD",
                   covers = "three_step()", quantity = "step-three means, equal variances",
                   tolerance = 5e-5, measure = "absolute", key = classes,
                   note = note("c09_distal_ml_equal.lst"))
  expect_reference(equal$variances[1L], 4.0231, reference = "Latent GOLD",
                   covers = "three_step()", quantity = "step-three common variance, equal variances",
                   tolerance = 5e-5, measure = "absolute", note = note("c09_distal_ml_equal.lst"))
})
