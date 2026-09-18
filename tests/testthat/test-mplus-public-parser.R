test_that("the public Mplus parser handles known synthetic output and label reversal", {
  synthetic <- c("Number of observations 100", "Number of Free Parameters 17",
                 " H0 Value -123.456", " Akaike (AIC) 280.912", " Bayesian (BIC) 325.199",
                 " Entropy 0.900", "Latent Class 1", "", " Means",
                 paste0("    Y", 1:4, " 2.000 0.100 20.000 0.000"), "", " Variances",
                 paste0("    Y", 1:4, " 1.000 0.100 10.000 0.000"), "", "Latent Class 2",
                 "", " Means", paste0("    Y", 1:4, " -2.000 0.100 -20.000 0.000"),
                 "", " Variances", paste0("    Y", 1:4, " 3.000 0.100 30.000 0.000"), "",
                 "BASED ON THE ESTIMATED MODEL", "", " Latent", " Classes", "",
                 " 1 40.00000 0.40000", " 2 60.00000 0.60000", rep("", 5L),
                 "BASED ON THEIR MOST LIKELY LATENT CLASS MEMBERSHIP", "", "Class Counts",
                 "", " Latent", " Classes", "", " 1 39 0.39000", " 2 61 0.61000", "")
  parsed <- parse_public_mplus(synthetic)
  expect_equal(parsed$means, matrix(rep(c(-2, 2), 4L), nrow = 2L))
  expect_equal(parsed$variances, matrix(rep(c(3, 1), 4L), nrow = 2L))
  expect_equal(parsed$proportions, c(0.6, 0.4))
  expect_equal(parsed$effective_counts, c(60, 40))
  expect_equal(parsed$modal_counts, c(61, 39))
  expect_equal(parsed$log_likelihood, -123.456)
  expect_equal(parsed$aic, 280.912)
  expect_equal(parsed$bic_individual, 325.199)
  expect_equal(parsed$n_parameters, 17)
  expect_equal(parsed$n_observations, 100)
  expect_equal(parsed$entropy, 0.9)
  expect_error(parse_public_mplus("missing output"))
})

