# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Classification summaries of latent profile fits against tidyLPA (Rosenberg et
# al., 2018), whose get_fit() reports, from its own mclust fit of the same
# model, the relative entropy (Entropy), the smallest and largest modal class
# proportion (n_min, n_max) and the sample size (n). Its prob_min and prob_max
# are not compared: they are the diagonal of the column-normalised matrix
# (Mplus's classification probabilities for the most likely class by latent
# class), not the average posteriors of get_avepp() (on the faithful data
# 0.99979 against the AvePP 0.99961).
# Two data sets: the Old Faithful eruptions (two profiles, mclust's VVI) and
# 600 rows drawn by simulab::simulate_lpa() from three profiles with unit
# variances (mclust's EEI).
#
# tidyLPA is given mclust's EM control at a relative change of 1e-12 (its
# default, 1e-5, stops 9e-4 below the maximum of the simulab data); latents
# stops at 1e-10. Tolerances, fixed before the first run: 1e-7 relative on the
# likelihood, 1e-4 relative on the entropy (which moves with estimates that
# agree to about 1e-5), exact (1e-12) on the modal
# proportions, ratios of counts. tidyLPA also draws a bootstrap likelihood
# ratio test that is not compared; its draws are seeded.

faithful_fit <- lpa(datasets::faithful, c("eruptions", "waiting"), n_profiles = 2, n_starts = 5, seed = 1)
simulated_lpa <- as.data.frame(simulab::simulate_lpa(600, means = rbind(c(0, 0, 0), c(3, 3, 0), c(0, 3, 3)),
                                                     sds = 1, proportions = c(0.5, 0.3, 0.2), seed = 11))
simulated_vars <- c("indicator_1", "indicator_2", "indicator_3")
simulated_fit <- lpa(simulated_lpa, simulated_vars, n_profiles = 3, model = "EEI", n_starts = 10, seed = 1)

tidylpa_fit <- function(data, n_profiles, variances) {
  printed <- utils::capture.output(profiles <- withr::with_seed(1L, tidyLPA::estimate_profiles(
    data, n_profiles = n_profiles, variances = variances, covariances = "zero",
    control = mclust::emControl(tol = c(1e-12, 1e-12)))))
  cbind(level = "individuals", as.data.frame(tidyLPA::get_fit(profiles)))
}

# The summaries tidyLPA reports, read from latents' accessors.
classification_summaries <- function(fit) {
  classification <- subset(get_classification(fit), level == "individuals")
  data.frame(level = "individuals", Entropy = subset(get_entropy(fit), level == "individuals")$relative_entropy,
             n_min = min(classification$proportion_modal), n_max = max(classification$proportion_modal))
}

check_against_tidylpa <- function(fit, reference, label) {
  expect_reference(as.numeric(logLik(fit)), reference$LogLik, reference = "tidyLPA",
                   covers = c("lpa()", "logLik.multilpa()"),
                   quantity = sprintf("log likelihood, %s", label), tolerance = 1e-7, measure = "relative")
  summaries <- classification_summaries(fit)
  expect_reference(summaries[c("level", "Entropy")], reference[c("level", "Entropy")],
                   reference = "tidyLPA", covers = "get_entropy()",
                   quantity = sprintf("relative entropy, %s", label),
                   tolerance = 1e-4, measure = "relative", key = "level")
  expect_reference(summaries[c("level", "n_min", "n_max")], reference[c("level", "n_min", "n_max")], reference = "tidyLPA",
                   covers = "get_classification()",
                   quantity = sprintf("smallest and largest modal class proportion, %s", label),
                   tolerance = 1e-12, measure = "relative", key = "level")
  expect_reference(nobs(fit), reference$n, reference = "tidyLPA", covers = "nobs.multilpa()",
                   quantity = sprintf("number of observations, %s", label), tolerance = 0)
}

test_that("entropy and modal proportions of the faithful fit are tidyLPA's", {
  skip_without_reference("tidyLPA")
  check_against_tidylpa(faithful_fit, tidylpa_fit(datasets::faithful, 2, "varying"), "faithful, VVI")
})

test_that("entropy and modal proportions of a simulab three-profile fit are tidyLPA's", {
  skip_without_reference("tidyLPA")
  check_against_tidylpa(simulated_fit, tidylpa_fit(simulated_lpa[simulated_vars], 3, "equal"),
                        "simulab three profiles, EEI")
})

test_that("diagnostics() reports tidyLPA's entropy and modal proportions", {
  skip_without_reference("tidyLPA")
  reference <- tidylpa_fit(datasets::faithful, 2, "varying")
  diagnosed <- diagnostics(faithful_fit, plots = FALSE)
  entropy <- data.frame(level = "individuals", relative_entropy = reference$Entropy)
  expect_reference(subset(as.data.frame(diagnosed), level == "individuals", select = c(level, relative_entropy)),
                   entropy, reference = "tidyLPA", covers = c("diagnostics()", "as.data.frame.multilpa_diagnostics()"),
                   quantity = "relative entropy, diagnostics()", tolerance = 1e-4, measure = "relative", key = "level")
  expect_reference(subset(as.data.frame(summary(diagnosed)), level == "individuals", select = c(level, relative_entropy)),
                   entropy, reference = "tidyLPA",
                   covers = c("summary.multilpa_diagnostics()", "as.data.frame.multilpa_diagnostics()"),
                   quantity = "relative entropy, summary(diagnostics())", tolerance = 1e-4, measure = "relative",
                   key = "level")
  classification <- get_results(diagnosed, "classification")
  expect_reference(data.frame(level = "individuals", n_min = min(classification$proportion_modal),
                              n_max = max(classification$proportion_modal)),
                   reference[c("level", "n_min", "n_max")], reference = "tidyLPA",
                   covers = c("diagnostics()", "get_results.multilpa_diagnostics()"),
                   quantity = "extreme modal proportions, diagnostics()", tolerance = 1e-12, measure = "relative",
                   key = "level")
  checks <- subset(get_checks(faithful_fit), check == "relative entropy", select = c(check, value))
  expect_reference(checks, data.frame(check = "relative entropy", value = reference$Entropy),
                   reference = "tidyLPA", covers = "get_checks()",
                   quantity = "relative entropy, get_checks()", tolerance = 1e-4,
                   measure = "relative", key = "check")
})
