# Latent transition models against the external references JStats pins:
# depmixS4 (Visser & Speekenbrink, 2010) panel hidden Markov fits, two Mplus
# LTA runs, and the plain-LTA row of Table 5 of Muthen & Asparouhov (2022).
# fit_transitions() with one group class is exactly this model: first-order,
# homogeneous transitions, measurement invariant over occasions.
#
# Tolerances. Mplus prints three decimals, so a value it reports can sit 5e-4
# from the unrounded one, and a criterion computed from a rounded likelihood
# another 1e-3 away. depmixS4 stops on a relative likelihood change of 1e-8,
# so its estimates are not exact maxima; 1e-4 on the likelihood and 1e-3 on
# the parameters are the gaps that stopping rule leaves at these sample sizes.

skip_if_not_installed("jsonlite")

fit_panel <- function(data, vars, k, categorical = character()) {
  fit_transitions(data, vars, "id", n_profiles = k, time = "time",
                  categorical = categorical, n_starts = 20L, seed = 1L,
                  tol = 1e-12, max_iter = 5000L)
}

test_that("Gaussian LTA reproduces depmixS4 on both JStats scenarios", {
  reference <- jstats_fixture("gaussian-lta-ref.json")
  invisible(lapply(c("scenario1", "scenario2"), function(scenario) {
    expected <- jstats_gaussian(reference[[scenario]])
    data <- jstats_panel(expected$data, "y")
    vars <- paste0("y", seq_len(expected$d))
    fit <- fit_panel(data, vars, expected$k)
    p <- match_profiles(fit$means[, 1L], expected$mu[, 1L])

    expect_true(fit$converged, label = scenario)
    # depmixS4's `df` is the free-parameter count; its `npar_total` also
    # counts the probabilities fixed by the sum-to-one constraints.
    expect_identical(fit$n_parameters, as.integer(expected$df), label = scenario)
    expect_lt(abs(fit$log_likelihood - expected$logLikelihood), 1e-4)
    expect_lt(abs(fit$aic - expected$aic), 2e-4)
    # depmixS4's BIC counts N x T observations: multilpa's individual level.
    expect_lt(abs(fit$bic_individual - expected$bic_depmixS4), 2e-4)
    expect_lt(max(abs(fit$means[p, ] - expected$mu)), 1e-3)
    expect_lt(max(abs(sqrt(fit$variances[p, , drop = FALSE]) - expected$sd)), 1e-3)
    expect_lt(max(abs(fit$initial_probabilities[1L, p] - expected$pi)), 1e-3)
    expect_lt(max(abs(fit$transition_probabilities[p, p, 1L] - expected$tau)), 1e-3)
  }))
})

test_that("binary LTA reproduces depmixS4 where depmixS4 reached the maximum", {
  expected <- jstats_fixture("lta-ref.json")$scenario1
  data <- jstats_panel(expected$data, "u")
  vars <- paste0("u", seq_len(expected$m))
  fit <- fit_panel(data, vars, expected$k, categorical = vars)
  success <- success_probabilities(fit)
  p <- match_profiles(success[, 1L], expected$rho[, 1L])

  expect_true(fit$converged)
  expect_identical(fit$n_parameters, as.integer(expected$df))
  expect_lt(abs(fit$log_likelihood - expected$logLikelihood), 1e-4)
  expect_lt(abs(fit$aic - expected$aic), 2e-4)
  expect_lt(abs(fit$bic_individual - expected$bic), 2e-4)
  expect_lt(max(abs(success[p, ] - expected$rho)), 1e-3)
  expect_lt(max(abs(fit$initial_probabilities[1L, p] - expected$pi)), 1e-3)
  expect_lt(max(abs(fit$transition_probabilities[p, p, 1L] - expected$tau)), 1e-3)
})

test_that("on the three-profile binary scenario multilpa finds a higher maximum than depmixS4", {
  # depmixS4's recorded solution for this scenario is a local maximum. The
  # claim is checked with a likelihood that shares no code with the package:
  # brute-force summation over all 3^4 latent paths, evaluated at both
  # solutions. The same function must first reproduce depmixS4's own number
  # at depmixS4's own estimates, or it proves nothing.
  expected <- jstats_fixture("lta-ref.json")$scenario2
  data <- jstats_panel(expected$data, "u")
  vars <- paste0("u", seq_len(expected$m))
  fit <- fit_panel(data, vars, expected$k, categorical = vars)

  at_reference <- path_log_likelihood(data, vars, expected$pi, expected$tau,
                                      binary_emission(expected$rho))
  at_multilpa <- path_log_likelihood(data, vars, fit$initial_probabilities[1L, ],
                                     fit$transition_probabilities[, , 1L],
                                     binary_emission(success_probabilities(fit)))

  expect_identical(fit$n_parameters, as.integer(expected$df))
  expect_lt(abs(at_reference - expected$logLikelihood), 1e-6)
  expect_lt(abs(at_multilpa - fit$log_likelihood), 1e-6)
  expect_gt(fit$log_likelihood - expected$logLikelihood, 1e-3)
})

test_that("binary LTA reproduces the two Mplus LTA runs pinned by JStats", {
  reference <- jstats_fixture("mplus-ref.json")
  invisible(lapply(c("lta_k2_t2_m3", "lta_k2_t3_m2"), function(scenario) {
    run <- reference[[scenario]]
    expected <- run$expected
    data <- jstats_panel(run$data, "u")
    vars <- paste0("u", seq_len(run$m))
    fit <- fit_panel(data, vars, run$k, categorical = vars)
    success <- success_probabilities(fit)
    p <- match_profiles(success[, 1L], expected$rho[, 1L])
    criteria <- get_data(fit, "information_criteria")

    expect_true(fit$converged, label = scenario)
    expect_identical(fit$n_parameters, as.integer(expected$npar), label = scenario)
    expect_lt(abs(fit$log_likelihood - expected$logLikelihood), 5e-4)
    # Mplus counts subjects for BIC and sample-size-adjusted BIC.
    expect_lt(abs(fit$aic - expected$aic), 1.5e-3)
    expect_lt(abs(fit$bic_groups - expected$bic), 1.5e-3)
    expect_lt(abs(criteria$sabic_groups - expected$abic), 1.5e-3)
    expect_lt(max(abs(fit$initial_probabilities[1L, p] - expected$pi)), 6e-4)
    expect_lt(max(abs(fit$transition_probabilities[p, p, 1L] - expected$tau)), 6e-4)
    expect_lt(max(abs(success[p, ] - expected$rho)), 1e-3)
  }))
})

test_that("plain LTA reproduces Table 5 of Muthen and Asparouhov (2022)", {
  # Eid/Langeheine mood panel, N = 494 as a frequency-weighted pattern table.
  # Table 5's regular-LTA row: stationary transitions, measurement invariant.
  anchor <- jstats_fixture("ri-lta-mplus-anchor.json")
  expected <- anchor$regular_lta
  data <- jstats_patterns(anchor$patterns, anchor$dims$T, anchor$dims$M, "u")
  vars <- paste0("u", seq_len(anchor$dims$M))
  fit <- fit_panel(data, vars, anchor$dims$K, categorical = vars)
  # Table 5 reports only the transition matrix among the estimates; its
  # diagonal identifies the labelling.
  diagonal <- diag(fit$transition_probabilities[, , 1L])
  p <- match_profiles(diagonal, diag(expected$transition_stationary))
  criteria <- get_data(fit, "information_criteria")

  expect_identical(fit$n_groups, as.integer(anchor$dims$N))
  expect_true(fit$converged)
  expect_identical(fit$n_parameters, as.integer(expected$nFreeParameters))
  expect_lt(abs(fit$log_likelihood - expected$logLikelihood), 5e-4)
  expect_lt(abs(fit$aic - expected$aic), 1.5e-3)
  expect_lt(abs(fit$bic_groups - expected$bic), 1.5e-3)
  expect_lt(abs(criteria$sabic_groups - expected$sabic), 1.5e-3)
  expect_lt(max(abs(fit$transition_probabilities[p, p, 1L] -
                      expected$transition_stationary)), 6e-4)
})

test_that("the brute-force path likelihood reproduces a Gaussian reference", {
  # Calibrates path_log_likelihood() on a continuous model, so the local
  # maximum finding above does not rest on a function checked only once.
  expected <- jstats_gaussian(jstats_fixture("gaussian-lta-ref.json")$scenario2)
  data <- jstats_panel(expected$data, "y")
  vars <- paste0("y", seq_len(expected$d))
  value <- path_log_likelihood(data, vars, expected$pi, expected$tau,
                               gaussian_emission(expected$mu, expected$sd))
  expect_lt(abs(value - expected$logLikelihood), 1e-6)
})
