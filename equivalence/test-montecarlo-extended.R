# Extended Monte Carlo equivalence (montecarlo/extended/): two-level profiles
# with a full within-profile covariance, invariant (EEE) and profile-specific
# (VVV), against Mplus and Latent GOLD; and latent transition models with
# missing indicators and with variances shared by the states, against Latent
# GOLD's Markov module. 4 configurations x 10 simulated datasets.

extended_dir <- file.path(montecarlo_dir, "extended")
extended_cases <- function() readRDS(file.path(extended_dir, "kit", "cases.rds"))

extended_fits <- montecarlo_fit_all(extended_cases(), n_starts = 30L)

test_that("the extended kit holds 4 configurations of 10, each with its outputs", {
  cases <- extended_cases()
  expect_identical(length(cases), 40L)
  expect_false(any(vapply(extended_fits, inherits, logical(1), "try-error")))
  returned <- file.path(extended_dir, "returned")
  invisible(lapply(names(cases), function(name) {
    expect_true(file.exists(file.path(returned, paste0(name, ".lst"))), label = name)
    if ("mplus" %in% cases[[name]]$spec$programs) {
      expect_true(file.exists(file.path(returned, paste0(name, ".out"))), label = name)
    }
  }))
})

test_that("two-level full-covariance models reproduce Mplus", {
  cases <- extended_cases()
  two <- names(cases)[vapply(cases, \(case) "mplus" %in% case$spec$programs, logical(1))]
  expect_identical(length(two), 20L)
  invisible(lapply(two, function(name) {
    case <- cases[[name]]
    fit <- extended_fits[[name]]
    m <- montecarlo_mplus(name, case, root = extended_dir)
    expect_true(m$terminated && m$replicated, label = name)
    expect_identical(as.integer(fit$n_parameters), m$n_parameters, label = name)
    expect_lt(abs(fit$log_likelihood - m$log_likelihood), 5e-4, label = paste(name, "LL"))
    expect_identical(fit$n_informative, m$n_observations, label = name)
    expect_false(anyDuplicated(m$id) > 0L)
    rows <- match(m$id, case$data$id)
    expect_lt(montecarlo_relabel(fit$subject_posteriors[rows, , drop = FALSE],
                                 m$posteriors$profiles)$gap, 1e-3,
              label = paste(name, "posteriors"))
    matched <- montecarlo_relabel(fit$group_posteriors[fit$group_index[rows], , drop = FALSE],
                                  m$posteriors$groups)
    expect_lt(matched$gap, 1e-3, label = paste(name, "group posteriors"))
  }))
})

test_that("every extended case reproduces Latent GOLD, every row included", {
  cases <- extended_cases()
  invisible(lapply(names(cases), function(name) {
    case <- cases[[name]]
    fit <- extended_fits[[name]]
    lg <- montecarlo_latentgold(name, case, root = extended_dir)
    expect_identical(as.integer(fit$n_parameters), lg$n_parameters, label = name)
    expect_lt(abs(fit$log_likelihood - lg$log_likelihood), 5e-5 + 2e-4,
              label = paste(name, "LL"))
    expect_false(anyDuplicated(lg$id) > 0L)
    expect_setequal(lg$id, case$data$id)
    rows <- match(lg$id, case$data$id)
    expect_lt(montecarlo_relabel(fit$subject_posteriors[rows, , drop = FALSE],
                                 lg$profiles)$gap, 1e-3, label = paste(name, "posteriors"))
  }))
})
