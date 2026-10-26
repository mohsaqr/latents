# Verify protocol bookkeeping without running the complete external kit.

ji_protocol_environment <- function() {
  environment <- new.env(parent = asNamespace("latents"))
  sys.source(file.path("features", "joint-items", "helpers.R"), envir = environment)
  environment
}

test_that("joint verification controls never loosen the original requested precision", {
  environment <- ji_protocol_environment()
  expect_identical(environment$ji_verification_controls(1e-10, 1000L),
                   list(n_starts = 1L, tol = 1e-13, max_iter = 5000L))
  expect_identical(environment$ji_verification_controls(1e-14, 200000L),
                   list(n_starts = 1L, tol = 1e-14, max_iter = 200000L))
  seen <- NULL
  own <- list(token = 12)
  environment$starting_values <- function(fit) {
    expect_identical(fit, own)
    list(token = fit$token)
  }
  environment$ji_fit_once <- function(case, seed, n_starts, tol, max_iter,
                                      local_dependence, start) {
    seen <<- list(seed = seed, starts = n_starts, tol = tol, max_iter = max_iter,
                   pairs = local_dependence, start = start)
    own
  }
  expect_identical(environment$ji_refine(list(), own, 7L, local_dependence = FALSE), own)
  expect_identical(seen, list(seed = 7L, starts = 1L, tol = 1e-13,
                             max_iter = 5000L, pairs = FALSE, start = list(token = 12)))
})

test_that("joint verification retains all original seeds rather than replacing their evidence", {
  environment <- ji_protocol_environment()
  calls <- list()
  environment$ji_fit_once <- function(case, seed, n_starts, tol, max_iter,
                                      local_dependence) {
    calls[[length(calls) + 1L]] <<- list(seed = seed, starts = n_starts,
                                       tol = tol, max_iter = max_iter)
    list(log_likelihood = -100 - seed, converged = TRUE)
  }
  environment$ji_refine <- function(case, fit, seed, tol, max_iter, local_dependence) {
    expect_equal(fit$log_likelihood, -100 - seed)
    fit$log_likelihood <- fit$log_likelihood + .25
    fit
  }
  fit <- environment$ji_fit(list(), n_starts = 3L, seeds = c(7L, 9L))
  expect_equal(length(calls), 2L)
  expect_equal(vapply(calls, `[[`, integer(1), "seed"), c(7L, 9L))
  expect_equal(vapply(fit$ji_original_fits, `[[`, numeric(1), "log_likelihood"),
               c(-107, -109))
  expect_equal(fit$ji_seed_maxima, c(-106.75, -108.75))
  expect_identical(fit$ji_seeds, c(7L, 9L))
  expect_equal(fit$ji_original_controls$n_starts, 3L)
  expect_equal(fit$ji_original_controls$tol, 1e-10)
  expect_equal(fit$ji_verification_controls$tol, 1e-13)
})

test_that("actual joint seed fits resume from their own parameters", {
  environment <- ji_protocol_environment()
  case <- readRDS(file.path("features", "joint-items", "kit", "cases.rds"))$bin_small_01
  fit <- environment$ji_fit(case, n_starts = 2L, seeds = 1:2)
  original <- vapply(fit$ji_original_fits, `[[`, numeric(1), "log_likelihood")
  verified <- vapply(fit$ji_verified_fits, `[[`, numeric(1), "log_likelihood")
  expect_true(all(verified >= original - 1e-8))
  expect_identical(verified, fit$ji_seed_maxima)
  expect_true(all(vapply(fit$ji_verified_fits, function(x) isTRUE(x$converged), logical(1))))
  expect_equal(length(fit$ji_original_fits), 2L)
  expect_equal(fit$ji_original_controls$tol, 1e-10)
  expect_equal(fit$ji_verification_controls$n_starts, 1L)
})
