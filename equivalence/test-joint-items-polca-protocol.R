# Declared reference budgets and diagnostics; no complete external-kit rerun.

polca_protocol_environment <- function() {
  environment <- new.env(parent = asNamespace("latents"))
  sys.source("helper-montecarlo.R", envir = environment)
  invisible(lapply(c("helpers.R", "reference.R", "compare.R"), function(file) {
    sys.source(file.path("features", "joint-items", file), envir = environment)
  }))
  environment
}

test_that("all eligible poLCA cases use the same two declared iteration budgets", {
  environment <- polca_protocol_environment()
  cases <- readRDS(file.path("features", "joint-items", "kit", "cases.rds"))
  eligible <- Filter(function(case) is.null(case$spec$h) && is.null(case$spec$weights) &&
                       !identical(case$config, "huge"), cases)
  expect_length(eligible, 49L)
  calls <- list()
  environment$ji_polca_fit <- function(formula, frame, k, seed, nrep, maxiter, tol) {
    calls[[length(calls) + 1L]] <<- list(seed = seed, nrep = nrep,
                                       maxiter = maxiter, tol = tol)
    list(llik = if (maxiter == 20000L) -100 else -99,
         numiter = if (maxiter == 20000L) 20000L else 50000L,
         attempts = rep(if (maxiter == 20000L) -100 else -99, nrep))
  }
  environment$ji_polca_rows <- function(case, fit, polca, polish) {
    list(log_likelihood = polca$llik, own = NULL,
         rows = data.frame(case = case$name, log_likelihood = polca$llik))
  }
  results <- lapply(eligible, function(case) environment$ji_polca(case, list()))
  expect_length(calls, 98L)
  expect_equal(vapply(calls, `[[`, integer(1), "maxiter"), rep(c(20000L, 200000L), 49L))
  expect_true(all(vapply(calls, function(x) identical(x$nrep, 10L) &&
                          identical(x$tol, 1e-12), logical(1))))
  expected_seed <- rep(vapply(eligible, `[[`, numeric(1), "seed"), each = 2L)
  expect_equal(unname(vapply(calls, `[[`, numeric(1), "seed")), unname(expected_seed))
  expect_true(all(vapply(results, function(x) x$log_likelihood == -99 &&
                          x$original$log_likelihood == -100, logical(1))))
  expect_true(all(vapply(results, function(x) identical(names(x$fits),
                          c("original", "verification")), logical(1))))
  expect_true(all(vapply(results, function(x) nrow(x$attempts) == 20L &&
                          identical(x$diagnostics$iteration_cap_reached,
                                    c(TRUE, FALSE)), logical(1))))
})

test_that("poLCA resets the same case seed for both stock calls and restores caller RNG", {
  skip_if_not_installed("poLCA")
  environment <- polca_protocol_environment()
  draws <- list()
  testthat::local_mocked_bindings(poLCA = function(...) {
    arguments <- list(...)
    draws[[length(draws) + 1L]] <<- stats::runif(8L)
    list(llik = -50, numiter = 2L, attempts = rep(-50, arguments$nrep))
  }, .package = "poLCA")
  environment$ji_polca_rows <- function(case, fit, polca, polish) {
    list(log_likelihood = polca$llik, rows = data.frame(log_likelihood = polca$llik))
  }
  case <- list(name = "seed_probe", seed = 175L,
               data = data.frame(a = rep(1:2, 5L), b = rep(2:1, 5L)),
               spec = list(vars = c("a", "b"), pairs = character(), k = 2L))
  withr::local_seed(2059L)
  before <- .Random.seed
  result <- environment$ji_polca(case, list(), nrep = 3L)
  expect_length(draws, 2L)
  expect_identical(draws[[1L]], draws[[2L]])
  expect_identical(.Random.seed, before)
  expect_equal(result$controls$nrep, c(3L, 3L))
  expect_equal(nrow(result$attempts), 6L)
})

test_that("an ordinary one-class stock poLCA fit is unchanged by the larger cap", {
  skip_if_not_installed("poLCA")
  environment <- polca_protocol_environment()
  data <- expand.grid(a = 1:2, b = 1:2, c = 1:2)
  data <- data[rep(seq_len(nrow(data)), seq_len(nrow(data)) + 5L), ]
  case <- list(name = "ordinary_one_class", seed = 267L, data = data,
               spec = list(vars = c("a", "b", "c"), pairs = "a:b", k = 1L))
  fit <- lca(data, case$spec$vars, 1L, local_dependence = "a:b",
              n_starts = 1L, seed = 1L)
  result <- environment$ji_polca(case, fit, nrep = 1L)
  expect_identical(result$fits$original$llik, result$fits$verification$llik)
  expect_identical(result$fits$original$probs, result$fits$verification$probs)
  expect_identical(result$rows, result$original$rows)
  expect_true(all(result$rows$difference <= result$rows$tolerance))
  expect_false(any(result$diagnostics$iteration_cap_reached))
  expect_equal(result$diagnostics$numiter, c(1L, 1L))
  expect_true(all(result$attempts$closed_form))
  expect_true(all(is.na(result$attempts$attempt)))
  expect_equal(result$attempts$log_likelihood, result$diagnostics$log_likelihood)
})
