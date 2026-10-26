# Pure control/selection checks: no external programs or Monte Carlo kit runs.

search_protocol_environment <- function(fit_one) {
  environment <- new.env(parent = asNamespace("latents"))
  sys.source(file.path("features", "direct-effects", "helpers.R"), envir = environment)
  environment$lpa <- fit_one
  environment
}

search_protocol_case <- function(n = 50L, model = "VVI") {
  list(name = "synthetic", data = data.frame(y = seq_len(n)),
       spec = list(vars = "y", k = 2L, model = model, direct_effects = "y ~ z"))
}

test_that("both search controls are fixed before reading any reference", {
  environment <- search_protocol_environment(function(...) stop("No fitting in this test"))
  small <- environment$de_search_controls(search_protocol_case(), 1:3, 20L)
  large <- environment$de_search_controls(search_protocol_case(10000L), 1:3, 20L)
  expect_setequal(small$search, c("kmeans", "emEM"))
  expect_equal(nrow(small), 6L)
  expect_true(all(small$n_starts == 20L))
  expect_equal(small$n_screen[small$search == "emEM"], rep(100L, 3))
  expect_true(all(is.na(small$n_screen[small$search == "kmeans"])))
  expect_equal(nrow(large), 4L)
  expect_equal(unique(large$seed), 1:2)
  expect_true(all(large$n_starts == 10L))
  expect_equal(large$n_screen[large$search == "emEM"], rep(50L, 2))
  shared <- c("tol", "max_iter", "min_variance", "boundary_policy")
  expect_true(all(vapply(small[shared], function(x) length(unique(x)) == 1L, logical(1))))
  expect_true(all(vapply(shared, function(name) identical(small[[name]][1],
                                                         large[[name]][1]), logical(1))))
  equal <- environment$de_search_controls(search_protocol_case(model = "EEI"))
  expect_true(all(equal$default_search == "emEM"))
  expect_true(all(small$default_search == "kmeans"))
})

test_that("selection reports defaults and rejects boundary or unconverged higher attempts", {
  seen <- list()
  fit_one <- function(...) {
    arguments <- list(...)
    seen[[length(seen) + 1L]] <<- arguments
    search <- arguments$search
    seed <- arguments$seed
    list(log_likelihood = if (search == "kmeans") -100 else c(-90, -10, -20)[seed],
         converged = !(search == "emEM" && seed == 3L),
         boundary = search == "emEM" && seed == 2L)
  }
  environment <- search_protocol_environment(fit_one)
  fit <- environment$de_fit(search_protocol_case())
  expect_equal(length(seen), 6L)
  expect_equal(fit$log_likelihood, -90)
  expect_equal(fit$de_default_fit$log_likelihood, -100)
  expect_equal(fit$de_default_seed_log_likelihoods, rep(-100, 3))
  expect_equal(fit$seed_log_likelihoods, c(-90, -100, -100))
  attempts <- fit$de_search_attempts
  expect_equal(nrow(attempts), 6L)
  expect_equal(sum(attempts$selected), 1L)
  expect_equal(sum(attempts$is_default), 3L)
  expect_false(any(attempts$admissible & (attempts$boundary | !attempts$converged)))
  expect_true(all(vapply(seen, function(x) identical(x$min_variance, 1e-6) &&
                          identical(x$boundary, "allow"), logical(1))))
  expect_true(all(vapply(seen, function(x) if (x$search == "emEM") {
    identical(x$n_screen, 100L)
  } else is.null(x$n_screen), logical(1))))
})

test_that("failed attempts remain visible and no admissible fit fails validation", {
  fit_one <- function(...) {
    arguments <- list(...)
    if (arguments$search == "kmeans") stop("recorded attempt failure")
    list(log_likelihood = -90, converged = TRUE, boundary = FALSE)
  }
  environment <- search_protocol_environment(fit_one)
  fit <- environment$de_fit(search_protocol_case())
  expect_null(fit$de_default_fit)
  expect_true(all(is.na(fit$de_default_seed_log_likelihoods)))
  expect_equal(sum(!is.na(fit$de_search_attempts$error)), 3L)
  expect_true(all(grepl("recorded attempt failure",
                       fit$de_search_attempts$error[fit$de_search_attempts$is_default])))
  environment$lpa <- function(...) list(log_likelihood = -10, converged = TRUE,
                                         boundary = TRUE)
  expect_error(environment$de_fit(search_protocol_case()), "neither configured search")
})
