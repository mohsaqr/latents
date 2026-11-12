# Audit 2026-11-01, slice `profile-fitting` (R/multilpa.R, R/lpa.R,
# R/multilca.R, R/categorical.R, R/extra-indicators.R, R/input-checks.R,
# R/utils.R). One block per repair, then invariant checks of the ordinal and
# negative-binomial likelihood pieces against independent computations.

skip_on_cran()

audit_1101_pf_data <- function(n = 240L) {
  set.seed(11)
  state <- sample.int(2L, n, replace = TRUE)
  data.frame(
    x = stats::rnorm(n, c(0, 2)[state]),
    o = pmin(4, pmax(1, round(stats::rnorm(n, c(1.8, 3.2)[state], 0.8)))),
    k = stats::rpois(n, c(2, 6)[state]),
    b = c("p", "q")[stats::rbinom(n, 1, c(0.2, 0.8)[state]) + 1L],
    # An indicator whose name is also a value of `count_dispersion`.
    equal = stats::rbinom(n, 1, c(0.25, 0.75)[state]))
}

# Fix 1 -----------------------------------------------------------------------
# `ordinal` or `count` that is not a character vector raised an unclassed
# stopifnot() error from .latents_check_extra_arguments().

test_that("a non-character ordinal or count argument is a classed bad argument", {
  data <- audit_1101_pf_data()
  expect_error(lpa(data, c("x", "o"), 2, ordinal = 2, n_starts = 1),
               class = "latents_bad_argument")
  expect_error(lpa(data, c("x", "o"), 2, ordinal = NA_character_, n_starts = 1),
               class = "latents_bad_argument")
  expect_error(lpa(data, c("x", "k"), 2, count = 3, n_starts = 1),
               class = "latents_bad_argument")
  expect_error(lca(data, c("b", "o"), 2, ordinal = NA_character_, n_starts = 1),
               class = "latents_bad_argument")
})

# Fix 2 -----------------------------------------------------------------------
# lca()/multilca()/enumerate_lca() read `extra$count` with `$`, which matches
# partially: with `count_dispersion` (or `count_model`) given and `count`
# absent, its value was taken as count indicator names. An indicator named
# like that value dropped out of the categorical set and was fitted as a
# Gaussian indicator instead, silently changing the model.

test_that("lca() keeps an indicator categorical whatever count_* values are given", {
  data <- audit_1101_pf_data()
  # Since 2026-11-11 a count argument that changes nothing says so.
  expect_warning(
    fit <- lca(data, c("equal", "b", "o"), n_classes = 2, count_dispersion = "equal",
               n_starts = 1, seed = 1),
    class = "latents_unused_argument")
  expect_setequal(names(fit$response_probabilities), c("equal", "b", "o"))
  expect_identical(ncol(fit$means), 0L)
  reference <- lca(data, c("equal", "b", "o"), n_classes = 2, n_starts = 1, seed = 1)
  expect_equal(fit$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  expect_identical(fit$n_parameters, reference$n_parameters)
})

# Fix 3 -----------------------------------------------------------------------
# multilpa() binds argument values to internal symbols before validating, so a
# non-data-frame passed as an expression was named `.data` in the refusal.

test_that("multilpa() names the caller's data, not an internal symbol, when refusing it", {
  data <- audit_1101_pf_data()
  condition <- tryCatch(multilpa(as.matrix(data), c("x", "k"), id = NULL, n_profiles = 2),
                        error = identity)
  expect_s3_class(condition, "latents_bad_data")
  expect_no_match(conditionMessage(condition), "`.data`", fixed = TRUE)
  matrix_data <- as.matrix(data)
  condition <- tryCatch(multilpa(matrix_data, c("x", "k"), id = NULL, n_profiles = 2),
                        error = identity)
  expect_s3_class(condition, "latents_bad_data")
  expect_match(conditionMessage(condition), "`matrix_data`", fixed = TRUE)
})

# Fix 4 -----------------------------------------------------------------------
# The ordinal intercept scores summed with na.rm = TRUE, which would hide an NA
# rather than report it. With missing ordinal cells the scores must still be
# finite and equal the gradient of the observed-data log likelihood (Fisher's
# identity), computed here independently by numerical differentiation.

test_that("ordinal and count scores equal the numerical gradient with missing cells", {
  skip_if_not_installed("numDeriv")
  ns <- asNamespace("latents")
  data <- audit_1101_pf_data(200L)
  data$o[c(3, 17, 40, 41, 99)] <- NA
  data$k[c(5, 60)] <- NA
  extra <- ns$.latents_prepare_extra(data, ordinal = "o", count = "k", missing = "fiml")
  mixing <- c(0.4, 0.6)
  log_likelihood <- function(theta) {
    point <- ns$.latents_extra_decode(theta, extra, 2L)
    joint <- sweep(ns$.latents_extra_log_density(extra, point, nrow(data), 2L), 2L,
                   log(mixing), "+")
    sum(apply(joint, 1L, function(row) max(row) + log(sum(exp(row - max(row))))))
  }
  start <- list(ordinal_intercepts = list(c(0.3, -0.2, -0.9)),
                ordinal_locations = matrix(c(-1.1, 0), 2L),
                count_means = matrix(c(2.2, 5.5), 2L))
  theta <- unname(ns$.latents_extra_coordinates(start, extra, 2L, "unconstrained"))
  point <- ns$.latents_extra_decode(theta, extra, 2L)
  joint <- sweep(ns$.latents_extra_log_density(extra, point, nrow(data), 2L), 2L,
                 log(mixing), "+")
  posteriors <- exp(joint - apply(joint, 1L, max))
  posteriors <- posteriors / rowSums(posteriors)
  scores <- ns$.latents_extra_scores(extra, posteriors, point)
  expect_false(anyNA(scores))
  expect_equal(as.vector(scores), numDeriv::grad(log_likelihood, theta),
               tolerance = 1e-6)
})

# Invariants of the likelihood pieces ------------------------------------------

test_that("ordinal category probabilities are the adjacent-category logit, summing to one", {
  ns <- asNamespace("latents")
  intercepts <- c(0.5, -0.2, 0.1)
  locations <- c(1.3, -0.4, 0)
  probability <- exp(ns$.latents_ordinal_log_probabilities(intercepts, locations))
  hand <- t(vapply(locations, function(eta) {
    linear <- c(0, intercepts) + eta * (0:3)
    exp(linear) / sum(exp(linear))
  }, numeric(4)))
  expect_equal(probability, hand, tolerance = 1e-14)
  expect_equal(rowSums(probability), rep(1, 3), tolerance = 1e-14)
})

test_that("the ordinal M-step reaches the maximum a general optimizer finds", {
  ns <- asNamespace("latents")
  counts <- matrix(c(10, 5, 2, 3, 6, 9, 4, 4, 4), 3L, byrow = TRUE)
  newton <- ns$.latents_ordinal_maximize(counts, c(0, 0), c(0, 0, 0))
  objective <- function(theta) {
    -sum(counts * ns$.latents_ordinal_log_probabilities(theta[1:2], c(theta[3:4], 0)))
  }
  general <- stats::optim(c(0, 0, 0, 0), objective, method = "BFGS",
                          control = list(reltol = 1e-14, maxit = 1000))
  expect_equal(c(newton$intercepts, newton$locations[1:2]), general$par, tolerance = 1e-5)
  expect_equal(newton$locations[[3]], 0)
  expect_lte(objective(c(newton$intercepts, newton$locations[1:2])), general$value + 1e-9)
})

test_that("negative-binomial dispersion derivatives match finite differences of dnbinom", {
  skip_if_not_installed("numDeriv")
  ns <- asNamespace("latents")
  cases <- expand.grid(y = c(0, 3, 40), mu = c(0.7, 6, 80), alpha = c(0.02, 0.4, 3))
  derivatives <- ns$.latents_negative_binomial_dispersion_derivatives(
    cases$y, cbind(cases$mu), cases$alpha)
  invisible(lapply(seq_len(nrow(cases)), function(row) {
    density <- function(log_alpha) {
      stats::dnbinom(cases$y[row], size = exp(-log_alpha), mu = cases$mu[row], log = TRUE)
    }
    at <- log(cases$alpha[row])
    expect_equal(derivatives$score[row, 1L], numDeriv::grad(density, at),
                 tolerance = 1e-6)
    expect_equal(derivatives$curvature[row, 1L],
                 drop(numDeriv::hessian(density, at)), tolerance = 1e-5)
  }))
})

test_that("integer sampling weights fit ordinal and count indicators as repeated rows", {
  data <- audit_1101_pf_data(200L)
  set.seed(12)
  data$w <- sample(1:2, nrow(data), replace = TRUE)
  repeated <- data[rep(seq_len(nrow(data)), data$w), ]
  weighted <- lpa(data, c("x", "o", "k"), 2, ordinal = "o", count = "k", weights = "w",
                  n_starts = 2, seed = 1, tol = 1e-12, max_iter = 5000)
  plain <- lpa(repeated, c("x", "o", "k"), 2, ordinal = "o", count = "k",
               n_starts = 2, seed = 1, tol = 1e-12, max_iter = 5000)
  order_of <- function(fit) order(fit$means[, "x"])
  expect_equal(unname(weighted$count_means[order_of(weighted), ]),
               unname(plain$count_means[order_of(plain), ]), tolerance = 1e-5)
  located <- function(fit) {
    location <- fit$ordinal_locations[order_of(fit), "o"]
    location - location[[2]]
  }
  expect_equal(unname(located(weighted)), unname(located(plain)), tolerance = 1e-5)
  expect_equal(weighted$log_likelihood * sum(data$w) / nrow(data), plain$log_likelihood,
               tolerance = 1e-8)
})
