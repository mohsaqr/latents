# The shared estimation kernel (validation/ENGINE_DESIGN.md): the EM driver,
# the start runner, the quasi-Newton finish, the measurement block, the
# Markov structure and the inference service. Every engine is checked
# against its recorded outputs by equivalence/engine-golden/; these tests
# pin the kernel's own contract.

kernel_toy_em <- function(rate = 0.5) {
  # A one-parameter fixed point theta -> rate * theta + 1, whose "likelihood"
  # -(theta - 1 / (1 - rate))^2 rises monotonically along the iteration.
  target <- 1 / (1 - rate)
  list(evaluate = function(state) list(log_likelihood = -(state - target)^2),
       maximize = function(state, expectation) rate * state + 1,
       target = target)
}

test_that("the EM driver converges, counts steps and records the history", {
  toy <- kernel_toy_em()
  fit <- .latents_em(0, toy$evaluate, toy$maximize, max_iter = 500, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$state, toy$target, tolerance = 1e-5)
  expect_length(fit$history, fit$iterations + 1L)
  # Invariant: the likelihood never falls along a monotone EM path.
  expect_true(all(diff(fit$history) >= 0))
  limited <- .latents_em(0, toy$evaluate, toy$maximize, max_iter = 3, tol = 1e-12)
  expect_false(limited$converged)
  expect_identical(limited$iterations, 3L)
  none <- .latents_em(0, toy$evaluate, toy$maximize, max_iter = 0, tol = 1e-12)
  expect_identical(none$iterations, 0L)
  expect_identical(none$state, 0)
})

test_that("a decrease beyond rounding raises the engine's condition", {
  evaluate <- function(state) list(log_likelihood = -state)
  maximize <- function(state, expectation) state + 1
  expect_error(.latents_em(0, evaluate, maximize, 10, 1e-8),
               "decreased beyond numerical roundoff")
  expect_error(.latents_em(0, evaluate, maximize, 10, 1e-8,
                           .latents_em_settings(decrease_condition = "classed")),
               class = "latents_em_decrease")
  expect_error(.latents_em(0, evaluate, maximize, 10, 1e-8,
                           .latents_em_settings(decrease_condition = "Engine message.")),
               "Engine message.", fixed = TRUE)
  unguarded <- .latents_em(0, evaluate, maximize, 4, 1e-8,
                           .latents_em_settings(guard_decrease = FALSE))
  expect_identical(unguarded$iterations, 4L)
})

test_that("a non-finite step ends a run that asks for it, keeping that state", {
  evaluate <- function(state) list(log_likelihood = if (state > 2) NaN else state)
  maximize <- function(state, expectation) state + 1
  fit <- .latents_em(0, evaluate, maximize, 10, 1e-8,
                     .latents_em_settings(stop_on_nonfinite = TRUE))
  expect_identical(fit$state, 3)
  expect_identical(fit$expectation$log_likelihood, 2)
  expect_false(fit$converged)
})

test_that("an unsolved inner step blocks convergence until the stall limit", {
  toy <- kernel_toy_em()
  unsolved <- function(state, expectation) {
    structure(toy$maximize(state, expectation), solved = FALSE)
  }
  fit <- .latents_em(toy$target, toy$evaluate, unsolved, max_iter = 100, tol = 1e-8,
                     .latents_em_settings(max_stalled = 5L))
  expect_false(fit$converged)
  expect_true(fit$settled)
  expect_identical(fit$iterations, 5L)
})

test_that("SQUAREM reaches the same fixed point in fewer evaluations", {
  toy <- kernel_toy_em(rate = 0.95)
  accelerator <- list(flatten = function(state) state,
                      restore = function(state, values) values,
                      valid = function(state) is.finite(state))
  plain <- .latents_em(0, toy$evaluate, toy$maximize, 5000, 1e-14)
  fast <- .latents_em(0, toy$evaluate, toy$maximize, 5000, 1e-14,
                      .latents_em_settings(accelerate = accelerator))
  expect_equal(fast$state, plain$state, tolerance = 1e-6)
  expect_lt(fast$iterations, plain$iterations)
})

test_that("the start runner records failures and refuses when all fail", {
  run <- function(start) {
    if (start == 2L) stop("bad start")
    list(expectation = list(log_likelihood = -start), converged = TRUE)
  }
  starts <- .latents_run_starts(3L, run, "likelihood")
  expect_identical(starts$valid, c(TRUE, FALSE, TRUE))
  expect_identical(starts$scores, c(-1, -Inf, -3))
  expect_identical(starts$best_start, 1L)
  expect_identical(starts$attempts[[2L]]$error, "bad start")
  expect_error(.latents_run_starts(2L, function(start) stop("no"), "likelihood"),
               class = "latents_all_starts_failed")
})

test_that("the quasi-Newton finish keeps a point only when it improves", {
  value <- function(v) sum((v - c(1, -2))^2)
  gradient <- function(v) 2 * (v - c(1, -2))
  finish <- .latents_quasi_newton(c(0, 0), c(-Inf, -Inf), value, gradient, 1e-10,
                                  list(maxit = 100L, factr = 10, pgtol = 0))
  expect_true(finish$improved)
  expect_equal(finish$par, c(1, -2), tolerance = 1e-6)
  bounded <- .latents_quasi_newton(c(0, 0), c(-Inf, 0), value, gradient, 1e-10,
                                   list(maxit = 100L, factr = 10, pgtol = 0))
  expect_equal(bounded$par, c(1, 0), tolerance = 1e-6)
  # Started at the optimum it stays there.
  settled <- .latents_quasi_newton(c(1, -2), c(-Inf, -Inf), value, gradient, 1e-10,
                                   list(maxit = 100L, factr = 10, pgtol = 0))
  expect_true(settled$converged)
  expect_identical(settled$par, c(1, -2))
})

test_that("the measurement block's two normalizers agree to rounding", {
  x <- matrix(c(-1, 0.5, 2, 3, -0.2, 1.1), 3L)
  means <- matrix(c(0, 1, 2, -1), 2L)
  variances <- matrix(c(1, 2, 0.5, 3), 2L)
  separate <- .latents_gaussian_log_density(x, means, variances, "separate")
  joint <- .latents_gaussian_log_density(x, means, variances, "joint")
  expect_equal(separate, joint, tolerance = 1e-14)
  by_hand <- vapply(1:2, function(k) {
    rowSums(stats::dnorm(x, rep(means[k, ], each = 3L),
                         rep(sqrt(variances[k, ]), each = 3L), log = TRUE))
  }, numeric(3))
  expect_equal(separate, by_hand, tolerance = 1e-12)
  # A single row keeps its matrix shape.
  expect_identical(dim(.latents_gaussian_log_density(x[1L, , drop = FALSE],
                                                     means, variances)), c(1L, 2L))
  offset <- .latents_measurement_log_density(
    x, list(means = means, variances = variances))
  expect_equal(offset$log_density + offset$offset, separate, tolerance = 1e-14)
  expect_identical(apply(offset$log_density, 1L, max), rep(0, 3L))
})

test_that("the homogeneous chain is the shared Markov structure exactly", {
  layout <- list(slot = matrix(1:6, 2L), within = matrix(TRUE, 2L, 3L),
                 n_occasions = 3L)
  emission <- lapply(1:3, function(t) matrix(c(-1, -2, -0.5, -3) * t, 2L))
  initial <- c(0.3, 0.7)
  transition <- matrix(c(0.8, 0.3, 0.2, 0.7), 2L)
  homogeneous <- .multilpa_forward_backward(emission, layout, initial, transition)
  # By hand: the forward recursion over every path.
  paths <- as.matrix(expand.grid(1:2, 1:2, 1:2))
  by_group <- vapply(1:2, function(g) {
    log(sum(apply(paths, 1L, function(p) {
      initial[p[1L]] * transition[p[1L], p[2L]] * transition[p[2L], p[3L]] *
        exp(emission[[1L]][g, p[1L]] + emission[[2L]][g, p[2L]] +
              emission[[3L]][g, p[3L]])
    })))
  }, numeric(1))
  expect_equal(homogeneous$log_scaled, by_group, tolerance = 1e-12)
  expect_error(.multilpa_forward_backward(emission, layout, initial,
                                          transition * c(1, 0)))
})

test_that("the inference service refuses what it cannot invert", {
  expect_error(.inference_invert(diag(c(1, -1))), class = "latents_singular_information")
  expect_error(.inference_invert(-diag(2), "ratio_le_eps"),
               class = "latents_singular_information")
  expect_false(.inference_positive_definite(diag(c(1, 1e-12)), "min_le_eps_maxabs"))
  expect_equal(.inference_invert(diag(c(2, 4))), diag(c(0.5, 0.25)))
  expect_error(.inference_require_rank(matrix(1, 3L, 2L), "too few"),
               class = "latents_too_few_groups")
  expect_identical(.inference_critical(0.95, "upper"),
                   .inference_critical(0.95, "one_minus"))
  # Invariant: a probability vector's softmax Jacobian has zero column sums.
  jacobian <- .inference_simplex_jacobian(c(0.2, 0.3, 0.5), "first")
  expect_equal(colSums(jacobian), c(0, 0), tolerance = 1e-15)
  covariance <- matrix(c(2, 0.3, 0.3, 1), 2L)
  linear <- matrix(c(1, 0, 1, 1, 0, 2), 3L)
  expect_equal(.inference_delta_se(linear, covariance, "quadratic"),
               .inference_delta_se(linear, covariance, "rowsums"), tolerance = 1e-14)
  expect_equal(.inference_delta_covariance(linear, covariance),
               linear %*% covariance %*% t(linear))
})

test_that("the score information matches a known Hessian", {
  # log L = -(a^2 + a b + 2 b^2): information [[2, 1], [1, 4]].
  score <- function(theta) c(-(2 * theta[1L] + theta[2L]), -(theta[1L] + 4 * theta[2L]))
  theta <- c(a = 0.3, b = -0.2)
  expected <- matrix(c(2, 1, 1, 4), 2L, dimnames = list(c("a", "b"), c("a", "b")))
  expect_equal(.inference_score_information(score, theta, 1e-5, "central", "plus1"),
               expected, tolerance = 1e-8)
  expect_equal(.inference_score_information(score, theta, 1e-4, "five_point", "max1"),
               expected, tolerance = 1e-8)
})
