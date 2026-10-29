# Review 2026-10-26, round 2, area kernel: regression tests for K2-1 to K2-7.
# Each test failed on the code before its fix (checked by running this file
# against a copy of the package with the old lines restored).
# The comparisons with independent references are in tests/equivalence/test-review-1026-kernel-r2.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

# K2-1 -----------------------------------------------------------------------

test_that("K2-1: the VEI shape sums carry no variance bound", {
  # Variances (1e-4, 1e14) and (3e-4, 2e14), all far above the 1e-6 bound.
  # The unnormalised shape sum of the first indicator is about 1e-7, which the
  # old code floored at min_variance: variances 12.98x too large, objective
  # worse by 183.27.
  weights <- c(60, 40)
  truth <- rbind(c(1e-4, 1e14), c(3e-4, 2e14))
  diagonals <- truth * weights
  free_vei <- function(scatter, w) {
    normalise <- function(s) s / exp(mean(log(s)))
    shape <- normalise(colSums(scatter))
    volumes <- NULL
    invisible(lapply(seq_len(5000), function(i) {
      volumes <<- rowSums(sweep(scatter, 2L, shape, "/")) / (w * ncol(scatter))
      shape <<- normalise(colSums(scatter / volumes))
    }))
    volumes <- rowSums(sweep(scatter, 2L, shape, "/")) / (w * ncol(scatter))
    outer(volumes, shape)
  }
  objective <- function(v) sum(weights * rowSums(log(v)) + rowSums(diagonals / v))
  fitted <- .multilpa_structure_variances(diagonals, weights, "VEI", 1e-6)
  reference <- free_vei(diagonals, weights)
  expect_lt(max(abs(fitted / reference - 1)), 1e-8)
  expect_lt(abs(objective(fitted) - objective(reference)), 1e-8)
})

# K2-4 -----------------------------------------------------------------------

test_that("K2-4: the EVE and VVE M-steps are equivariant at squared units of 1e-18", {
  # The old absolute floor .Machine$double.eps on the rotated scatter gave
  # relative errors of 0.531 (EVE) and 0.491 (VVE) at c^2 = 1e-18.
  withr::local_seed(2)
  make_scatter <- function(n) {
    a <- matrix(stats::rnorm(n * 3), n) %*%
      matrix(c(1, .6, .2, 0, 1, .4, 0, 0, .5), 3)
    crossprod(a)
  }
  scatter <- list(make_scatter(40), make_scatter(60) * 2)
  weights <- c(40, 60)
  invisible(lapply(c(TRUE, FALSE), function(equal_volume) {
    base <- .multilpa_eve_covariances(scatter, weights, 3, 1000, 1e-12, equal_volume)
    scaled <- .multilpa_eve_covariances(lapply(scatter, `*`, 1e-18), weights, 3, 1000,
                                        1e-12, equal_volume)
    errors <- vapply(1:2, function(k) {
      max(abs(scaled[[k]] / 1e-18 - base[[k]])) / max(abs(base[[k]]))
    }, numeric(1))
    expect_lt(max(errors), 1e-10)
  }))
})

# K2-6 -----------------------------------------------------------------------

test_that("K2-6: with a stop scale, the EM stop does not move with the units", {
  # Rescaling the data by c adds -n d log(c) to every log likelihood. The
  # kernel's rule with `stop_scale` reads only changes, so a constant added to
  # the likelihood leaves the path, and the step at which it stops, alone.
  evaluate_with <- function(constant) {
    function(state) list(log_likelihood = constant - 1000 * sum((state$p - 2)^2))
  }
  maximize <- function(state, expectation) list(p = state$p + 0.3 * (2 - state$p))
  run <- function(constant, settings) {
    .latents_em(list(p = c(0, 5)), evaluate_with(constant), maximize,
                max_iter = 10000, tol = 1e-8, settings = settings)
  }
  scaled <- .latents_em_settings(stop_scale = 200)
  near <- run(0, scaled)
  far <- run(-1e9, scaled)
  expect_true(near$converged && far$converged)
  expect_identical(near$iterations, far$iterations)
  expect_equal(near$state$p, far$state$p, tolerance = 0)
  # The older rule, kept for an engine that supplies no scale, stops earlier
  # when |log L| is large.
  legacy_far <- run(-1e9, .latents_em_settings())
  expect_lt(legacy_far$iterations, far$iterations)
  expect_error(.latents_em_settings(stop_scale = -1))
})

# K2-7 -----------------------------------------------------------------------

test_that("K2-7: SQUAREM counts the M-steps it keeps", {
  # An extrapolation that is never valid: every cycle falls back to two plain
  # steps. The old count added three per cycle (30 for 20 M-steps).
  calls <- 0L
  evaluate <- function(state) list(log_likelihood = -sum((state$p - 2)^2))
  maximize <- function(state, expectation) {
    calls <<- calls + 1L
    list(p = state$p + 0.5 * (2 - state$p))
  }
  never_valid <- list(flatten = function(s) s$p, restore = function(s, v) list(p = v),
                      valid = function(s) FALSE)
  fit <- .latents_em(list(p = c(0, 5)), evaluate, maximize, max_iter = 30, tol = 1e-14,
                     settings = .latents_em_settings(accelerate = never_valid))
  expect_identical(fit$iterations, calls)
  expect_identical(length(fit$history) - 1L, fit$iterations)
  # Accepted extrapolations keep three M-steps each, and the count still
  # matches the history.
  calls <- 0L
  always_valid <- list(flatten = function(s) s$p, restore = function(s, v) list(p = v),
                       valid = function(s) TRUE)
  accelerated <- .latents_em(list(p = c(0, 5)), evaluate, maximize, max_iter = 30,
                             tol = 1e-14,
                             settings = .latents_em_settings(accelerate = always_valid))
  expect_identical(length(accelerated$history) - 1L, accelerated$iterations)
  expect_lte(accelerated$iterations, calls)
})

# K2-6, wired into the engines ----------------------------------------------

test_that("K2-6: a multilpa() fit stops at the same EM step whatever the units", {
  # Overlapping two-profile VVI data (n = 3000 in the review's measurement;
  # n = 1000 here). Measured after the fix, distance at stop from a 1e-15
  # reference in Wald standard errors:
  #   plain EM:  0.048 at units x1, x1e-6 and x1e6 alike (126 steps each);
  #              before, 0.062 (x1), 0.082 (x0.01), 0.137 (x1e3), 0.189 (x1e6).
  #   SQUAREM:   0.0063 (x1), 0.018 (x1e-6), 0.00027 (x1e6); the remaining
  #              spread is the extrapolation steplength, whose norm mixes
  #              coordinates of different units, not the stopping rule.
  withr::local_seed(42)
  n <- 1000
  z <- stats::rbinom(n, 1, 0.4)
  x <- matrix(stats::rnorm(n * 3), n, 3) + outer(z, c(0.7, 0.5, 0.6))
  x[z == 1, ] <- x[z == 1, ] * 1.4
  data <- as.data.frame(x)
  names(data) <- c("a", "b", "c")
  fit_at <- function(units) {
    suppressMessages(suppressWarnings(multilpa(
      data * units, names(data), id = NULL, n_profiles = 2, model = "VVI",
      n_starts = 1, seed = 1, acceleration = "none", max_iter = 100000)))
  }
  base <- fit_at(1)
  invisible(lapply(c(1e-6, 1e6), function(units) {
    scaled <- fit_at(units)
    expect_identical(scaled$iterations, base$iterations)
    expect_equal(unname(scaled$means) / units, unname(base$means), tolerance = 1e-7)
    expect_equal(scaled$log_likelihood + n * 3 * log(units), base$log_likelihood,
                 tolerance = 1e-9)
  }))
})

# K2-2 -----------------------------------------------------------------------

test_that("K2-2: the default variance bound is a millionth of the smallest indicator variance", {
  x <- cbind(a = c(1, 2, 4, 8, NA), b = c(0.1, 0.3, 0.2, 0.25, 0.15))
  expect_equal(.latents_resolve_min_variance(NULL, x),
               1e-6 * min(stats::var(x[, "a"], na.rm = TRUE), stats::var(x[, "b"])))
  # Rows of weight zero do not enter; a number is kept as given.
  weights <- c(1, 1, 1, 0, 1)
  expect_equal(.latents_resolve_min_variance(NULL, x, weights),
               1e-6 * min(stats::var(c(1, 2, 4)), stats::var(c(0.1, 0.3, 0.2, 0.15))))
  expect_identical(.latents_resolve_min_variance(3e-7, x), 3e-7)
  # It moves with the units, also at extreme magnitudes.
  expect_equal(.latents_resolve_min_variance(NULL, x * 1e-100) / 1e-200,
               .latents_resolve_min_variance(NULL, x), tolerance = 1e-12)
  expect_equal(.latents_resolve_min_variance(NULL, x + 1e12),
               .latents_resolve_min_variance(NULL, x), tolerance = 1e-6)
})

test_that("K2-2: iris in metres fits as iris in centimetres", {
  # Before: VVV in metres gave a log likelihood 1.00 below the centimetre
  # fit after the change of units, VVI at x0.001 414.4 below, EEE 464.9, each
  # with only `latents_boundary`. After: within 1e-8 at x0.01, 1e-3, 1e-4.
  data <- iris[, 1:4]
  names(data) <- c("a", "b", "c", "d")
  fit_at <- function(units, model, ...) {
    suppressMessages(multilpa(data * units, names(data), id = NULL, n_profiles = 3,
                              model = model, n_starts = 5, seed = 1, ...))
  }
  invisible(lapply(c("VVV", "VVI"), function(model) {
    base <- fit_at(1, model)
    expect_equal(base$min_variance,
                 1e-6 * min(vapply(data, stats::var, numeric(1))))
    metres <- expect_no_warning(fit_at(0.01, model))
    expect_equal(metres$log_likelihood + nrow(data) * 4 * log(0.01),
                 base$log_likelihood, tolerance = 1e-8)
    expect_equal(metres$min_variance, base$min_variance * 1e-4)
  }))
  # An explicit number keeps its absolute meaning (and its old behaviour).
  absolute <- suppressWarnings(fit_at(0.001, "VVI", min_variance = 1e-6))
  expect_identical(absolute$min_variance, 1e-6)
  expect_true(absolute$boundary)
})

test_that("K2-2: lta(), the additive and the cross-level families resolve the same default", {
  withr::local_seed(3)
  long <- data.frame(id = rep(seq_len(40), each = 3), time = rep(1:3, 40),
                     y1 = stats::rnorm(120, 0, 0.01), y2 = stats::rnorm(120, 5, 2))
  long$y1 <- long$y1 + rep(stats::rbinom(40, 1, 0.5), each = 3) * 0.02
  smallest <- 1e-6 * stats::var(long$y1)
  quiet <- function(expr) suppressMessages(suppressWarnings(expr))
  transitions <- quiet(lta(long, c("y1", "y2"), "id", n_profiles = 2, time = "time",
                           n_starts = 2, seed = 1))
  expect_equal(transitions$min_variance, smallest)
  additive <- quiet(multilpa(long, c("y1", "y2"), "id", n_group_classes = 2,
                             family = "additive", n_starts = 2, seed = 1))
  expect_equal(additive$min_variance, smallest)
  cross <- quiet(multilpa(long, c("y1", "y2"), "id", n_profiles = 2,
                          n_group_classes = 2, family = "full_cross_level",
                          n_starts = 2, seed = 1))
  expect_equal(cross$min_variance, smallest)
})

# K2-3 -----------------------------------------------------------------------

test_that("K2-3: person standardization of two-observation groups is refused or flagged", {
  withr::local_seed(9)
  make <- function(sizes) {
    g <- rep(seq_along(sizes), sizes)
    z <- stats::rbinom(length(g), 1, 0.5)
    x <- matrix(stats::rnorm(length(g) * 2), length(g), 2) + outer(z, c(1.5, 1.5))
    data.frame(g = g, a = x[, 1], b = x[, 2])
  }
  fit <- function(data, ...) {
    suppressMessages(multilpa(data, c("a", "b"), id = "g", n_profiles = 2,
                              n_group_classes = 1, centering = "person",
                              scaling = "person", n_starts = 2, seed = 1, ...))
  }
  # Every group of two: before, LL +1270.88 with every variance at the bound.
  expect_error(fit(make(rep(2, 150))), class = "latents_bad_data")
  # Some groups of two: before, no warning, and a profile pulled onto the
  # atoms (28.9% of its mass on 10.7% of the rows).
  mixed <- make(c(rep(5, 100), rep(2, 30)))
  expect_warning(fit(mixed), class = "latents_atomic_groups")
  message <- tryCatch(fit(mixed), latents_atomic_groups = conditionMessage)
  expect_match(message, "30 of 130 groups (10.7% of rows)", fixed = TRUE)
  # Groups of weight zero are not judged, and the other schemes are untouched.
  mixed$w <- ifelse(mixed$g > 100, 0, 1)
  expect_no_warning(suppressWarnings(fit(mixed, weights = "w"),
                                     classes = setdiff(c("latents_boundary",
                                                         "latents_thin_profiles"), "")))
  expect_no_condition(suppressMessages(multilpa(
    make(c(rep(5, 100), rep(2, 30))), c("a", "b"), id = "g", n_profiles = 2,
    n_group_classes = 1, centering = "person", scaling = "grand", n_starts = 2,
    seed = 1)), class = "latents_atomic_groups")
})
