# The Vuong-Lo-Mendell-Rubin test (Vuong 1989; Lo, Mendell and Rubin 2001)
# as Latent GOLD implements it (Vermunt 2024).
# The comparisons with independent references are in tests/equivalence/test-lmr.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

test_that("for regular nested models the weights are one per added parameter", {
  # Information equal to the variance of the score, and the smaller model's
  # scores the larger one's mapped by the Jacobian of the nesting: Vuong's
  # matrix then has p1 - p0 unit eigenvalues and zeros, so the reference
  # distribution is chi-square on p1 - p0 degrees of freedom.
  set.seed(8)
  scores <- matrix(stats::rnorm(400 * 5), 400, 5)
  jacobian <- matrix(stats::rnorm(5 * 3), 5, 3)
  information <- crossprod(scores)
  alternative <- list(scores = scores, information = information)
  null <- list(scores = scores %*% jacobian,
               information = t(jacobian) %*% information %*% jacobian)
  lambda <- sort(.latents_vuong_weights(alternative, null), decreasing = TRUE)
  expect_equal(lambda[1:2], c(1, 1), tolerance = 1e-8)
  # The rest are zero to the rounding of an unsymmetric eigendecomposition.
  expect_lt(max(abs(lambda[-(1:2)])), 1e-6)
})

test_that("the comparison table lists each neighbour pair, with its reason when unavailable", {
  skip_on_cran()
  data <- subset(course_engagement, student <= 60)
  grid <- quietly(enumerate_classes(data, c("browse", "lectures", "forum_read"),
                                    "student", n_profiles = 1:3,
                                    n_group_classes = 1:2, model = "VVI",
                                    n_starts = 2, seed = 1, tol = 1e-10))
  test <- get_results(grid, "lmr")
  expect_s3_class(test, "data.frame")
  expect_named(test, c("direction", "null_candidate", "alternative_candidate",
                       "n_profiles", "n_group_classes", "model", "n_units",
                       "parameters_added", "lr", "vlmr_p", "adjusted_lr",
                       "almr_p", "reference_mean", "reference_sd", "note"))
  expect_setequal(unique(test$direction), c("profiles", "group_classes"))
  # Every comparison adds one profile or one group class to its null model.
  table <- as.data.frame(grid)
  null_rows <- match(test$null_candidate, table$candidate)
  alternative_rows <- match(test$alternative_candidate, table$candidate)
  added <- (table$n_profiles[alternative_rows] - table$n_profiles[null_rows]) +
    (table$n_group_classes[alternative_rows] - table$n_group_classes[null_rows])
  expect_true(all(added == 1L))
  computed <- !is.na(test$lr)
  expect_true(all(test$vlmr_p[computed] >= 0 & test$vlmr_p[computed] <= 1))
  expect_true(all(test$adjusted_lr[computed] < test$lr[computed]))
  expect_true(all(test$n_units[computed] == 60L))
  # The table needs every candidate's information, so "all" leaves it out.
  expect_false("lmr" %in% names(get_results(grid, "all")))
  single <- quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                                  n_profiles = 2, model = "VVI", n_starts = 1,
                                  seed = 1))
  expect_error(get_results(single, "lmr"), class = "latents_no_lmr")
})

test_that("Imhof handles negative statistics and retains genuine small weights", {
  expect_equal(.latents_imhof_upper(-1, 1), 1, tolerance = 1e-10)
  expect_equal(.latents_imhof_upper(-1, -1), stats::pchisq(1, 1), tolerance = 1e-10)
  expect_equal(.latents_imhof_upper(-2, c(1, -1)),
               1 - .latents_imhof_upper(2, c(1, -1)), tolerance = 1e-10)
  invisible(lapply(c(1e9, 1e10, 1e11, 1e12), function(ratio) {
    expected <- stats::pf(ratio, 1, 1, lower.tail = FALSE)
    expect_lt(abs(.latents_imhof_upper(0, c(1, -ratio)) - expected), 1e-12)
    expect_lt(abs(.latents_imhof_upper(0, c(ratio, -1)) - (1 - expected)), 1e-12)
  }))
  expect_equal(.latents_imhof_upper(-1, c(0, 0)), 1)
  expect_equal(.latents_imhof_upper(0, c(0, 0)), 0)
})

test_that("Vuong weights are invariant to information coordinate units", {
  set.seed(8)
  scores <- matrix(stats::rnorm(400 * 5), 400, 5)
  jacobian <- matrix(stats::rnorm(5 * 3), 5, 3)
  alternative <- list(scores = scores, information = crossprod(scores))
  null_scores <- scores %*% jacobian
  null <- list(scores = null_scores, information = crossprod(null_scores))
  rescale <- function(parts, scale) {
    list(scores = sweep(parts$scores, 2L, scale, "*"),
         information = parts$information * tcrossprod(scale))
  }
  changed <- .latents_vuong_weights(rescale(alternative, c(1e-8, 1, 1e8, 1, 1)),
                                    rescale(null, c(1e8, 1e-8, 1)))
  expect_equal(sort(changed, decreasing = TRUE)[1:2], c(1, 1), tolerance = 1e-8)
  expect_lt(max(abs(sort(changed, decreasing = TRUE)[-(1:2)])), 1e-6)
})

test_that("public LMR retains its reference distribution when indicator units change", {
  set.seed(21)
  n <- 180
  class <- sample(1:2, n, TRUE)
  data <- data.frame(a = stats::rnorm(n, c(-1.5, 1.5)[class]),
                     b = stats::rnorm(n, c(-1.5, 1.5)[class]))
  fit_grid <- function(data, model, size = 1) {
    quietly(enumerate_lpa(data, c("a", "b"), n_profiles = 1:2,
                          model = model, n_starts = 3, seed = 1, tol = 1e-12,
                          min_variance = 1e-6 * size^2))
  }
  invisible(lapply(c("EEI", "EEE"), function(model) {
    reference <- get_results(fit_grid(data, model), "lmr")
    invisible(lapply(c(1e-8, 1e8), function(size) {
      result <- get_results(fit_grid(transform(data, a = a * size, b = b * size), model, size), "lmr")
      expect_true(is.na(result$note))
      expect_equal(result$lr, reference$lr, tolerance = 1e-7)
      expect_equal(result$reference_mean, reference$reference_mean, tolerance = 1e-6)
      expect_equal(result$reference_sd, reference$reference_sd, tolerance = 1e-6)
      expect_equal(result$vlmr_p, reference$vlmr_p, tolerance = 1e-7)
    }))
  }))
})

test_that("numerical null eigenvalues do not create tails in exact nested references", {
  set.seed(8)
  scores <- matrix(stats::rnorm(400 * 5), 400, 5)
  jacobian <- matrix(stats::rnorm(5 * 4), 5, 4)
  null_scores <- scores %*% jacobian
  alternative <- list(scores = scores, information = crossprod(scores))
  null <- list(scores = null_scores, information = crossprod(null_scores))
  lambda <- .latents_vuong_weights(alternative, null)
  expect_equal(sum(lambda != 0), 1)
  expect_equal(lambda[lambda != 0], 1, tolerance = 1e-12)
  invisible(lapply(c(0, 1e-12, 0.5, 5), function(q) {
    expect_equal(.latents_imhof_upper(q, lambda),
                 stats::pchisq(q, 1, lower.tail = FALSE), tolerance = 1e-10)
  }))
  # A real small eigenvalue is well above backward error and must survive.
  tiny <- list(scores = cbind(rep(1, 20), rep(c(-1, 1), 10) * 1e-5),
               information = diag(20, 2))
  no_parameters <- list(scores = matrix(0, 20, 1), information = matrix(1))
  expect_equal(sort(.latents_vuong_weights(tiny, no_parameters), decreasing = TRUE),
               c(1, 1e-10, 0), tolerance = 1e-15)
})

test_that("LMR reports a reversed likelihood as an unavailable comparison", {
  set.seed(21)
  n <- 100
  class <- sample(1:2, n, TRUE)
  data <- data.frame(a = stats::rnorm(n, c(-2, 2)[class]),
                     b = stats::rnorm(n, c(-2, 2)[class]))
  grid <- quietly(enumerate_lpa(data, c("a", "b"), n_profiles = 1:2,
                               model = "EEI", n_starts = 2, seed = 1, tol = 1e-12))
  # Inject a reversed likelihood at the scoring adapter, after the genuine
  # fits have passed data/likelihood validation. This exercises the public
  # table's invalid-comparison policy without forging observed-data scores.
  original <- .latents_vuong_parts
  null_likelihood <- get_candidate(grid, 1)$log_likelihood
  testthat::local_mocked_bindings(.latents_vuong_parts = function(fit) {
    parts <- original(fit)
    if (fit$n_profiles == 2L) parts$log_likelihood <- null_likelihood - 1
    parts
  })
  result <- get_results(grid, "lmr")
  expect_equal(nrow(result), 1L)
  expect_true(is.na(result$vlmr_p))
  expect_match(result$note, "lower likelihood")
})

test_that("LMR uses the conditioned chart for translated direct-effect predictors", {
  set.seed(77)
  n <- 240
  class <- sample(1:2, n, TRUE)
  z <- stats::rnorm(n, 1, 2)
  data <- data.frame(z = z, a = stats::rnorm(n, c(-2, 2)[class]) + 0.7 * z,
                     b = stats::rnorm(n, c(-2, 2)[class]))
  fit <- function(data, size = 1) {
    quietly(enumerate_lpa(data, c("a", "b"), n_profiles = 1:2, model = "EEI",
                          direct_effects = "a ~ z", n_starts = 3, seed = 1,
                          tol = 1e-12, min_variance = 1e-6 * size^2))
  }
  base <- get_results(fit(data), "lmr")
  shifted <- get_results(fit(transform(data, z = z + 1e6)), "lmr")
  scaled <- get_results(fit(transform(data, a = a * 1e8, b = b * 1e8,
                                      z = z + 1e6), 1e8), "lmr")
  invisible(lapply(list(shifted, scaled), function(result) {
    expect_true(is.na(result$note))
    expect_equal(result$lr, base$lr, tolerance = 1e-7)
    expect_equal(result$reference_mean, base$reference_mean, tolerance = 1e-6)
    expect_equal(result$reference_sd, base$reference_sd, tolerance = 1e-6)
  }))
})

test_that("LMR adjustment counts informative positive-weight independent units", {
  set.seed(21)
  n <- 180
  class <- sample(1:2, n, TRUE)
  data <- data.frame(a = stats::rnorm(n, c(-1.5, 1.5)[class]),
                     b = stats::rnorm(n, c(-1.5, 1.5)[class]), weight = 1)
  empty <- data.frame(a = rep(NA_real_, 30), b = rep(NA_real_, 30), weight = 1)
  fit <- function(data, ...) {
    quietly(enumerate_lpa(data, c("a", "b"), n_profiles = 1:2, model = "EEI",
                          missing = "fiml", n_starts = 3, seed = 1, tol = 1e-12, ...))
  }
  base <- get_results(fit(data), "lmr")
  padded <- get_results(fit(rbind(data, empty)), "lmr")
  columns <- c("n_units", "lr", "adjusted_lr", "reference_mean", "reference_sd",
               "vlmr_p", "almr_p")
  expect_equal(padded[columns], base[columns], tolerance = 1e-6)
  legacy <- get_results(fit(rbind(data, empty), drop_empty = FALSE), "lmr")
  expect_equal(legacy$n_units, n + 30L)
  expect_equal(legacy$adjusted_lr,
               legacy$lr / (1 + 1 / (legacy$parameters_added * log(n + 30L))))
  # Zero-weight units are excluded even when their indicators are observed.
  zero <- transform(data[1:30, ], weight = 0)
  weighted <- get_results(fit(rbind(data, zero), weights = "weight"), "lmr")
  expect_equal(weighted$n_units, n)
  expect_equal(weighted$adjusted_lr / weighted$lr,
               base$adjusted_lr / base$lr, tolerance = 1e-10)
})

test_that("Imhof remains accurate across signs at nearly zero statistics", {
  invisible(lapply(c(-1e-12, 1e-12), function(q) {
    expect_equal(.latents_imhof_upper(q, 2),
                 stats::pchisq(q / 2, 1, lower.tail = FALSE), tolerance = 1e-10)
    expect_equal(.latents_imhof_upper(q, -2),
                 if (q < 0) stats::pchisq(-q / 2, 1) else 0, tolerance = 1e-10)
    expect_equal(.latents_imhof_upper(q, c(1, -1)),
                 0.5, tolerance = 1e-10)
  }))
  cases <- expand.grid(ratio = c(1e-12, 1e12, 1e24), m = c(1, 3), k = c(1, 4))
  invisible(lapply(seq_len(nrow(cases)), function(index) {
    case <- cases[index, ]
    actual <- with(case, .latents_imhof_upper(0, c(rep(1, m), rep(-ratio, k))))
    expected <- with(case, stats::pf(ratio * k / m, m, k, lower.tail = FALSE))
    expect_lt(abs(actual - expected), 1e-10)
  }))
})

test_that("Imhof extreme weight ranges avoid overflow or fail with a typed reason", {
  weights <- c(rep(1e-100, 20), -1e100)
  expected <- stats::pf(1e200 / 20, 20, 1, lower.tail = FALSE)
  expect_lt(abs(.latents_imhof_upper(0, weights) - expected), 1e-10)
  expect_error(.latents_imhof_upper(0, c(1e308, -1e-308)),
                class = "latents_inaccurate_probability")
})

test_that("Vuong numerical zeros account for conditioned information whitening", {
  # A near-collinear chart preserves the exact four-dimensional null span.
  # Its cross-product information amplifies Cholesky roundoff; this fixture
  # reproduces locally the platform-dependent extra weight seen on R 4.1.
  set.seed(83)
  scores <- matrix(stats::rnorm(400 * 5), 400, 5)
  alternative <- list(scores = scores, information = crossprod(scores))
  invisible(lapply(c(1e-2, 1e-3, 1e-4), function(separation) {
    null_scores <- cbind(scores[, 1], scores[, 1] + separation * scores[, 2],
                          scores[, 3:4])
    null <- list(scores = null_scores, information = crossprod(null_scores))
    lambda <- tryCatch(.latents_vuong_weights(alternative, null),
                       latents_inaccurate_probability = identity)
    if (inherits(lambda, "condition")) {
      expect_s3_class(lambda, "latents_inaccurate_probability")
      return(invisible(NULL))
    }
    expect_equal(sum(lambda != 0), 1)
    expect_equal(lambda[lambda != 0], 1, tolerance = 1e-12)
    expect_equal(.latents_imhof_upper(0, lambda), 1, tolerance = 1e-12)
    expect_equal(.latents_imhof_upper(1e-12, lambda),
                 stats::pchisq(1e-12, 1, lower.tail = FALSE), tolerance = 1e-10)
  }))
})

test_that("uncertain genuine Vuong weights are refused rather than replaced by zeros", {
  set.seed(83)
  scores <- qr.Q(qr(matrix(stats::rnorm(400 * 5), 400, 5)))
  alternative <- list(scores = scores, information = diag(5))
  invisible(lapply(list(c(1e-4, 1e-8), c(1e-3, 1e-10)), function(case) {
    transform <- diag(4)
    transform[1, 2] <- 1
    transform[2, 2] <- case[1]
    delta <- case[2]
    null <- list(scores = scores[, 1:4] %*% transform * sqrt(1 + delta),
                 information = crossprod(transform))
    # Exactly one +1 and four -delta weights: the upper tail at zero is
    # pf(4 delta, 1, 4), below one by 1.5e-4 or 1.5e-5. Uncertain signs
    # cannot justify silently reporting probability one.
    expect_error(.latents_vuong_weights(alternative, null),
                  class = "latents_inaccurate_probability")
    plain <- list(scores = scores[, 1:4] * sqrt(1 + delta), information = diag(4))
    lambda <- .latents_vuong_weights(alternative, plain)
    expect_equal(sort(lambda), c(rep(-delta, 4), rep(0, 4), 1), tolerance = 1e-12)
    expect_equal(.latents_imhof_upper(0, lambda),
                 stats::pf(4 * delta, 1, 4, lower.tail = FALSE), tolerance = 1e-8)
  }))
})

test_that("nonfinite or overwhelming Vuong uncertainty cannot erase finite weights", {
  set.seed(83)
  scores <- qr.Q(qr(matrix(stats::rnorm(80), 20, 4)))
  transform <- diag(4)
  transform[1, 2] <- 1
  transform[2, 2] <- 1e-7
  null <- list(scores = matrix(0, 20, 1), information = matrix(1))
  invisible(lapply(c(8e153, 1.1e154), function(size) {
    alternative <- list(scores = scores %*% transform * size,
                         information = crossprod(transform))
    expect_error(.latents_vuong_weights(alternative, null),
                  class = "latents_inaccurate_probability")
  }))
})
