# Properties of the Vuong-Lo-Mendell-Rubin test (Vuong 1989; Lo, Mendell and
# Rubin 2001): the Imhof tail, Vuong's weights and the public table. The
# comparisons with Latent GOLD, Mplus and an independent reference are in
# equivalence/test-features-vlmr.R.

# Base R draws: simulab is not among the package's Suggests.
vlmr_small_data <- function(n = 150, seed = 9) {
  set.seed(seed)
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  data.frame(a = stats::rnorm(n, c(0, 1.6)[class]), b = stats::rnorm(n, c(0, 1.4)[class]),
             c = stats::rnorm(n, c(0, 1.2)[class]))
}

vlmr_small_grid <- function(data, size = 1, model = "EEI") {
  quietly(enumerate_lpa(data, c("a", "b", "c"), n_profiles = 1:3, model = model,
                        n_starts = 5, seed = 1, tol = 1e-12,
                        min_variance = 1e-6 * size^2))
}

vlmr_columns <- c("n_units", "parameters_added", "lr", "vlmr_p", "adjusted_lr", "almr_p",
                  "reference_mean", "reference_sd")

test_that("Imhof probabilities lie in [0, 1] and fall as the statistic grows", {
  set.seed(31)
  weights <- lapply(seq_len(40), function(index) {
    size <- sample(1:6, 1L)
    stats::rnorm(size) * 10^stats::runif(size, -2, 2)
  })
  invisible(lapply(weights, function(lambda) {
    q <- c(0, sort(stats::rexp(6, 1 / (sum(abs(lambda)) + 1))))
    p <- vapply(q, .latents_imhof_upper, numeric(1), lambda = lambda)
    expect_true(all(p >= 0 & p <= 1))
    # P(Q > q) is non-increasing in q.
    expect_true(all(diff(p) <= 1e-12))
  }))
})

test_that("Imhof reproduces exact chi-square, F and scaled probabilities", {
  set.seed(32)
  invisible(lapply(seq_len(15), function(index) {
    scale <- 10^stats::runif(1L, -6, 6)
    df <- sample(1:12, 1L)
    q <- scale * stats::rexp(1L, 1 / df)
    expect_lt(abs(.latents_imhof_upper(q, rep(scale, df)) -
                    stats::pchisq(q / scale, df, lower.tail = FALSE)), 1e-9)
    # P(a X_m - b Y_k > 0) = P(F(m, k) > b k / (a m)).
    m <- sample(1:5, 1L)
    k <- sample(1:5, 1L)
    ratio <- 10^stats::runif(1L, -3, 3)
    expect_lt(abs(.latents_imhof_upper(0, scale * c(rep(1, m), rep(-ratio, k))) -
                    stats::pf(ratio * k / m, m, k, lower.tail = FALSE)), 1e-9)
  }))
})

test_that("Imhof is invariant to joint scaling of the weights and the statistic", {
  set.seed(33)
  invisible(lapply(seq_len(10), function(index) {
    lambda <- stats::rnorm(5) * 10^stats::runif(5, -1, 1)
    q <- abs(stats::rnorm(1L, sum(lambda), 2))
    base <- .latents_imhof_upper(q, lambda)
    invisible(lapply(c(1e-6, 1e6), function(size) {
      expect_lt(abs(.latents_imhof_upper(q * size, lambda * size) - base), 1e-10)
    }))
  }))
})

test_that("Wynn's epsilon algorithm accelerates an alternating series to its limit", {
  partial <- cumsum((-1)^(0:14) / seq_len(15))
  expect_lt(abs(.latents_wynn(partial) - log(2)), 1e-10)
  expect_gt(abs(partial[15] - log(2)), 1e-2)
  expect_identical(.latents_wynn(c(1, 1, 1)), 1)
})

test_that("regular nested scores give one unit weight per added parameter", {
  # Information equal to the score cross-product and null scores that are the
  # alternative's mapped by a Jacobian: chi-square on p1 - p0 degrees.
  invisible(lapply(1:6, function(seed) {
    set.seed(seed)
    p1 <- sample(4:8, 1L)
    p0 <- sample(1:(p1 - 1L), 1L)
    scores <- matrix(stats::rnorm(300 * p1), 300, p1)
    jacobian <- matrix(stats::rnorm(p1 * p0), p1, p0)
    null_scores <- scores %*% jacobian
    lambda <- .latents_vuong_weights(list(scores = scores, information = crossprod(scores)),
                                     list(scores = null_scores,
                                          information = crossprod(null_scores)))
    expect_length(lambda, p1 + p0)
    expect_equal(sort(lambda[lambda != 0]), rep(1, p1 - p0), tolerance = 1e-10)
  }))
})

test_that("Vuong weights do not depend on the order or the units of the coordinates", {
  set.seed(34)
  scores <- matrix(stats::rnorm(200 * 5), 200, 5)
  alternative <- list(scores = scores, information = crossprod(scores) * 1.3 + diag(5))
  null_scores <- scores[, 1:3] + matrix(stats::rnorm(600, sd = 0.3), 200, 3)
  null <- list(scores = null_scores, information = crossprod(null_scores) + diag(3))
  base <- sort(.latents_vuong_weights(alternative, null))
  # Relabelling classes permutes a model's parameters.
  permute <- function(parts, order) {
    list(scores = parts$scores[, order], information = parts$information[order, order])
  }
  expect_equal(sort(.latents_vuong_weights(permute(alternative, c(4, 2, 5, 1, 3)),
                                           permute(null, c(3, 1, 2)))), base,
               tolerance = 1e-10)
  # A change of units is a congruence of the information.
  rescale <- function(parts, scale) {
    list(scores = sweep(parts$scores, 2L, scale, "*"),
         information = parts$information * tcrossprod(scale))
  }
  expect_equal(sort(.latents_vuong_weights(rescale(alternative, 10^c(-6, 3, 0, 6, -2)),
                                           rescale(null, 10^c(4, -4, 1)))), base,
               tolerance = 1e-8)
  # Exchanging the two models negates every weight.
  expect_equal(sort(-.latents_vuong_weights(null, alternative)), base, tolerance = 1e-10)
})

test_that("singular or nonpositive information is refused with its class", {
  set.seed(35)
  scores <- matrix(stats::rnorm(100 * 3), 100, 3)
  good <- list(scores = scores[, 1:2], information = crossprod(scores[, 1:2]))
  collinear <- cbind(scores, scores[, 1] + scores[, 2])
  expect_error(.latents_vuong_weights(list(scores = collinear,
                                           information = crossprod(collinear)), good),
               class = "latents_inaccurate_probability")
  information <- crossprod(scores)
  information[2, 2] <- 0
  expect_error(.latents_vuong_weights(list(scores = scores, information = information), good),
               class = "latents_singular_information")
  indefinite <- crossprod(scores)
  indefinite[1, 2] <- indefinite[2, 1] <- 2 * sqrt(indefinite[1, 1] * indefinite[2, 2])
  expect_error(.latents_vuong_weights(list(scores = scores, information = indefinite), good),
               class = "latents_singular_information")
  expect_error(.latents_imhof_upper(0, c(1e308, -1e-308)),
               class = "latents_inaccurate_probability")
  lower <- list(scores = scores, information = crossprod(scores), log_likelihood = -10,
                n_parameters = 3, n_units = 100)
  higher <- list(scores = scores[, 1:2], information = crossprod(scores[, 1:2]),
                 log_likelihood = -9, n_parameters = 2, n_units = 100)
  expect_error(.latents_vuong_row(lower, higher), class = "latents_invalid_lmr")
})

test_that("the comparison row applies the adjustment and the moments of its weights", {
  set.seed(36)
  scores <- matrix(stats::rnorm(250 * 4), 250, 4)
  alternative <- list(scores = scores, information = crossprod(scores) * 1.2,
                      log_likelihood = -100, n_parameters = 4, n_units = 250)
  null <- list(scores = scores[, 1:2], information = crossprod(scores[, 1:2]),
               log_likelihood = -104.5, n_parameters = 2, n_units = 250)
  row <- .latents_vuong_row(alternative, null)
  lambda <- .latents_vuong_weights(alternative, null)
  expect_named(row, vlmr_columns)
  expect_equal(row$lr, 9)
  expect_equal(row$adjusted_lr, 9 / (1 + 1 / (2 * log(250))), tolerance = 1e-14)
  expect_equal(row$reference_mean, sum(lambda), tolerance = 1e-14)
  expect_equal(row$reference_sd, sqrt(2 * sum(lambda^2)), tolerance = 1e-14)
  expect_equal(row$vlmr_p, .latents_imhof_upper(9, lambda), tolerance = 1e-14)
  expect_gt(row$almr_p, row$vlmr_p)
})

test_that("the public table is invariant to row order, units and origin", {
  data <- vlmr_small_data()
  base <- get_results(vlmr_small_grid(data), "lmr")
  expect_named(base, c("direction", "null_candidate", "alternative_candidate", "n_profiles",
                       "n_group_classes", "model", vlmr_columns, "note"))
  expect_identical(nrow(base), 2L)
  expect_true(all(is.na(base$note)))
  expect_true(all(base$vlmr_p >= 0 & base$vlmr_p <= 1 & base$almr_p >= base$vlmr_p))
  expect_equal(base$adjusted_lr,
               base$lr / (1 + 1 / (base$parameters_added * log(base$n_units))),
               tolerance = 1e-14)
  # EM stops on the likelihood (tol 1e-12), so a different path through the
  # same maximum moves the estimates, and with them the weights, by ~1e-6.
  set.seed(37)
  shuffled <- get_results(vlmr_small_grid(data[sample(nrow(data)), ]), "lmr")
  expect_equal(shuffled[vlmr_columns], base[vlmr_columns], tolerance = 1e-5)
  moved <- get_results(vlmr_small_grid(data * 1e-4 + 1e6, size = 1e-4), "lmr")
  expect_true(all(is.na(moved$note)))
  expect_equal(moved[vlmr_columns], base[vlmr_columns], tolerance = 1e-5)
})

test_that("unavailable comparisons keep their row with the reason", {
  # Two identical items: probabilities on their bounds, and a rank-deficient
  # one-class score.
  set.seed(5)
  data <- data.frame(u1 = rep(1:2, each = 20), u2 = rep(1:2, each = 20),
                     u3 = sample(1:2, 40, TRUE))
  # Three classes on three binary items are not identified (Goodman 1974).
  grid <- quietly(enumerate_lca(data, c("u1", "u2", "u3"), n_classes = 1:3, n_starts = 3,
                                seed = 1),
                  classes = c(.multilpa_expected_warnings, "latents_unidentified"))
  table <- get_results(grid, "lmr")
  expect_identical(nrow(table), 2L)
  expect_true(all(is.na(as.matrix(table[vlmr_columns]))))
  expect_true(all(nzchar(table$note)))
  expect_match(table$note[2L], "bound")
  single <- quietly(enumerate_lpa(vlmr_small_data(), c("a", "b"), n_profiles = 2,
                                  model = "EEI", n_starts = 1, seed = 1))
  expect_error(get_results(single, "lmr"), class = "latents_no_lmr")
})

test_that("a one-class null true model gives a large VLMR p-value on average", {
  skip_on_cran()
  # Under the null the p-value is roughly uniform; across 12 datasets from
  # one class its mean stays well inside (0.2, 0.8).
  p <- vapply(1:12, function(seed) {
    set.seed(seed)
    data <- data.frame(a = stats::rnorm(120), b = stats::rnorm(120), c = stats::rnorm(120))
    grid <- quietly(enumerate_lpa(data, c("a", "b", "c"), n_profiles = 1:2, model = "EEI",
                                  n_starts = 5, seed = 1, tol = 1e-10))
    get_results(grid, "lmr")$vlmr_p
  }, numeric(1))
  expect_true(all(p >= 0 & p <= 1))
  expect_gt(mean(p), 0.2)
  expect_lt(mean(p), 0.8)
})
