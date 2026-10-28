# Properties of the ML step-three distal outcome model (Vermunt 2010; Bakk,
# Tekle and Vermunt 2013) that hold whatever the reference: equivariance
# under the outcome's units, invariance to row order and class labels, the
# limits of perfect classification and equal variances, and the refusals.
# The comparisons with Latent GOLD, Mplus and an independent implementation
# are in equivalence/test-features-ml-three-step.R.

features_ml_data <- function(seed = 3L, n = 300L, separation = 2.4, k = 2L) {
  set.seed(seed)
  truth <- sample.int(k, n, replace = TRUE, prob = c(0.45, 0.35, 0.2)[seq_len(k)])
  centers <- c(-1, 1, 0)[seq_len(k)] * separation / 2
  offsets <- c(0, 0, 1.6)[seq_len(k)] * separation / 2
  data.frame(id = seq_len(n),
             a = stats::rnorm(n, centers[truth]),
             b = stats::rnorm(n, centers[truth] + offsets[truth]),
             y = stats::rnorm(n, c(10, 14, 11)[truth], c(2, 3, 1.5)[truth]),
             truth = truth)
}

features_ml_fit <- function(data, k = 2L) {
  quietly(lpa(data, c("a", "b"), n_profiles = k, n_starts = 4, seed = 1, tol = 1e-12,
              max_iter = 10000))
}

features_ml_table <- function(...) as.data.frame(three_step(..., method = "ml"))

test_that("an affine map of the outcome maps the estimates and keeps the tests", {
  data <- features_ml_data()
  fit <- features_ml_fit(data)
  mapped <- transform(data, y = -37 - 2.5 * y)
  invisible(lapply(c("varying", "equal"), \(variance) {
    base <- features_ml_table(fit, data, "y", variance = variance)
    moved <- features_ml_table(fit, mapped, "y", variance = variance)
    expect_equal(moved$estimate, -37 - 2.5 * base$estimate, tolerance = 1e-8)
    expect_equal(moved$standard_error, 2.5 * base$standard_error, tolerance = 1e-6)
    expect_equal(moved$variance, 6.25 * base$variance, tolerance = 1e-8)
    expect_equal(moved$effective_n, base$effective_n, tolerance = 1e-8)
    # A negative slope swaps the interval's ends.
    expect_equal(moved$conf_low, -37 - 2.5 * base$conf_high, tolerance = 1e-6)
    base_pairs <- features_ml_table(fit, data, "y", variance = variance, contrast = "pairs")
    moved_pairs <- features_ml_table(fit, mapped, "y", variance = variance, contrast = "pairs")
    expect_equal(moved_pairs$statistic, -base_pairs$statistic, tolerance = 1e-6)
    expect_equal(moved_pairs$p_value, base_pairs$p_value, tolerance = 1e-6)
  }))
})

test_that("the order of the rows does not change the step-three model", {
  data <- features_ml_data(k = 3L, n = 360L)
  fit <- features_ml_fit(data, k = 3L)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  set.seed(8)
  order <- sample.int(nrow(data))
  shuffled <- list(posteriors = pieces$posteriors[order, , drop = FALSE],
                   modal = pieces$modal[order], n_classes = pieces$n_classes,
                   group_index = seq_along(order))
  invisible(lapply(c("varying", "equal"), \(variance) {
    run <- function(p, values, contrast = "none") {
      as.data.frame(.multilpa_three_step_ml(p, values, "individuals", variance, 0.95,
                                            contrast, "BH", "independent"))
    }
    expect_equal(run(shuffled, data$y[order]), run(pieces, data$y), tolerance = 1e-7)
    expect_equal(run(shuffled, data$y[order], "pairs"), run(pieces, data$y, "pairs"),
                 tolerance = 1e-7)
  }))
  # Refitting the shuffled data reaches the same classes, possibly renumbered.
  refit <- features_ml_fit(data[order, ], k = 3L)
  sorted <- function(x) x[order(x$estimate), c("estimate", "standard_error", "variance")]
  expect_equal(sorted(features_ml_table(refit, data[order, ], "y")),
               sorted(features_ml_table(fit, data, "y")), tolerance = 1e-5,
               ignore_attr = TRUE)
})

test_that("renumbering the classes of the measurement model permutes the result", {
  data <- features_ml_data(k = 3L, n = 360L)
  fit <- features_ml_fit(data, k = 3L)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  relabel <- c(3L, 1L, 2L)
  renumbered <- list(posteriors = pieces$posteriors[, relabel], modal = match(pieces$modal, relabel),
                     n_classes = 3L, group_index = pieces$group_index)
  run <- function(p, contrast = "none") {
    as.data.frame(.multilpa_three_step_ml(p, data$y, "individuals", "varying", 0.95,
                                          contrast, "BH", "independent"))
  }
  base <- run(pieces)
  moved <- run(renumbered)
  columns <- c("estimate", "standard_error", "variance", "effective_n")
  expect_equal(moved[columns], base[relabel, columns], tolerance = 1e-7, ignore_attr = TRUE)
  # Pair (a, b) of the renumbered classes is pair (relabel[a], relabel[b]).
  base_pairs <- run(pieces, "pairs")
  moved_pairs <- run(renumbered, "pairs")
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  matched <- match(key(relabel[moved_pairs$class], relabel[moved_pairs$reference_class]),
                   key(base_pairs$class, base_pairs$reference_class))
  sign <- ifelse(relabel[moved_pairs$class] == base_pairs$class[matched], 1, -1)
  expect_equal(moved_pairs$estimate, sign * base_pairs$estimate[matched], tolerance = 1e-7)
  expect_equal(moved_pairs$standard_error, base_pairs$standard_error[matched], tolerance = 1e-7)
})

test_that("equal variances give one shared variance and are recorded", {
  data <- features_ml_data(k = 3L, n = 360L)
  fit <- features_ml_fit(data, k = 3L)
  equal <- three_step(fit, data, "y", method = "ml", variance = "equal")
  expect_length(unique(equal$variance), 1L)
  expect_identical(attr(equal, "variance"), "equal")
  varying <- three_step(fit, data, "y", method = "ml")
  expect_identical(attr(varying, "variance"), "varying")
  expect_gt(length(unique(round(varying$variance, 8))), 1L)
  # The pooled variance is the responsibility-weighted average of the class
  # variances only at the same means; it lies within their range here.
  expect_true(equal$variance[1L] > min(varying$variance) &&
                equal$variance[1L] < max(varying$variance))
})

test_that("with perfect classification ML gives the modal means and their sandwich", {
  data <- features_ml_data(separation = 40)
  fit <- features_ml_fit(data)
  ml <- features_ml_table(fit, data, "y", vcov_type = "independent")
  modal <- as.data.frame(three_step(fit, data, "y", method = "modal",
                                    vcov_type = "independent"))
  expect_equal(ml$estimate, modal$estimate, tolerance = 1e-10)
  groups <- split(data$y, fit$subject_profiles)
  n <- nrow(data)
  # The sandwich of a normal mean is sqrt(sum of squared deviations) / n_c,
  # scaled by the CR1 factor n / (n - 1).
  closed <- vapply(groups, \(v) sqrt(sum((v - mean(v))^2)) / length(v), numeric(1)) *
    sqrt(n / (n - 1))
  expect_equal(ml$standard_error, unname(closed), tolerance = 1e-7)
  expect_equal(ml$variance, unname(vapply(groups, \(v) mean((v - mean(v))^2), numeric(1))),
               tolerance = 1e-10)
  equal <- features_ml_table(fit, data, "y", variance = "equal")
  pooled <- sum(vapply(groups, \(v) sum((v - mean(v))^2), numeric(1))) / n
  expect_equal(equal$variance, rep(pooled, 2L), tolerance = 1e-10)
})

test_that("intervals and pair tests follow from the estimates and errors", {
  data <- features_ml_data(k = 3L, n = 360L)
  fit <- features_ml_fit(data, k = 3L)
  classes <- features_ml_table(fit, data, "y", ci_level = 0.9)
  expect_equal(classes$conf_high - classes$estimate,
               stats::qnorm(0.95) * classes$standard_error, tolerance = 1e-12)
  pairs <- features_ml_table(fit, data, "y", contrast = "pairs", adjust = "holm")
  expect_identical(nrow(pairs), 3L)
  expect_equal(pairs$estimate, classes$estimate[pairs$class] -
                 classes$estimate[pairs$reference_class], tolerance = 1e-10)
  expect_equal(pairs$statistic, pairs$estimate / pairs$standard_error, tolerance = 1e-12)
  expect_equal(pairs$p_value, 2 * stats::pnorm(-abs(pairs$statistic)), tolerance = 1e-12)
  expect_equal(pairs$p_adjusted, stats::p.adjust(pairs$p_value, "holm"), tolerance = 1e-12)
  expect_identical(attr(pairs, "adjust"), "holm")
})

test_that("the step-three estimate is a maximum with a vanishing score", {
  data <- features_ml_data(k = 3L, n = 360L)
  fit <- features_ml_fit(data, k = 3L)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  invisible(lapply(c("varying", "equal"), \(variance) {
    step <- .multilpa_step_ml_fit(pieces, data$y, variance)
    expect_true(step$converged)
    at <- .multilpa_step_ml_coordinates(step, variance)
    best <- .multilpa_step_ml_log_likelihood(step, data$y, variance, at)
    expect_equal(best, step$log_likelihood, tolerance = 1e-10)
    nudged <- vapply(seq_along(at), \(j) {
      vapply(c(-1e-3, 1e-3), \(h) {
        moved <- at
        moved[j] <- moved[j] + h
        .multilpa_step_ml_log_likelihood(step, data$y, variance, moved)
      }, numeric(1))
    }, numeric(2))
    expect_true(all(nudged < best))
    expect_lt(max(abs(colSums(.multilpa_step_ml_scores(step, data$y, variance)))), 1e-6)
  }))
})

test_that("a single-level fit's clustered and independent errors coincide", {
  data <- features_ml_data()
  fit <- features_ml_fit(data)
  expect_identical(features_ml_table(fit, data, "y")$standard_error,
                   features_ml_table(fit, data, "y", vcov_type = "independent")$standard_error)
})

test_that("weighted fits and impossible outcomes are refused by class", {
  data <- features_ml_data()
  fit <- features_ml_fit(data)
  weighted <- quietly(lpa(transform(data, w = rep(1:2, length.out = nrow(data))),
                          c("a", "b"), n_profiles = 2, n_starts = 2, seed = 1, weights = "w"))
  expect_error(three_step(weighted, data, "y", method = "ml"),
               class = "latents_unsupported_weights")
  expect_error(three_step(fit, transform(data, y = 5), "y", method = "ml"),
               class = "latents_bad_outcome")
  expect_error(three_step(fit, data, "a", method = "ml"), class = "multilpa_indicator_reused")
  expect_error(three_step(fit, data, "y", method = "modal", variance = "equal"),
               class = "latents_bad_argument")
  expect_error(three_step(fit, data, "y", method = "proportional", variance = "varying"),
               class = "latents_bad_argument")
  expect_error(three_step(fit, data, "y", method = "ml", variance = "pooled"),
               "should be one of")
  # Classed since 2026-10-27: a missing column names it, and NA is refused.
  expect_error(three_step(fit, data, "absent", method = "ml"), "`absent`",
               class = "latents_bad_data")
  with_missing <- data
  with_missing$y[4L] <- NA
  expect_error(three_step(fit, with_missing, "y", method = "ml"), "must not be missing",
               class = "latents_bad_argument")
})

test_that("standard errors keep their relative size at large origins and small units", {
  skip_on_cran()
  data <- features_ml_data(k = 3L, n = 360L)
  fit <- features_ml_fit(data, k = 3L)
  invisible(lapply(c("varying", "equal"), \(variance) {
    base <- features_ml_table(fit, data, "y", variance = variance)
    offset <- features_ml_table(fit, transform(data, y = y + 1e6), "y", variance = variance)
    small <- features_ml_table(fit, transform(data, y = y * 1e-4), "y", variance = variance)
    expect_equal(offset$standard_error, base$standard_error, tolerance = 1e-6)
    expect_equal(offset$estimate - 1e6, base$estimate, tolerance = 1e-7)
    expect_equal(small$standard_error * 1e4, base$standard_error, tolerance = 1e-6)
    expect_equal(small$estimate * 1e4, base$estimate, tolerance = 1e-8)
    pairs <- features_ml_table(fit, data, "y", variance = variance, contrast = "pairs")
    small_pairs <- features_ml_table(fit, transform(data, y = y * 1e-4), "y",
                                     variance = variance, contrast = "pairs")
    expect_equal(small_pairs$statistic, pairs$statistic, tolerance = 1e-6)
  }))
})

test_that("the step-three likelihood is the fixed-error mixture of Vermunt (2010)", {
  # Written out here from the formula: pi_t the average posterior, D[t, s]
  # the posterior mass of class t among units assigned to s over the mass of
  # t, and the outcome normal within each true class.
  data <- features_ml_data(k = 3L, n = 360L, separation = 1.6)
  fit <- features_ml_fit(data, k = 3L)
  posteriors <- fit$subject_posteriors
  modal <- max.col(posteriors, ties.method = "first")
  expect_identical(as.integer(fit$subject_profiles), modal)
  errors <- vapply(1:3, \(s) colSums(posteriors * (modal == s)), numeric(3)) /
    colSums(posteriors)
  prior <- colMeans(posteriors)
  log_likelihood <- function(means, sds) {
    terms <- vapply(1:3, \(t) prior[t] * errors[t, modal] *
                      stats::dnorm(data$y, means[t], sds[t]), numeric(nrow(data)))
    sum(log(rowSums(terms)))
  }
  invisible(lapply(c("varying", "equal"), \(variance) {
    step <- .multilpa_step_ml_fit(.multilpa_level_assignments(fit, "individuals"), data$y,
                                  variance)
    sds <- sqrt(step$variances)
    expect_equal(log_likelihood(step$means, sds), step$log_likelihood, tolerance = 1e-10)
    # The estimate is stationary for this likelihood: a step of 1e-4 standard
    # deviations in any mean changes it only at second order.
    slopes <- vapply(1:3, \(t) {
      h <- 1e-4 * sds[t]
      up <- step$means
      up[t] <- up[t] + h
      down <- step$means
      down[t] <- down[t] - h
      (log_likelihood(up, sds) - log_likelihood(down, sds)) / (2 * h) * sds[t]
    }, numeric(1))
    expect_lt(max(abs(slopes)), 1e-4)
  }))
})

# Posteriors and a heavy-tailed outcome on which the step-three likelihood
# has several maxima: from the proportional start alone EM stopped 2.57 below
# the best one (equivalence/features/ml-three-step/repro-local-maximum.R).
features_ml_heavy <- function(seed = 5L, n = 300L) {
  set.seed(seed)
  class <- sample.int(3L, n, TRUE)
  logits <- matrix(stats::rnorm(n * 3L, 0, 0.4), n, 3L)
  logits[cbind(seq_len(n), class)] <- logits[cbind(seq_len(n), class)] + 1.3
  posteriors <- exp(logits) / rowSums(exp(logits))
  list(pieces = list(posteriors = posteriors, modal = max.col(posteriors, "first"),
                     n_classes = 3L),
       y = c(0, 0.6, 1.2)[class] + stats::rt(n, 3))
}

test_that("the step-three model reaches the best of many random starts", {
  skip_on_cran()
  case <- features_ml_heavy()
  pieces <- case$pieces
  errors <- vapply(1:3, \(s) colSums(pieces$posteriors * (pieces$modal == s)), numeric(3)) /
    colSums(pieces$posteriors)
  prior <- colMeans(pieces$posteriors)
  log_likelihood <- function(theta) {
    terms <- vapply(1:3, \(t) prior[t] * errors[t, pieces$modal] *
                      stats::dnorm(case$y, theta[t], exp(theta[3L + t])), numeric(length(case$y)))
    sum(log(rowSums(terms)))
  }
  set.seed(11)
  best <- max(vapply(seq_len(20L), \(i) {
    start <- c(stats::quantile(case$y, stats::runif(3L), names = FALSE), rep(log(stats::sd(case$y)), 3L))
    stats::optim(start, log_likelihood, method = "BFGS",
                 control = list(fnscale = -1, maxit = 2000, reltol = 1e-14))$value
  }, numeric(1)))
  step <- .multilpa_step_ml_fit(pieces, case$y, "varying")
  expect_gt(step$log_likelihood, -880.5)
  expect_gt(step$log_likelihood, best - 1e-6 * (1 + abs(best)))
})

test_that("a class carried by one unit gets no standard error, with a warning", {
  set.seed(4)
  n <- 60L
  posteriors <- matrix(c(0.9, 0.1, 0), n, 3L, byrow = TRUE)
  posteriors[seq_len(n / 2), ] <- matrix(c(0.1, 0.9, 0), n / 2, 3L, byrow = TRUE)
  posteriors[n, ] <- c(0, 0, 1)
  pieces <- list(posteriors = posteriors, modal = max.col(posteriors, "first"), n_classes = 3L)
  y <- stats::rnorm(n, c(0, 1, 3)[pieces$modal])
  expect_warning(result <- .multilpa_three_step_ml(pieces, y, "individuals", "equal", 0.95,
                                                   "none", "BH", "independent"),
                 class = "latents_small_classes")
  expect_true(is.na(result$standard_error[3L]) && is.na(result$conf_low[3L]))
  expect_true(all(is.finite(result$standard_error[1:2]) & result$standard_error[1:2] > 0))
})

test_that("a class carried by two units keeps its standard error", {
  # The single-unit rule must not catch a class whose weight is two units':
  # their scores are not zero at the class mean, and the sandwich is finite.
  set.seed(4)
  n <- 60L
  posteriors <- matrix(c(0.9, 0.1, 0), n, 3L, byrow = TRUE)
  posteriors[seq_len(n / 2), ] <- matrix(c(0.1, 0.9, 0), n / 2, 3L, byrow = TRUE)
  posteriors[c(n - 1L, n), ] <- matrix(c(0, 0, 1), 2L, 3L, byrow = TRUE)
  pieces <- list(posteriors = posteriors, modal = max.col(posteriors, "first"), n_classes = 3L)
  y <- stats::rnorm(n, c(0, 1, 3)[pieces$modal])
  warnings <- character()
  result <- withCallingHandlers(
    .multilpa_three_step_ml(pieces, y, "individuals", "equal", 0.95, "none", "BH",
                            "independent"),
    warning = \(w) {
      warnings <<- c(warnings, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  expect_false("latents_small_classes" %in% warnings)
  expect_true(all(is.finite(result$standard_error) & result$standard_error > 0))
})

test_that("a start that empties a class is dropped, not carried into the fit", {
  # sweep_041 of the equivalence kit (rounded): four classes on 50 units, one
  # unit alone in class 3 and no other unit with weight there. Its
  # proportional spread is the variance floor, and a tail start at that
  # spread empties the class; EM then stopped on an NA comparison. Varying
  # variances must end in the classed refusal, equal variances in a table.
  y <- c(
    585.4202, 549.9566, 439.4126, 432.4698, 462.8022, 635.3091, 545.8242, 524.5443,
    496.6691, 438.3774, 540.6471, 457.5127, 535.6670, 346.2529, 494.7356, 666.0787,
    527.8625, 470.2603, 435.3889, 530.8389, 595.3411, 569.7874, 645.0332, 422.9532,
    525.8397, 653.6453, 573.2465, 437.3385, 506.6312, 678.1226, 647.6448, 353.9873,
    468.2812, 466.6810, 424.9223, 447.0059, 617.4322, 475.8268, 644.3070, 486.9187,
    589.2064, 322.1966, 556.1361, 572.4241, 493.1043, 467.5396, 593.4719, 613.0681,
    636.4896, 556.4009)
  modal <- c(
    1L, 2L, 1L, 4L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 4L, 1L, 1L, 1L, 4L, 4L, 1L,
    2L, 1L, 2L, 4L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L, 1L, 4L, 1L, 1L, 1L, 4L,
    1L, 4L, 1L, 1L, 1L, 1L, 1L, 2L, 3L, 1L)
  posteriors <- matrix(c(
    1.00000e+00, 2.24894e-04, 9.99722e-01, 3.25200e-02, 9.99514e-01, 1.01586e-03,
    1.00000e+00, 1.00000e+00, 9.99998e-01, 9.94113e-01, 3.59705e-01, 9.98430e-01,
    9.91396e-01, 1.47356e-03, 9.96420e-01, 9.99999e-01, 9.99987e-01, 7.77247e-02,
    2.14927e-02, 9.98618e-01, 5.43447e-03, 1.00000e+00, 1.78931e-04, 6.29235e-04,
    1.00000e+00, 3.80967e-02, 9.99999e-01, 7.29607e-01, 9.99988e-01, 1.00000e+00,
    9.99998e-01, 1.00000e+00, 9.99877e-01, 3.39332e-04, 9.99740e-01, 3.07167e-04,
    1.00000e+00, 9.98005e-01, 9.99856e-01, 4.89069e-01, 1.00000e+00, 4.47859e-03,
    9.99978e-01, 1.00000e+00, 9.99999e-01, 9.99759e-01, 9.47126e-01, 4.50270e-03,
    3.49519e-18, 1.00000e+00, 1.17586e-10, 9.99775e-01, 4.99838e-10, 2.00982e-15,
    1.57054e-04, 9.98984e-01, 2.97084e-14, 1.67556e-07, 2.43256e-06, 1.37660e-04,
    6.40295e-01, 7.27258e-08, 8.60423e-03, 1.73900e-20, 3.58012e-03, 1.17491e-06,
    1.31169e-05, 5.31348e-11, 2.34112e-14, 4.31281e-10, 9.94566e-01, 1.86470e-07,
    9.99821e-01, 5.30686e-22, 2.84520e-11, 9.61903e-01, 1.34506e-06, 2.52954e-10,
    1.02142e-05, 9.11495e-10, 1.86223e-06, 7.42860e-08, 1.22633e-04, 8.39133e-25,
    2.60280e-04, 1.39149e-21, 2.22843e-08, 1.97876e-03, 1.43874e-04, 1.05301e-11,
    1.21518e-10, 4.65743e-15, 4.63149e-07, 9.95667e-09, 1.23621e-06, 1.46898e-04,
    5.28745e-02, 9.95497e-01, 8.94031e-15, 1.16515e-11, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
    0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.00000e+00, 0.00000e+00,
    5.17320e-27, 7.68654e-28, 2.78364e-04, 9.67480e-01, 3.28513e-04, 1.61494e-31,
    7.77196e-15, 3.68268e-11, 7.06414e-14, 5.74898e-03, 6.99332e-29, 1.57022e-03,
    4.75147e-15, 9.98526e-01, 1.59963e-09, 4.88664e-18, 2.30132e-27, 9.22275e-01,
    9.78507e-01, 1.38152e-03, 5.28668e-24, 1.25932e-15, 1.87639e-27, 9.99371e-01,
    4.40120e-37, 4.53887e-20, 1.12812e-17, 2.70393e-01, 1.42196e-06, 3.00776e-23,
    1.24613e-26, 2.40791e-26, 2.08304e-10, 9.99661e-01, 5.68702e-13, 9.99693e-01,
    1.68216e-22, 1.61222e-05, 3.38123e-07, 5.10931e-01, 4.11534e-18, 9.95521e-01,
    2.17591e-05, 2.16475e-23, 9.25520e-08, 9.37424e-05, 1.64484e-08, 8.40496e-18,
    1.10755e-54, 1.72432e-22), 50L, 4L)
  pieces <- list(posteriors = posteriors, modal = modal, n_classes = 4L)
  expect_error(suppressWarnings(.multilpa_three_step_ml(
    pieces, y, "individuals", "varying", 0.95, "none", "BH", "independent")),
    class = "latents_singular_information")
  equal <- suppressWarnings(.multilpa_three_step_ml(
    pieces, y, "individuals", "equal", 0.95, "none", "BH", "independent"))
  expect_true(all(is.finite(equal$estimate)))
})

test_that("the step-three starts neither read nor move the caller's random stream", {
  set.seed(7)
  n <- 120L
  posteriors <- matrix(c(0.85, 0.1, 0.05), n, 3L, byrow = TRUE)
  posteriors[seq_len(40L), ] <- matrix(c(0.1, 0.8, 0.1), 40L, 3L, byrow = TRUE)
  posteriors[81:120, ] <- matrix(c(0.05, 0.15, 0.8), 40L, 3L, byrow = TRUE)
  pieces <- list(posteriors = posteriors, modal = max.col(posteriors, "first"), n_classes = 3L)
  y <- stats::rt(n, 3) + c(0, 1, 3)[pieces$modal]
  set.seed(1)
  before <- .Random.seed
  first <- .multilpa_step_ml_fit(pieces, y, "varying")
  expect_identical(.Random.seed, before)
  set.seed(99)
  second <- .multilpa_step_ml_fit(pieces, y, "varying")
  expect_identical(first$means, second$means)
  expect_identical(first$log_likelihood, second$log_likelihood)
})

test_that("a class carried by one group has no clustered standard error", {
  # Perfect classification; class 3 is five observations of group 1 alone, so
  # its group-summed score is zero at its mean. Counted by observations it
  # looks carried by five units; the sandwich sums over groups.
  set.seed(7)
  classes <- c(rep(3L, 5), rep(1:2, length.out = 115))
  groups <- c(rep(1L, 5), rep(2:24, each = 5))
  posteriors <- diag(3)[classes, ]
  pieces <- list(posteriors = posteriors, modal = classes, units = seq_along(classes),
                 group_index = groups, n_classes = 3L)
  values <- stats::rnorm(120, c(0, 2, 5)[classes])
  expect_warning(
    clustered <- .multilpa_three_step_ml(pieces, values, "individuals", "equal", 0.95,
                                         "none", "BH", "cluster"),
    class = "latents_small_classes")
  expect_true(is.na(clustered$standard_error[3]))
  expect_true(all(is.finite(clustered$standard_error[1:2])))
})
