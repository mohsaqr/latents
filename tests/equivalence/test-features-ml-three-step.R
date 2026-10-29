# Comparisons with independent references moved from tests/testthat/test-features-ml-three-step.R,
# so that they run in CI's equivalence job rather than under R CMD check.

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
