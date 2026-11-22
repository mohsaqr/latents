# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

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
  pieces <- latents:::.multilpa_level_assignments(fit, "individuals")
  steps <- lapply(c(varying = "varying", equal = "equal"),
                  \(variance) latents:::.multilpa_step_ml_fit(pieces, data$y, variance))
  ml_note <- "Package side: the internal step-three ML fit (.multilpa_step_ml_fit) behind three_step(method = 'ml')."
  # The same likelihood summed twice; only the order of the additions differs.
  expect_reference(vapply(steps, \(step) step$log_likelihood, numeric(1)),
                   vapply(steps, \(step) log_likelihood(step$means, sqrt(step$variances)), numeric(1)),
                   reference = "independent: Vermunt (2010) fixed-error step-three likelihood, base R",
                   covers = "three_step()", quantity = "step-three log likelihood by variance model",
                   tolerance = 1e-10, measure = "relative", note = ml_note)
  # The estimate is stationary for this likelihood: a step of 1e-4 standard
  # deviations in any mean changes it only at second order. Each entry is the
  # central-difference slope in standard-deviation units; the original bound
  # on its size, 1e-4 absolute.
  slopes <- vapply(steps, \(step) {
    sds <- sqrt(step$variances)
    vapply(1:3, \(t) {
      h <- 1e-4 * sds[t]
      up <- step$means
      up[t] <- up[t] + h
      down <- step$means
      down[t] <- down[t] - h
      (log_likelihood(up, sds) - log_likelihood(down, sds)) / (2 * h) * sds[t]
    }, numeric(1))
  }, numeric(3))
  rownames(slopes) <- sprintf("class_%d", 1:3)
  expect_reference(slopes, matrix(0, 3L, 2L, dimnames = dimnames(slopes)),
                   reference = "independent: Vermunt (2010) fixed-error step-three likelihood, base R",
                   covers = "three_step()",
                   quantity = "standardised slope of the likelihood at the step-three means",
                   tolerance = 1e-4, measure = "absolute", note = ml_note)
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
  step <- latents:::.multilpa_step_ml_fit(pieces, case$y, "varying")
  expect_gt(step$log_likelihood, -880.5)
  # The best of 20 BFGS maximisations of the likelihood written above. The
  # original bound was one-sided, 1e-6 * (1 + |best|) below the best; 1e-6
  # relative on both sides is no looser.
  expect_reference(step$log_likelihood, best,
                   reference = "independent: stats::optim() over 20 random starts of the step-three likelihood",
                   covers = "three_step()", quantity = "step-three log likelihood, heavy-tailed outcome",
                   tolerance = 1e-6, measure = "relative",
                   note = "Package side: the internal step-three ML fit (.multilpa_step_ml_fit) behind three_step(method = 'ml').")
})
