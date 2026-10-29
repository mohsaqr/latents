# Comparisons with independent references moved from tests/testthat/test-lta-three-step.R,
# so that they run in CI's equivalence job rather than under R CMD check.

markov_fixture <- function(n = 200L, occasions = 4L, seed = 12L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  w <- stats::rnorm(n)
  z <- matrix(stats::rnorm(n * occasions), n)
  state <- matrix(NA_integer_, n, occasions)
  state[, 1L] <- ifelse(stats::runif(n) < stats::plogis(0.3 + 0.8 * w), 1L, 2L)
  # Sequential by construction: each occasion's state depends on the last.
  for (t in seq_len(occasions)[-1L]) {  # simulation only
    logit <- ifelse(state[, t - 1L] == 1L, -1.2 + 0.8 * z[, t], -1.0 - 0.6 * z[, t])
    moves <- stats::runif(n) < stats::plogis(logit)
    state[, t] <- ifelse(moves, 3L - state[, t - 1L], state[, t - 1L])
  }
  data <- data.frame(id = rep(seq_len(n), each = occasions),
                     time = rep(seq_len(occasions), n), w = rep(w, each = occasions),
                     z = as.vector(t(z)))
  s <- as.vector(t(state))
  data$y1 <- stats::rnorm(length(s), c(-1.2, 1.2)[s], 1)
  data$y2 <- stats::rnorm(length(s), c(-1, 1)[s], 1)
  data
}

markov_fit <- function(data, ...) {
  lta(data, c("y1", "y2"), "id", n_profiles = 2, time = "time", n_starts = 2,
      seed = 1, ...)
}

test_that("the step-three likelihood equals a sum over every path", {
  data <- markov_fixture()
  fit <- markov_fit(data)
  result <- r3step(fit, data, initial_covariates = "w", transition_covariates = "z")
  expect_named(result, c("block", "from", "to", "term", "estimate", "standard_error",
                         "statistic", "p_value", "p_adjusted", "conf_low", "conf_high"))
  errors <- attr(result, "classification_errors")
  expect_equal(unname(rowSums(errors)), c(1, 1), tolerance = 1e-12)
  # Each occasion's modal profile from its own indicators: pooled shares times
  # the fitted normal densities of both indicators, written out here.
  shares <- colMeans(fit$subject_posteriors)
  joint <- vapply(1:2, function(k) {
    shares[k] * stats::dnorm(data$y1, fit$means[k, "y1"], sqrt(fit$variances[k, "y1"])) *
      stats::dnorm(data$y2, fit$means[k, "y2"], sqrt(fit$variances[k, "y2"]))
  }, numeric(nrow(data)))
  posteriors <- joint / rowSums(joint)
  expect_equal(unname(latents:::.lta_occasion_posteriors(fit)), unname(posteriors),
               tolerance = 1e-10)
  modal <- max.col(posteriors, ties.method = "first")
  expect_equal(unname(errors), unname(t(vapply(1:2, function(k) {
    vapply(1:2, function(j) sum(posteriors[modal == j, k]), numeric(1)) / sum(posteriors[, k])
  }, numeric(2)))), tolerance = 1e-10)
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  coefficient <- function(block, from, term, theta) {
    theta[result$block == block & (is.na(result$from) | result$from == from) &
            result$term == term]
  }
  path_likelihood <- function(theta) {
    sum(vapply(split(seq_len(nrow(data)), data$id), function(rows) {
      w <- data$w[rows[1L]]
      first <- stats::plogis(coefficient("initial", NA, "(Intercept)", theta) +
                               coefficient("initial", NA, "w", theta) * w)
      initial <- c(first, 1 - first)
      move <- function(from, t) {
        stats::plogis(coefficient("transition", paste0("profile_", from), "(Intercept)", theta) +
                        coefficient("transition", paste0("profile_", from), "z", theta) *
                        data$z[rows[t]])
      }
      log(sum(apply(paths, 1L, function(path) {
        moves <- vapply(2:4, function(t) {
          p <- move(path[t - 1L], t)
          if (path[t] == path[t - 1L]) 1 - p else p
        }, numeric(1))
        initial[path[1L]] * prod(moves) * prod(errors[cbind(path, modal[rows])])
      })))
    }, numeric(1)))
  }
  expect_equal(attr(result, "log_likelihood"), path_likelihood(result$estimate),
               tolerance = 1e-9)
  # The estimate is a maximum of that likelihood: its gradient vanishes.
  gradient <- vapply(seq_along(result$estimate), function(j) {
    step <- replace(numeric(nrow(result)), j, 1e-5)
    (path_likelihood(result$estimate + step) - path_likelihood(result$estimate - step)) / 2e-5
  }, numeric(1))
  expect_lt(max(abs(gradient)), 1e-3)
})
