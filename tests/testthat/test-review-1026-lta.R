# Calibration and invariance checks of the latent transition engines against
# a brute-force sum over every state path (Rabiner, 1989): for a person with
# T occasions and K profiles, L = sum over the K^T paths s of
# pi(s_1) prod_t A_t(s_{t-1}, s_t) prod_t f(y_t | s_t), with an occasion that
# has no observed indicator contributing f = 1.
# The comparisons with independent references are in tests/equivalence/test-review-1026-lta.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

review_lta_data <- function(n = 30L, n_occasions = 3L, seed = 1L, stay = 0.75) {
  withr::with_seed(seed, {
    do.call(rbind, lapply(seq_len(n), function(i) {
      states <- Reduce(function(previous, t) {
        if (stats::runif(1L) < stay) previous else 3L - previous
      }, seq_len(n_occasions - 1L), init = sample.int(2L, 1L), accumulate = TRUE)
      data.frame(id = i, time = seq_len(n_occasions), z = stats::rnorm(1L),
                 w = stats::rnorm(n_occasions),
                 y1 = stats::rnorm(n_occasions, c(0, 2)[states], 0.7),
                 y2 = stats::rnorm(n_occasions, c(0, 1.5)[states], 0.7))
    }))
  })
}

review_softmax <- function(eta) exp(eta - max(eta)) / sum(exp(eta - max(eta)))

test_that("relabelled general-fit logits imply the permuted probabilities", {
  data <- review_lta_data(n = 30L, n_occasions = 3L, seed = 13L)
  fit <- suppressWarnings(lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                              time = "time", n_starts = 1, seed = 1,
                              transition_covariates = "w", initial_covariates = "z"))
  estimates <- .lta_bootstrap_estimates(fit)
  implied <- function(values, w, z) {
    pick <- function(name) values[[name]]
    moves <- t(vapply(1:2, function(k) {
      l <- 3L - k
      prefix <- sprintf("transition.group_class_1.profile_%d->profile_%d", k, l)
      eta <- numeric(2L)
      eta[l] <- pick(paste0(prefix, ".(Intercept)")) + w * pick(paste0(prefix, ".w"))
      review_softmax(eta)
    }, numeric(2)))
    start <- review_softmax(c(pick("initial.group_class_1.profile_1.(Intercept)") +
                                z * pick("initial.group_class_1.profile_1.z"), 0))
    list(moves = moves, start = start)
  }
  original <- implied(estimates, 0.4, -1.1)
  relabelled <- implied(.lta_relabel_profiles(estimates, c(2L, 1L)), 0.4, -1.1)
  expect_equal(relabelled$moves, original$moves[2:1, 2:1], tolerance = 1e-12)
  expect_equal(relabelled$start, original$start[2:1], tolerance = 1e-12)
})

test_that("an occasion with no observed indicator leaves the likelihood unchanged", {
  data <- review_lta_data(n = 30L, n_occasions = 3L, seed = 17L)
  trailing <- rbind(data, data.frame(id = 1:8, time = 4L, z = 0, w = 0,
                                     y1 = NA_real_, y2 = NA_real_))
  fit <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
             n_starts = 1, seed = 1, missing = "fiml")
  padded <- lta(trailing, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                time = "time", n_starts = 1, seed = 1, missing = "fiml")
  expect_equal(padded$log_likelihood, fit$log_likelihood, tolerance = 1e-8)
  expect_equal(padded$transition_probabilities, fit$transition_probabilities,
               tolerance = 1e-5)
  expect_identical(nobs(padded), nobs(fit))
  # Under the grid, a skipped position and an explicit all-missing row are
  # the same missing emission.
  skipped <- subset(data, !(id %in% 1:6 & time == 2))
  explicit <- data
  explicit[explicit$id %in% 1:6 & explicit$time == 2, c("y1", "y2")] <- NA
  grid_skip <- lta(skipped, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                   time = "time", n_starts = 1, seed = 1, missing = "fiml",
                   occasions = "grid")
  grid_na <- lta(explicit, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                 time = "time", n_starts = 1, seed = 1, missing = "fiml",
                 occasions = "grid")
  expect_equal(grid_na$log_likelihood, grid_skip$log_likelihood, tolerance = 1e-8)
  expect_equal(BIC(grid_na), grid_na$bic, tolerance = 1e-12)
  expect_equal(BIC(grid_na), BIC(grid_skip), tolerance = 1e-6)
})

test_that("nested transition comparisons refuse different fitting bounds", {
  data <- review_lta_data(n = 25L, n_occasions = 3L, seed = 29L)
  null <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
              n_starts = 1, seed = 1, variance_model = "equal")
  alternative <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                     time = "time", n_starts = 1, seed = 1, min_variance = 1e-4)
  same_bounds <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                     time = "time", n_starts = 1, seed = 1)
  expect_true(.lta_nested(null, same_bounds))
  expect_false(.lta_nested(null, alternative))
  expect_error(bootstrap_lrt(null, alternative, data = data, iter = 2L, seed = 1),
               class = "latents_bad_argument")
})
