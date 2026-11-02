skip_on_cran()

# Audit 2026-11-01, slice `lta`: regression tests for the fixes and the
# independent checks of the transition engine.

audit_panel <- function(n = 60L, n_occasions = 3L, seed = 1L) {
  set.seed(seed)
  move <- matrix(c(0.8, 0.2, 0.25, 0.75), 2L, byrow = TRUE)
  first <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  states <- Reduce(function(previous, t) {
    vapply(previous, function(k) sample(1:2, 1L, prob = move[k, ]), integer(1))
  }, seq_len(n_occasions - 1L), first, accumulate = TRUE)
  state <- as.vector(t(do.call(cbind, states)))
  data <- data.frame(id = rep(seq_len(n), each = n_occasions),
                     time = rep(seq_len(n_occasions), n), state = state)
  centres <- rbind(c(0, 0), c(2.5, 2))
  data$y1 <- stats::rnorm(nrow(data), centres[data$state, 1L])
  data$y2 <- stats::rnorm(nrow(data), centres[data$state, 2L])
  data$z <- stats::rnorm(nrow(data))
  data$w <- stats::rnorm(nrow(data))
  data$g <- rep(stats::rnorm(n), each = n_occasions)
  data
}

# The likelihood summed over every latent path, written from the model's
# definition: initial logits against the last profile, transition logits
# against staying, covariates read at the destination occasion.
audit_path_sum <- function(fit, data, vars, transition_covariates, initial_covariates) {
  n_profiles <- fit$n_profiles
  softmax <- function(eta) exp(eta - max(eta)) / sum(exp(eta - max(eta)))
  per_group <- vapply(unique(data$id), function(g) {
    rows <- data[data$id == g, ]
    rows <- rows[order(rows$time), ]
    n_occasions <- nrow(rows)
    paths <- as.matrix(expand.grid(rep(list(seq_len(n_profiles)), n_occasions)))
    block <- fit$measurement[[1L]]
    density <- vapply(seq_len(n_occasions), function(t) {
      vapply(seq_len(n_profiles), function(k) {
        prod(stats::dnorm(unlist(rows[t, vars]), block$means[k, vars],
                          sqrt(block$variances[k, vars])))
      }, numeric(1))
    }, numeric(n_profiles))
    by_class <- vapply(seq_len(length(fit$group_probabilities)), function(h) {
      start <- c(1, unlist(rows[1L, initial_covariates]))
      initial <- softmax(c(drop(start %*% matrix(fit$initial_coefficients[, , h],
                                                  length(start))), 0))
      step <- function(t, k) {
        z <- c(1, unlist(rows[t, transition_covariates]))
        eta <- drop(z %*% matrix(fit$transition_coefficients[, , k, h], length(z)))
        p <- softmax(c(eta, 0))
        out <- numeric(n_profiles)
        out[setdiff(seq_len(n_profiles), k)] <- p[-n_profiles]
        out[k] <- p[n_profiles]
        out
      }
      sum(apply(paths, 1L, function(path) {
        initial[path[1L]] * density[path[1L], 1L] *
          prod(vapply(seq_len(n_occasions)[-1L], function(t) {
            step(t, path[t - 1L])[path[t]] * density[path[t], t]
          }, numeric(1)))
      }))
    }, numeric(1))
    log(sum(fit$group_probabilities * by_class))
  }, numeric(1))
  sum(per_group)
}

test_that("the forward-backward likelihood equals the sum over every latent path", {
  data <- audit_panel(n = 40L)
  fit <- suppressWarnings(lta(data, c("y1", "y2"), "id", 2L, "time",
                              transition_covariates = "z", initial_covariates = "g",
                              n_group_classes = 2L, n_starts = 2L, seed = 1L))
  expect_s3_class(fit, "multilpa_lta")
  brute <- audit_path_sum(fit, data, c("y1", "y2"), "z", "g")
  expect_equal(fit$log_likelihood, brute, tolerance = 1e-10)
  # Rows in another order describe the same sequences.
  shuffled <- data[sample(nrow(data)), ]
  refit <- suppressWarnings(lta(shuffled, c("y1", "y2"), "id", 2L, "time",
                                transition_covariates = "z", initial_covariates = "g",
                                n_group_classes = 2L, n_starts = 2L, seed = 1L))
  expect_equal(audit_path_sum(refit, shuffled, c("y1", "y2"), "z", "g"),
               refit$log_likelihood, tolerance = 1e-10)
})

test_that("reported transition rows are probability distributions", {
  data <- audit_panel(n = 40L)
  fit <- suppressWarnings(lta(data, c("y1", "y2"), "id", 2L, "time",
                              transitions = "occasion", transition_covariates = "z",
                              n_starts = 2L, seed = 1L))
  table <- get_results(fit, "transitions")
  sums <- stats::aggregate(probability ~ group_class + occasion + from, table, sum)
  expect_equal(sums$probability, rep(1, nrow(sums)), tolerance = 1e-12)
  expect_true(all(table$probability >= 0))
})

test_that("bootstrap controls of a transition fit are refused with a class", {
  data <- audit_panel(n = 30L)
  general <- lta(data, c("y1", "y2"), "id", 2L, "time", transitions = "occasion",
                 n_starts = 2L, seed = 1L)
  expect_error(parameter_inference(general, data = data, method = "bootstrap", iter = 1),
               class = "latents_bad_argument")
  expect_error(parameter_inference(general, data = data, method = "bootstrap",
                                   iter = 10, n_starts = 0),
               class = "latents_bad_argument")
  homogeneous <- lta(data, c("y1", "y2"), "id", 2L, "time", n_starts = 2L, seed = 1L)
  expect_error(parameter_inference(homogeneous, data = data, method = "bootstrap",
                                   iter = 10, tol = -1),
               class = "latents_bad_argument")
})

test_that("bootstrap_lrt() refuses two fits whose shared covariate differs", {
  data <- audit_panel(n = 40L)
  changed <- data
  changed$z <- 3 * data$z + 1 + stats::rnorm(nrow(data), sd = 0.5)
  null <- lta(data, c("y1", "y2"), "id", 2L, "time", transition_covariates = "z",
              n_starts = 2L, seed = 1L)
  alternative <- lta(changed, c("y1", "y2"), "id", 2L, "time",
                     transition_covariates = c("z", "w"), n_starts = 2L, seed = 1L)
  # Neither data set reproduces both fits, so no replicate is ever drawn.
  expect_error(bootstrap_lrt(null, alternative, data = data, iter = 2L, seed = 1L),
               class = "latents_bad_inference_data")
  expect_error(bootstrap_lrt(null, alternative, data = changed, iter = 2L, seed = 1L),
               class = "latents_bad_inference_data")
})

test_that("get_group_tna() refuses a bad label with a class", {
  skip_if_not_installed("tna")
  data <- audit_panel(n = 30L)
  fit <- lta(data, c("y1", "y2"), "id", 2L, "time", n_starts = 2L, seed = 1L)
  expect_error(get_group_tna(fit, label = 3), class = "latents_bad_argument")
  expect_error(get_group_tna(fit, label = NA_character_), class = "latents_bad_argument")
  profiles <- multilpa(data, c("y1", "y2"), n_profiles = 2L, id = "id", time = "time",
                       n_starts = 2L, seed = 1L)
  expect_error(get_group_tna(profiles, label = c("a", "b")), class = "latents_bad_argument")
})

test_that("a transition enumeration refuses unknown tables and bad counts with a class", {
  data <- audit_panel(n = 30L)
  grid <- suppressWarnings(enumerate_classes(data, c("y1", "y2"), "id", n_profiles = 2L,
                                             time = "time", n_starts = 1L, seed = 1L))
  expect_s3_class(get_results(grid, "best"), "data.frame")
  expect_error(get_results(grid, "bset"), class = "latents_bad_argument")
  expect_error(get_results(grid, criterion = "bci"), class = "latents_bad_argument")
  expect_error(enumerate_classes(data, c("y1", "y2"), "id", n_profiles = c(2, NA),
                                 time = "time", n_starts = 1L, seed = 1L),
               class = "latents_bad_argument")
  expect_error(enumerate_classes(data, c("y1", "y2"), "id", n_profiles = "2",
                                 time = "time", n_starts = 1L, seed = 1L),
               class = "latents_bad_argument")
  expect_error(enumerate_classes(data, c("y1", "y2"), "id", n_profiles = integer(0),
                                 time = "time", n_starts = 1L, seed = 1L),
               class = "latents_bad_argument")
  expect_error(enumerate_classes(data, c("y1", "y2"), "id", n_profiles = 2L,
                                 n_group_classes = 0L, time = "time", n_starts = 1L,
                                 seed = 1L),
               class = "latents_bad_argument")
})
