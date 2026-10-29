# Comparisons with independent references moved from tests/testthat/test-review-lta.R,
# so that they run in CI's equivalence job rather than under R CMD check.

# Rows 2 and 11 sit inside their groups' spans. By default they are kept with
# every indicator missing, an occasion without measurement information whose
# transition covariates are known; `gaps = TRUE` drops them, so the grid skips
# those positions -- which lta() refuses with transition covariates, since a
# skipped position has no covariate value.
review_lta_data <- function(gaps = FALSE) {
  data <- data.frame(id = rep(1:6, each = 4L), time = rep(1:4, 6L),
                     y = sin(seq_len(24L)), z = cos(seq_len(24L)),
                     w = rep(seq_len(6L) / 7, each = 4L),
                     c = factor(rep(c("no", "yes"), 12L)),
                     o = ordered(rep(c("low", "mid", "high"), 8L),
                                 levels = c("low", "mid", "high")),
                     k = rep(c(0, 1, 3, 5), 6L), weight = rep(1:3, each = 8L))
  data[c(3L, 18L), "y"] <- NA_real_
  data[c(7L, 14L), "o"] <- NA
  data[c(8L, 17L), "k"] <- NA_real_
  if (gaps) return(data[-c(2L, 11L, 23L, 24L), ])
  data[c(2L, 11L), c("y", "c", "o", "k")] <- NA
  data[-c(23L, 24L), ]
}

review_lta_fit <- function(..., gaps = FALSE,
                           transition_covariates = if (gaps) character() else "z") {
  quietly(lta(review_lta_data(gaps), c("y", "c", "o", "k"), "id", 3L,
    time = "time", categorical = "c", ordinal = "o", count = "k", missing = "fiml",
    transition_covariates = transition_covariates, initial_covariates = "w",
    occasions = "grid",
    n_starts = 1L, max_iter = 0L, seed = 51L, ...),
          c("latents_unconverged", "latents_boundary", "latents_small_classes",
            "latents_thin_profiles"))
}

review_lta_exact <- function(theta, fit) {
  p <- .lta_unpack(theta, fit)
  softmax <- function(v) exp(v - max(v)) / sum(exp(v - max(v)))
  K <- fit$n_profiles
  vapply(seq_len(fit$n_groups), function(g) {
    T <- fit$layout$span[g]
    paths <- as.matrix(expand.grid(rep(list(seq_len(K)), T)))
    class_likelihood <- vapply(seq_len(fit$n_group_classes), function(h) {
      initial <- softmax(c(drop(fit$designs$initial[g, , drop = FALSE] %*%
        matrix(p$initial[, , h], ncol(fit$designs$initial))), 0))
      sum(apply(paths, 1L, function(path) {
        moves <- vapply(seq_len(T)[-1L], function(t) {
          previous <- path[t - 1L]
          if (fit$stayer[h]) return(as.numeric(path[t] == previous))
          pair <- if (fit$order >= 2L && t >= 3L)
            (path[t - 2L] - 1L) * K + previous else previous
          coefficients <- if (fit$order >= 2L && t >= 3L) p$transition2 else p$transition
          logits <- drop(fit$designs$transition[[t - 1L]][g, , drop = FALSE] %*%
            matrix(coefficients[, , pair, h], ncol(fit$designs$transition[[t - 1L]])))
          destinations <- c(setdiff(seq_len(K), previous), previous)
          softmax(c(logits, 0))[match(path[t], destinations)]
        }, numeric(1))
        density <- vapply(seq_len(T), function(t) {
          row <- fit$layout$slot[g, t]
          if (is.na(row)) return(1)
          b <- p$measurement[[if (length(p$measurement) == 1L) 1L else t]]
          state <- path[t]
          y <- fit$x[row, 1L]
          d <- if (is.na(y)) 1 else dnorm(y, b$means[state, 1L], sqrt(b$variances[state, 1L]))
          c <- fit$codes[row, 1L]
          if (!is.na(c)) d <- d * b$response_probabilities[[1L]][state, c]
          o <- fit$extra_data$ordinal[row, 1L]
          if (!is.na(o)) {
            logits <- c(0, b$ordinal_intercepts[[1L]]) +
              (seq_len(3L) - 1L) * b$ordinal_locations[state, 1L]
            d <- d * softmax(logits)[o]
          }
          k <- fit$extra_data$count[row, 1L]
          if (!is.na(k)) d <- d * if (.latents_negative_binomial(fit$extra_data))
            dnbinom(k, size = 1 / b$count_dispersion[state, 1L], mu = b$count_means[state, 1L]) else
              dpois(k, b$count_means[state, 1L])
          d
        }, numeric(1))
        initial[path[1L]] * prod(moves) * prod(density)
      }))
    }, numeric(1))
    log(sum(class_likelihood * p$group_probabilities))
  }, numeric(1))
}

test_that("all general LTA extensions match independent path likelihoods and scores", {
  skip_on_cran()
  variants <- list(
    list(n_group_classes = 2L),
    list(transitions = "occasion", measurement = "occasion"),
    list(order = 2L, transitions = "occasion", mover_stayer = TRUE),
    list(count_model = "negative_binomial", count_dispersion = "equal"),
    list(count_model = "negative_binomial", count_dispersion = "varying", weights = "weight"),
    # Skipped grid positions (no row, no covariate) with occasion transitions.
    list(gaps = TRUE, transitions = "occasion", n_group_classes = 2L))
  lapply(variants, function(arguments) {
    fit <- do.call(review_lta_fit, arguments)
    theta <- .lta_theta(fit)
    theta[] <- sin(seq_along(theta)) / 3
    # Evaluate away from an optimized point: no stationary score can hide a
    # missing block. Finite coordinates keep all exact-path masses positive.
    exact <- review_lta_exact(theta, fit)
    expectation <- .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs,
      .lta_unpack(theta, fit), fit$occasion_of_row, extra = fit$extra_data)
    expect_equal(unname(expectation$group_log_likelihood), exact, tolerance = 1e-11)
    numeric_scores <- vapply(seq_along(theta), function(k) {
      plus <- minus <- theta
      plus[k] <- plus[k] + 1e-5
      minus[k] <- minus[k] - 1e-5
      (review_lta_exact(plus, fit) - review_lta_exact(minus, fit)) / 2e-5
    }, numeric(fit$n_groups)) * (unname(fit$sampling_weights) %||% 1)
    expect_equal(unname(.lta_group_scores(theta, fit)), numeric_scores, tolerance = 2e-7)
    expect_equal(length(theta), fit$n_parameters)
    expect_equal(.lta_pack(.lta_unpack(theta, fit), fit$variance_model, fit$extra_data),
                 unname(theta), tolerance = 1e-12)
    expect_invisible(.lta_check_data(fit, review_lta_data(isTRUE(arguments$gaps))))
    changed <- review_lta_data(isTRUE(arguments$gaps))
    changed$k[1L] <- changed$k[1L] + 1
    expect_error(.lta_check_data(fit, changed), class = "latents_bad_inference_data")
  })
})
