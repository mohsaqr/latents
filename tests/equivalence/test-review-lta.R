# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

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

# The group log likelihoods summed over every profile path, written from the
# model's definition: softmax links, Gaussian, categorical, adjacent-category
# ordinal and Poisson / negative-binomial densities. The coefficient vector is
# split into its blocks by the package's internal .lta_unpack().
review_lta_exact <- function(theta, fit) {
  p <- latents:::.lta_unpack(theta, fit)
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
          if (!is.na(k)) d <- d * if (identical(fit$extra_data$count_model, "negative_binomial"))
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
  settings <- c("two group classes", "occasion transitions and measurement",
                "second order, occasion transitions, mover-stayer",
                "negative binomial, equal dispersion",
                "negative binomial, varying dispersion, weighted",
                "skipped grid positions, occasion transitions, two group classes")
  # One row per group (and parameter) and variant: the package's values and
  # the path-sum references, compared below in one table each.
  cells <- function(values, setting) {
    values <- as.matrix(values)
    data.frame(setting = setting, group = as.vector(row(values)),
               parameter = as.vector(col(values)), value = as.vector(values))
  }
  results <- Map(function(arguments, setting) {
    fit <- do.call(review_lta_fit, arguments)
    theta <- latents:::.lta_theta(fit)
    theta[] <- sin(seq_along(theta)) / 3
    # Evaluate away from an optimized point: no stationary score can hide a
    # missing block. Finite coordinates keep all exact-path masses positive.
    exact <- review_lta_exact(theta, fit)
    expectation <- latents:::.lta_expectation(fit$x, fit$codes, fit$layout, fit$designs,
      latents:::.lta_unpack(theta, fit), fit$occasion_of_row, extra = fit$extra_data)
    numeric_scores <- vapply(seq_along(theta), function(k) {
      plus <- minus <- theta
      plus[k] <- plus[k] + 1e-5
      minus[k] <- minus[k] - 1e-5
      (review_lta_exact(plus, fit) - review_lta_exact(minus, fit)) / 2e-5
    }, numeric(fit$n_groups)) * (unname(fit$sampling_weights) %||% 1)
    # Unit checks: the coefficient vector round-trips and the data guard holds.
    expect_equal(length(theta), fit$n_parameters)
    expect_equal(latents:::.lta_pack(latents:::.lta_unpack(theta, fit), fit$variance_model,
                                     fit$extra_data),
                 unname(theta), tolerance = 1e-12)
    expect_invisible(latents:::.lta_check_data(fit, review_lta_data(isTRUE(arguments$gaps))))
    changed <- review_lta_data(isTRUE(arguments$gaps))
    changed$k[1L] <- changed$k[1L] + 1
    expect_error(latents:::.lta_check_data(fit, changed), class = "latents_bad_inference_data")
    list(likelihood = cells(unname(expectation$group_log_likelihood), setting),
         exact = cells(exact, setting),
         scores = cells(unname(latents:::.lta_group_scores(theta, fit)), setting),
         numeric_scores = cells(numeric_scores, setting))
  }, variants, settings)
  stack <- function(part) do.call(rbind, unname(lapply(results, `[[`, part)))
  expect_reference(stack("likelihood"), stack("exact"),
                   reference = "independent: exact sum over every profile path, base R (coefficients split by the package's .lta_unpack)",
                   covers = "lta()", quantity = "group log likelihood away from the optimum",
                   tolerance = 1e-11, measure = "relative", key = c("setting", "group", "parameter"),
                   note = "Package side: the internal forward-backward likelihood (.lta_expectation) of lta().")
  # Central differences (step 1e-5) of the exact path likelihood. Scaled by
  # max(1, |reference|): a group's score on a parameter it barely informs is
  # near zero, where the differences' absolute error, not a relative one,
  # is what can be judged.
  expect_reference(stack("scores"), stack("numeric_scores"),
                   reference = "independent: central differences of the exact path likelihood, base R (coefficients split by .lta_unpack)",
                   covers = "lta()", quantity = "group scores away from the optimum",
                   tolerance = 2e-7, measure = "scaled", key = c("setting", "group", "parameter"),
                   note = "Package side: the internal analytic scores (.lta_group_scores) behind vcov() of lta().")
})
