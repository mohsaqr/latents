# Regression tests for the round-2 latent transition findings of the
# 2026-10-26 review (validation/review-2026-10-26/round-2/lta.md):
# R2-LTA-1 zero-weight units in the drop_empty = FALSE sample sizes;
# R2-LTA-2 Wald inference over identified coordinates only;
# R2-LTA-3 second-order logits with nothing observed from occasion 3 on;
# R2-LTA-4 the homogeneous engine with nothing observed after occasion 1;
# R2-LTA-5 the rank check reads only rows that enter the likelihood;
# R2-LTA-6 interior empty occasions count only the identified product;
# and the sampling-weighted entropy of a transition enumeration grid.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

r2_lta_data <- function(n = 80L, n_occasions = 4L, seed = 21L) {
  withr::with_seed(seed, {
    do.call(rbind, lapply(seq_len(n), function(i) {
      states <- Reduce(function(previous, t) {
        if (stats::runif(1L) < 0.7) previous else 3L - previous
      }, seq_len(n_occasions - 1L), init = sample.int(2L, 1L), accumulate = TRUE)
      data.frame(id = i, time = seq_len(n_occasions), z = stats::rnorm(1L),
                 w = stats::rnorm(n_occasions),
                 y1 = stats::rnorm(n_occasions, c(0, 2)[states], 0.7),
                 y2 = stats::rnorm(n_occasions, c(0, 1.5)[states], 0.7))
    }))
  })
}

r2_lta_fit <- function(data, ...) {
  quietly(lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
              n_starts = 3, seed = 1, missing = "fiml", ...))
}

# Ten groups gain a fifth occasion at which nothing is observed.
r2_lta_trailing <- function(data) {
  trailing <- data[data$time == 4L & data$id <= 20L, ]
  trailing$time <- 5L
  trailing[c("y1", "y2")] <- NA
  rbind(data, trailing)
}

# Hand count, K = 2, two Gaussian indicators with varying variances (8
# measurement parameters per block), one group class: initial logits p0,
# first-order transition logits 2 p, second-order logits 4 p2.
r2_lta_count <- function(p0 = 1L, p = 1L, p2 = 0L, blocks = 1L) {
  as.integer(blocks * 8L + p0 + 2L * p + 4L * p2)
}

# Numerical rank of the Hessian of the package's own log likelihood in the
# packed coordinates: the number of identified directions.
r2_lta_hessian_rank <- function(fit) {
  theta <- .lta_pack(.lta_parameters(fit), fit$variance_model, fit$extra_data)
  log_likelihood <- function(v) {
    .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs, .lta_unpack(v, fit),
                     fit$occasion_of_row, unname(fit$sampling_weights),
                     fit$extra_data)$log_likelihood
  }
  values <- eigen(-numDeriv::hessian(log_likelihood, theta), symmetric = TRUE,
                  only.values = TRUE)$values
  sum(values > 1e-6 * max(values))
}

# R2-LTA-1 -------------------------------------------------------------------

test_that("zero-weight groups are in no drop_empty = FALSE sample size", {
  data <- r2_lta_data(n = 50L, seed = 41L)
  data$sampling <- ifelse(data$id %in% 4:8, 0, 1)
  fit <- quietly(lta(data, c("y1", "y2"), "id", 2, "time", weights = "sampling",
                     drop_empty = FALSE, n_starts = 2, seed = 1))
  without <- quietly(lta(subset(data, sampling > 0), c("y1", "y2"), "id", 2, "time",
                         weights = "sampling", drop_empty = FALSE, n_starts = 2,
                         seed = 1))
  expect_equal(fit$log_likelihood, without$log_likelihood, tolerance = 1e-8)
  expect_equal(fit$bic, BIC(fit))
  expect_equal(fit$bic, without$bic, tolerance = 1e-8)
  expect_equal(fit$bic_individual, without$bic_individual, tolerance = 1e-8)
  expect_identical(fit$n_informative, 45L * 4L)
  expect_equal(get_results(fit, "fit")$bic, without$bic, tolerance = 1e-8)
})

# R2-LTA-2 -------------------------------------------------------------------

test_that("a flat trailing occasion leaves the other logits their standard errors", {
  data <- r2_lta_trailing(r2_lta_data())
  fit <- r2_lta_fit(data, transitions = "occasion")
  identified <- .lta_identified_coordinates(fit)
  expect_identical(sum(identified), fit$n_parameters)
  expect_identical(fit$n_parameters, r2_lta_count(p = 3L))
  table <- parameter_inference(fit)
  flat <- table$term == "occasion_5"
  expect_identical(sum(flat), 2L)
  expect_false(any(table$estimated[flat]))
  expect_true(all(is.na(table$standard_error[flat])))
  expect_true(all(table$estimated[!flat]))
  expect_true(all(is.finite(table$standard_error[!flat]) & table$standard_error[!flat] > 0))
  covariance <- vcov(fit)
  expect_true(all(is.na(covariance[!identified, ])))
  expect_true(all(is.finite(covariance[identified, identified])))
  coefficients <- get_results(fit, "transition_coefficients")
  expect_identical(coefficients$estimated, coefficients$term != "occasion_5")
  # The standard errors of the identified logits are those of the fit without
  # the empty occasion: the flat rows add nothing to the likelihood.
  base <- r2_lta_fit(r2_lta_data(), transitions = "occasion")
  reference <- parameter_inference(base)
  expect_equal(table$standard_error[!flat], reference$standard_error, tolerance = 1e-4)
})

test_that("the measurement block of an empty occasion is left out of inference", {
  data <- r2_lta_trailing(r2_lta_data())
  fit <- r2_lta_fit(data, measurement = "occasion")
  identified <- .lta_identified_coordinates(fit)
  expect_identical(sum(identified), fit$n_parameters)
  expect_identical(sum(!identified), 8L)
  expect_true(all(grepl("\\.occasion_5\\.", names(coef(fit))[!identified])))
  table <- parameter_inference(fit)
  expect_true(all(is.finite(table$standard_error)))
})

# R2-LTA-3 -------------------------------------------------------------------

test_that("second-order logits with nothing observed from occasion 3 on are not counted", {
  data <- subset(r2_lta_data(n = 150L, n_occasions = 3L), time <= 3)
  data[data$time == 3L, c("y1", "y2")] <- NA
  homogeneous <- r2_lta_fit(data, order = 2L)
  expect_identical(homogeneous$n_parameters, r2_lta_count(p = 1L))
  expect_identical(sum(.lta_identified_coordinates(homogeneous)), homogeneous$n_parameters)
  table <- parameter_inference(homogeneous)
  expect_false(any(table$estimated[table$block == "transition2"]))
  expect_true(all(is.finite(table$standard_error[table$block != "transition2"])))
  skip_if_not_installed("numDeriv")
  expect_identical(r2_lta_hessian_rank(homogeneous), homogeneous$n_parameters)
})

# R2-LTA-4 -------------------------------------------------------------------

test_that("the homogeneous engine refuses data with nothing observed after occasion 1", {
  data <- r2_lta_data()
  data[data$time > 1L, c("y1", "y2")] <- NA
  expect_error(r2_lta_fit(data), class = "latents_bad_transition")
  # The general engine has always refused the same data.
  expect_error(r2_lta_fit(data, transitions = "occasion"), class = "latents_bad_transition")
})

# R2-LTA-5 -------------------------------------------------------------------

test_that("a covariate varying only in rows outside the likelihood is refused", {
  data <- r2_lta_data()
  data$sampling <- ifelse(data$id <= 10L, 0, 1)
  data$moving <- ifelse(data$id <= 10L, data$w, 0.3)
  data$starting <- ifelse(data$id <= 10L, data$z, 0.3)
  expect_error(r2_lta_fit(data, weights = "sampling", transition_covariates = "moving"),
               class = "latents_bad_covariate")
  expect_error(r2_lta_fit(data, weights = "sampling", initial_covariates = "starting"),
               class = "latents_bad_covariate")
  # Varying only on the moves into an occasion at which nothing is observed.
  padded <- r2_lta_trailing(r2_lta_data())
  padded$late <- ifelse(padded$time == 5L, padded$w, 0.5)
  expect_error(r2_lta_fit(padded, transition_covariates = "late"),
               class = "latents_bad_covariate")
  # A covariate that varies where the likelihood reads it is still accepted,
  # with a trailing empty occasion and occasion transitions.
  padded$moving <- padded$w
  fit <- r2_lta_fit(padded, transitions = "occasion", transition_covariates = "moving")
  expect_identical(fit$n_parameters, r2_lta_count(p = 4L))
})

# R2-LTA-6 -------------------------------------------------------------------

test_that("an interior empty occasion counts only the identified product", {
  data <- r2_lta_data(n = 150L)
  data[data$time == 2L, c("y1", "y2")] <- NA
  fit <- withCallingHandlers(
    r2_lta_fit(data, transitions = "occasion"),
    latents_unidentified_occasion = function(condition) invokeRestart("muffleWarning"))
  # Three occasion intercepts, of which the move into empty occasion 2 joins
  # the move out of it: K (K - 1) = 2 logits fewer.
  expect_identical(fit$n_parameters, r2_lta_count(p = 2L))
  expect_identical(sum(.lta_identified_coordinates(fit)), fit$n_parameters)
  expect_equal(fit$aic, -2 * fit$log_likelihood + 2 * fit$n_parameters)
  table <- parameter_inference(fit)
  expect_identical(table$estimated, !(table$block == "transition" &
                                        table$term == "occasion_2"))
  skip_if_not_installed("numDeriv")
  expect_identical(r2_lta_hessian_rank(fit), fit$n_parameters)
})

test_that("counted parameters equal identified coordinates across specifications", {
  data <- r2_lta_data()
  data$moving <- data$w * 3 + 5
  specifications <- list(
    list(transitions = "occasion", transition_covariates = "moving"),
    list(order = 2L, transitions = "occasion"),
    list(mover_stayer = TRUE, initial_covariates = "z"),
    list(measurement = "occasion"))
  invisible(lapply(specifications, function(specification) {
    fit <- do.call(r2_lta_fit, c(list(data), specification))
    expect_identical(sum(.lta_identified_coordinates(fit)), fit$n_parameters)
    expect_true(all(.lta_identified_coordinates(fit)))
  }))
})

# Transition grid entropy -------------------------------------------------------

test_that("a weighted transition grid reports the sampling-weighted entropy", {
  data <- r2_lta_data(n = 60L)
  data$sampling <- withr::with_seed(3, stats::runif(60L, 0.2, 3))[data$id]
  grid <- quietly(enumerate_classes(data, c("y1", "y2"), "id", n_profiles = 2,
                                    n_group_classes = 1, time = "time",
                                    weights = "sampling", n_starts = 2, seed = 1))
  fit <- grid$fits[[1L]]
  posteriors <- .multilpa_classified_posteriors(fit)
  weights <- unname(fit$sampling_weights)[fit$group_index]
  expected <- .multilpa_relative_entropy(posteriors, weights)
  expect_equal(grid$table$entropy, expected)
  expect_false(isTRUE(all.equal(expected, .multilpa_relative_entropy(posteriors))))
})
