test_that("zero-count rows do not change the membership optimum or its score", {
  z <- seq(-2, 2, length.out = 81)
  design <- cbind(1, z)
  probability <- plogis(.3 + 2 * z)
  counts <- cbind(probability, 1 - probability)
  initial <- matrix(0, 2, 1)
  ordinary <- .multilpa_weighted_logits(design, counts, initial)
  padded <- .multilpa_weighted_logits(rbind(design, c(1, 1e12)),
                                      rbind(counts, c(0, 0)), initial)
  expect_true(padded$converged)
  expect_equal(drop(padded$coefficients), c(.3, 2), tolerance = 1e-7)
  expect_equal(padded$coefficients, ordinary$coefficients, tolerance = 1e-10)
  expect_equal(padded$scaled_score, ordinary$scaled_score, tolerance = 1e-12)
  empty <- .multilpa_weighted_logits(matrix(c(1, 1e12), 1, 2),
                                     matrix(0, 1, 2), initial)
  expect_identical(empty$coefficients, initial)
})

test_that("zero-weight predictor outliers do not change a public covariate fit", {
  set.seed(41)
  z <- rnorm(600)
  profile <- rbinom(600, 1, plogis(.4 + 2 * z))
  d <- data.frame(y = 4 * profile + rnorm(600, 0, .4), z = z, w = 1)
  fit <- function(value) lpa(rbind(d, data.frame(y = 2, z = value, w = 0)),
    "y", 2, profile_covariates = "z", weights = "w", variance_model = "equal",
    n_starts = 1, search = "kmeans", seed = 2, max_iter = 1000, tol = 1e-10)
  original <- fit(0)
  outlier <- fit(1e12)
  expect_true(outlier$converged)
  expect_false(outlier$extreme_logits)
  expect_equal(outlier$log_likelihood, original$log_likelihood, tolerance = 1e-8)
  expect_equal(outlier$profile_coefficients, original$profile_coefficients,
               tolerance = 1e-7)
  expect_equal(outlier$subject_posteriors[seq_len(600), ],
               original$subject_posteriors[seq_len(600), ], tolerance = 1e-7)
  expect_equal(parameter_inference(outlier)$standard_error,
               parameter_inference(original)$standard_error, tolerance = 1e-7)
})

test_that("an extreme common count density preserves priors and other indicators", {
  d <- data.frame(x = c(rep(c(-10.1, -9.9), 40), rep(c(9.9, 10.1), 10)),
                  y = rep(c(1, 2), 50))
  fit <- lpa(d, c("x", "y"), 2, count = "y", variance_model = "equal",
              n_starts = 1, seed = 1, tol = 1e-12, missing = "fiml")
  p <- predict(fit, newdata = data.frame(x = c(NA_real_, NA_real_), y = c(1, 1e16)))
  expect_equal(p$posterior, rep(max(fit$profile_probabilities), 2),
               tolerance = 1e-12)
  # A second count still carries information beside the huge common density.
  extra <- list(count = matrix(c(1e16, 2), 1, 2), ordinal = NULL)
  parameters <- list(means = matrix(numeric(), 2, 0),
    variances = matrix(numeric(), 2, 0), count_means = cbind(c(1.5, 1.5), c(1, 3)),
    profile_probabilities = matrix(c(.8, .2), 1), group_probabilities = 1)
  result <- .multilpa_expectation(matrix(numeric(), 1, 0), 1L, parameters,
                                   extra = extra)
  expected <- c(.8, .2) * dpois(2, c(1, 3))
  expect_equal(drop(result$subject_posteriors), expected / sum(expected),
               tolerance = 1e-12)
  expect_equal(result$log_likelihood,
    dpois(1e16, 1.5, log = TRUE) + log(sum(expected)), tolerance = 1e-14)
})

test_that("ordinary categorical measurements retain their original arithmetic", {
  codes <- matrix(c(1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L), 4, 2)
  probabilities <- list(matrix(c(.13, .87, .71, .29), 2, byrow = TRUE),
                        matrix(c(.42, .58, .83, .17), 2, byrow = TRUE))
  parameters <- list(means = matrix(numeric(), 2, 0),
    variances = matrix(numeric(), 2, 0), response_probabilities = probabilities)
  expected <- matrix(0, 4, 2)
  for (j in seq_len(ncol(codes))) {
    expected <- expected + t(log(probabilities[[j]][, codes[, j], drop = FALSE]))
  }
  x <- matrix(numeric(), 4, 0)
  for (extra in list(NULL, list(), list(ordinal = matrix(integer(), 4, 0),
                                       count = matrix(numeric(), 4, 0)))) {
    density <- .latents_measurement_log_density(x, parameters, codes, extra)
    expect_identical(density$log_density, expected)
    expect_identical(density$offset, rep(0, 4))
  }
})

test_that("ordinal whole numbers retain distinct native categories", {
  values <- c(1, 1e15, 1e15 + 1)
  d <- data.frame(y = rep(values, c(20, 30, 50)))
  fit <- lca(d, "y", 1, ordinal = "y", n_starts = 1, seed = 1)
  expect_identical(fit$extra_data$ordinal_categories, 3L)
  expect_identical(anyDuplicated(fit$extra_data$ordinal_levels$y), 0L)
  expect_identical(drop(fit$extra_data$ordinal), rep(1:3, c(20, 30, 50)))
  expect_equal(fit$log_likelihood, sum(c(20, 30, 50) * log(c(.2, .3, .5))),
               tolerance = 1e-10)
  density <- predict(fit, data.frame(y = values), type = "density")
  expect_equal(density$log_density, log(c(.2, .3, .5)), tolerance = 1e-9)
  # Simulation/refit factors use the same unambiguous labels.
  restored <- .latents_draw_extra_columns(fit$extra_data)
  prepared <- .latents_prepare_extra_like(fit, as.data.frame(restored))
  expect_identical(prepared$ordinal, fit$extra_data$ordinal)
})

test_that("indicator scaling depends on variation rather than its origin", {
  set.seed(101)
  d <- data.frame(y = rnorm(100))
  fit <- function(data) lpa(data, "y", 1, centering = "grand", scaling = "grand",
                             n_starts = 1, seed = 1)
  original <- fit(d)
  shifted_data <- transform(d, y = y + 1e9)
  shifted <- fit(shifted_data)
  expect_equal(unname(shifted$scaling_factors), matrix(sd(shifted_data$y), 1, 1),
               tolerance = 1e-7)
  expect_equal(shifted$variances, original$variances, tolerance = 1e-10)
  expect_equal(shifted$log_likelihood, original$log_likelihood, tolerance = 1e-8)
  expect_error(fit(data.frame(y = rep(1e9, 100))), class = "latents_bad_data")
  within <- matrix(c(1e9 + c(-1, 0, 1), -1e9 + c(-2, 0, 2)), 6, 1)
  expect_equal(drop(.multilpa_indicator_scales(within, rep(1:2, each = 3), 2,
                                                "person")), c(1, 2))
})

test_that("within-unit weight validation is invariant to a common scale", {
  group <- rep(1:2, each = 2)
  invisible(lapply(c(1e-300, 1e-12, 1, 1e300), function(scale) {
    expect_error(.latents_sampling_weights(data.frame(w = scale * c(1, 2, 3, 3)),
      "w", group, 2), class = "latents_bad_weights")
    expect_equal(.latents_sampling_weights(data.frame(w = scale * c(1, 1, 3, 3)),
      "w", group, 2), c(.5, 1.5), tolerance = 1e-14)
    expect_error(.latents_sampling_weights(data.frame(w = scale * c(0, 1, 3, 3)),
      "w", group, 2), class = "latents_bad_weights")
  }))
})

test_that("mixture structure posteriors preserve priors at extreme common densities", {
  prior <- c(.8, .2)
  density <- matrix(-1e16, 2, 2)
  observation <- list(nesting = "observation", w = matrix(1, 2, 1))
  parameters <- list(gamma = matrix(log(prior / prior[1]), 1))
  expected <- matrix(rep(prior, each = 2), 2)
  value <- .mixture_structure_expectation(observation, parameters, density)
  expect_equal(unname(value$tau), expected, tolerance = 1e-14)
  grouped <- observation
  grouped$nesting <- "group"
  grouped$w <- matrix(1, 1, 1)
  grouped$group_index <- c(1L, 1L)
  value <- .mixture_structure_expectation(grouped, parameters, density)
  expect_equal(drop(value$group_tau), prior, tolerance = 1e-14)
  expect_equal(unname(value$tau), expected, tolerance = 1e-14)
  # Conditional profile shares differ across group classes; equal measurements
  # leave both levels' priors intact.
  nested <- list(nesting = "two-level", n = 2L, n_groups = 1L,
    n_group_classes = 2L, group_index = c(1L, 1L),
    w = matrix(numeric(), 2, 0), v = matrix(1, 1, 1))
  shares <- rbind(c(.8, .2), c(.3, .7))
  parameters <- list(class_logits = cbind(0, log(shares[, 2] / shares[, 1])),
                      delta = matrix(c(0, log(.4 / .6)), 1))
  value <- .mixture_structure_expectation(nested, parameters, density)
  expect_equal(drop(value$rho), c(.6, .4), tolerance = 1e-14)
  expect_equal(value$tau, matrix(rep(drop(c(.6, .4) %*% shares), each = 2), 2),
               tolerance = 1e-14)
  expect_equal(rowSums(value$tau), c(1, 1), tolerance = 1e-14)
})

test_that("only positive-weight rows can identify membership predictors", {
  set.seed(11)
  d <- data.frame(y = rnorm(81), z = c(rep(0, 80), 1), w = c(rep(1, 80), 0))
  expect_error(lpa(d, "y", 2, profile_covariates = "z", weights = "w",
                    n_starts = 1, seed = 1), class = "latents_bad_covariate")
})

test_that("covariate logLik and nobs use the BIC group convention", {
  set.seed(13)
  d <- data.frame(g = rep(1:20, each = 4),
                  y = rep(c(-3, 3), each = 40) + rnorm(80, sd = .5), z = rnorm(80))
  padded <- rbind(d, data.frame(g = 21:30, y = NA_real_, z = 0))
  fit <- function(drop_empty) multilpa(padded, "y", "g", 2, 1,
    profile_covariates = "z", missing = "fiml", drop_empty = drop_empty,
    n_starts = 1, seed = 1)
  observed <- fit(TRUE)
  all <- fit(FALSE)
  expect_equal(nobs(observed), 20)
  expect_equal(attr(logLik(observed), "nobs"), 20)
  expect_equal(stats::BIC(observed), observed$bic)
  expect_equal(nobs(all), 30)
  expect_equal(stats::BIC(all), all$bic)
})

test_that("the additive study CLI validates its stage before loading helpers", {
  path <- testthat::test_path("..", "..", "validation", "additive-simulation.R")
  skip_if_not(file.exists(path), "supporting study is outside the installed package")
  expressions <- parse(path)
  first_source <- which(vapply(expressions, function(x) {
    is.call(x) && identical(x[[1L]], as.name("source"))
  }, logical(1)))[1L]
  resolve <- function(arguments) {
    environment <- new.env(parent = baseenv())
    environment$commandArgs <- function(trailingOnly) arguments
    eval(expressions[seq_len(first_source - 1L)], environment)
    c(stage = environment$stage, datasets = as.character(environment$n_datasets))
  }
  expect_identical(resolve(character()), c(stage = "pilot", datasets = "100"))
  expect_identical(resolve("full"), c(stage = "full", datasets = "1000"))
  expect_error(resolve("typo"), "stage")
  expect_error(resolve(c("pilot", "full")), "at most one")
})

test_that("the recovery pivot requires a complete unique grid of parameter cells", {
  path <- testthat::test_path("..", "..", "validation", "simulab-recovery",
                              "recovery-study.R")
  skip_if_not(file.exists(path), "supporting study is outside the installed package")
  expressions <- parse(path)
  definition <- Filter(function(x) {
    is.call(x) && identical(x[[1L]], as.name("<-")) &&
      identical(x[[2L]], as.name(".as_parameter_matrix"))
  }, as.list(expressions))[[1L]]
  pivot <- eval(definition[[3L]], asNamespace("stats"))
  table <- data.frame(row = c(2, 1, 2, 1), column = c("b", "a", "a", "b"),
                       value = c(0, 10, 12, 11))
  expected <- matrix(c(10, 12, 11, 0), 2, dimnames = list(NULL, c("a", "b")))
  expect_identical(pivot(table, "value", "row", "column", c("a", "b")), expected)
  expect_error(pivot(table[-1, ], "value", "row", "column", c("a", "b")),
                 "exactly one")
  expect_error(pivot(rbind(table, table[1, ]), "value", "row", "column", c("a", "b")),
                 "exactly one")
  expect_error(pivot(transform(table, value = NA_real_), "value", "row", "column",
                       c("a", "b")), "finite")
})

test_that("Poisson class evidence survives subtraction of huge log pmfs", {
  means <- c(1.5, 1.5 + 2^-50)
  extra <- list(count = matrix(1e16, 1, 1), ordinal = NULL)
  parameters <- list(means = matrix(numeric(), 2, 0),
    variances = matrix(numeric(), 2, 0), count_means = matrix(means, 2, 1),
    profile_probabilities = matrix(c(.8, .2), 1), group_probabilities = 1)
  value <- .multilpa_expectation(matrix(numeric(), 1, 0), 1L, parameters,
                                  extra = extra)
  log_ratio <- 1e16 * log1p((means[2] - means[1]) / means[1]) -
    (means[2] - means[1])
  expected <- c(.8, .2 * exp(log_ratio))
  expect_equal(drop(value$subject_posteriors), expected / sum(expected),
               tolerance = 1e-12)
})
