# Cross-model integration regressions for the 0.9.0--0.9.7 review.
review_discrete_rows <- function() {
  set.seed(712)
  n <- 480L
  z <- rep(1:2, each = n / 2)
  data.frame(g = rep(seq_len(n / 3), each = 3), x = rnorm(n),
             o = ordered(vapply(z, function(k) {
               sample(1:4, 1L, prob = if (k == 1L) c(.7, .2, .08, .02) else
                 c(.02, .08, .2, .7))
             }, integer(1))),
             k = rnbinom(n, mu = c(2, 14)[z], size = 2),
             a = factor(sample(c("a", "b"), n, TRUE)))
}

review_discrete_fit <- function(data, profiles = 2L, ...) {
  quietly(multilpa(data, c("o", "k"), "g", profiles, n_group_classes = 1,
                   ordinal = "o", count = "k", count_model = "negative_binomial",
                   n_starts = 2, seed = 7, tol = 1e-10, ...))
}

test_that("class wrappers partition structured discrete indicators", {
  data <- review_discrete_rows()
  arguments <- list(data = data, vars = c("o", "k", "a"), n_profiles = 2,
                    ordinal = "o", count = "k", count_model = "negative_binomial",
                    n_starts = 2, seed = 7)
  reference <- quietly(do.call(multilpa, c(arguments, list(id = NULL, categorical = "a"))))
  single <- do.call(lca, c(arguments[setdiff(names(arguments), "n_profiles")],
                          list(n_classes = 2)))
  grouped <- do.call(multilca, c(arguments, list(id = "g", n_group_classes = 1)))
  expect_identical(single$categorical, "a")
  expect_identical(grouped$categorical, "a")
  expect_equal(single$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  grid <- enumerate_lca(data, c("o", "k", "a"), n_classes = 1,
                        ordinal = "o", count = "k", n_starts = 1, seed = 7)
  expect_true(is.finite(grid$table$log_likelihood))
  expect_identical(candidate_fit(grid, 1)$categorical, "a")
})

test_that("profile permutation preserves ordinal and NB likelihood and labels", {
  data <- review_discrete_rows()
  fit <- review_discrete_fit(data)
  swapped <- .multilpa_permute_profiles(fit, 2:1)
  expect_equal(predict(swapped, data, type = "density")$log_density,
               predict(fit, data, type = "density")$log_density, tolerance = 1e-12)
  expect_equal(swapped$ordinal_locations[2L, ], 0)
  expect_equal(swapped$count_means, fit$count_means[2:1, , drop = FALSE], ignore_attr = TRUE)
  expect_equal(swapped$count_dispersion, fit$count_dispersion[2:1, , drop = FALSE],
               ignore_attr = TRUE)
  # Force an uninformative mixing signature: measurements must decide matching.
  fit$profile_probabilities[] <- 0.5
  swapped$profile_probabilities[] <- 0.5
  aligned <- .multilpa_align_labels(swapped, fit)
  expect_equal(aligned$ordinal_intercepts, fit$ordinal_intercepts, tolerance = 1e-12)
  expect_equal(aligned$ordinal_locations, fit$ordinal_locations, tolerance = 1e-12)
  expect_equal(aligned$count_means, fit$count_means)
  expect_equal(aligned$count_dispersion, fit$count_dispersion)
  expect_equal(aligned$subject_profiles, fit$subject_profiles)
})

test_that("covariate relabelling preserves extra measurement and likelihood", {
  data <- review_discrete_rows()
  fit <- review_discrete_fit(data, profile_covariates = "x")
  swapped <- .multilpa_cov_permute(fit, 2:1, 1L)
  expect_silent(.multilpa_cov_check_likelihood(swapped, data))
  expect_equal(swapped$count_means, fit$count_means[2:1, , drop = FALSE], ignore_attr = TRUE)
  expect_equal(swapped$ordinal_locations[2L, ], 0)
  expect_identical(.multilpa_match_order(.multilpa_cov_profile_signature(fit),
                                        .multilpa_cov_profile_signature(swapped, fit)), 2:1)
})

test_that("BLRT reproduces NB likelihood and refuses different indicator specifications", {
  data <- review_discrete_rows()
  null <- review_discrete_fit(data, 1L)
  alternative <- review_discrete_fit(data)
  result <- quietly(bootstrap_lrt(null, alternative, data, iter = 2,
                                  n_starts = 1, seed = 7))
  expect_true(is.finite(result$statistic))
  expect_equal(nrow(result$replicates), 2L)
  changed <- alternative
  changed$extra_data$count_model <- "poisson"
  expect_error(bootstrap_lrt(null, changed, data, iter = 2),
               class = "latents_incomparable_models")
  changed <- alternative
  changed$extra_data$count_dispersion <- "equal"
  expect_error(bootstrap_lrt(null, changed, data, iter = 2),
               class = "latents_incomparable_models")
  changed <- alternative
  changed$ordinal <- character()
  expect_error(bootstrap_lrt(null, changed, data, iter = 2),
               class = "latents_incomparable_models")
})

test_that("row alignment checks ordinal and count records within groups", {
  data <- review_discrete_rows()
  fit <- review_discrete_fit(data)
  expect_silent(.multilpa_check_alignment(fit, data))
  changed <- data
  changed$k[1L] <- changed$k[1L] + 1L
  expect_error(.multilpa_check_alignment(fit, changed), class = "latents_bad_inference_data")
  changed <- data
  changed$o[1L] <- if (changed$o[1L] == levels(changed$o)[1L]) levels(changed$o)[4L] else
    levels(changed$o)[1L]
  expect_error(.multilpa_check_alignment(fit, changed), class = "latents_bad_inference_data")
})

test_that("finite sampling weights remain finite at extreme common scales", {
  data <- data.frame(w = c(1, 2, 3))
  reference <- .latents_sampling_weights(data, "w", 1:3, 3L)
  data$w <- data$w * 1e307
  actual <- .latents_sampling_weights(data, "w", 1:3, 3L)
  expect_equal(actual, reference, tolerance = 1e-14)
  data$w <- rep(1e308, 3)
  expect_equal(.latents_sampling_weights(data, "w", 1:3, 3L), rep(1, 3))
  data$w <- rep(1e-310, 3)
  expect_equal(.latents_sampling_weights(data, "w", 1:3, 3L), rep(1, 3))
})

test_that("three-profile ordinal relabelling rebases every item exactly", {
  data <- data.frame(o = ordered(rep(1:4, 2)), p = ordered(rep(1:2, each = 4)),
                     k = 0:7)
  extra <- .latents_set_count_model(.latents_prepare_extra(data, c("o", "p"), "k", "error"),
                                     "negative_binomial", "equal")
  parameters <- list(n_profiles = 3L, ordinal = c("o", "p"),
                     ordinal_intercepts = list(o = c(.4, -.2, .7), p = .3),
                     ordinal_locations = rbind(c(-.8, .2), c(.5, -.6), c(0, 0)),
                     count_means = matrix(c(1, 4, 7), 3L),
                     count_dispersion = matrix(rep(.4, 3L), 3L))
  order <- c(2L, 3L, 1L)
  original <- .latents_extra_log_density(extra, parameters, 8L, 3L)
  reordered <- .multilpa_permute_extra(parameters, order)
  actual <- .latents_extra_log_density(extra, reordered, 8L, 3L)
  expect_equal(actual, original[, order], tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(unname(reordered$ordinal_locations[3L, ]), c(0, 0))
  restored <- .multilpa_permute_extra(reordered, match(1:3, order))
  expect_equal(restored$ordinal_intercepts, parameters$ordinal_intercepts,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(restored$ordinal_locations, parameters$ordinal_locations,
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("single-profile NB varying dispersion labels retain its standard error", {
  set.seed(431)
  data <- data.frame(k = rnbinom(350, mu = 4, size = 2))
  lapply(c("equal", "varying"), function(mode) {
    fit <- lpa(data, "k", 1, count = "k", count_model = "negative_binomial",
                count_dispersion = mode, n_starts = 1)
    table <- get_results(fit, "count_means")
    inference <- subset(parameter_inference(fit), parameter == "count_dispersion")
    expect_true(is.finite(table$dispersion_standard_error))
    expect_equal(table$dispersion_standard_error, inference$standard_error)
    expect_identical(inference$outcome, if (mode == "equal") "shared" else "profile_1")
  })
})
