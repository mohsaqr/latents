# Review 2026-10-26, round 2, multilevel area: zero-weight units in the
# cross-level, covariate, centring and criteria paths; the weighted entropy
# charge; a deterministic first start; the covariate search order; replication
# counts. Each block fails on the code before the fix (see the round-2 report,
# validation/review-2026-10-26/round-2/multilevel.md).

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

r2_two_level <- function(n_groups = 24L, n_rows = 6L, seed = 3L) {
  withr::with_seed(seed, {
    group_class <- sample(1:2, n_groups, TRUE)
    shares <- rbind(c(0.8, 0.2), c(0.25, 0.75))
    do.call(rbind, lapply(seq_len(n_groups), function(group) {
      profile <- sample(1:2, n_rows, TRUE, prob = shares[group_class[group], ])
      centre <- c(-1.5, 1.5)[profile]
      data.frame(g = group, y1 = stats::rnorm(n_rows, centre, 1),
                 y2 = stats::rnorm(n_rows, centre * 0.8, 1.2),
                 z = stats::rnorm(n_rows), wg = stats::rnorm(1))
    }))
  })
}

# The data, and the same data with five extra groups of weight zero whose
# values are moved far away: a zero-weight unit must change nothing.
r2_padded <- function() {
  data <- r2_two_level()
  weights <- withr::with_seed(9, sample(1:3, 24L, TRUE))
  data$w <- weights[data$g]
  extra <- subset(data, g <= 5)
  extra$g <- extra$g + 100L
  extra$w <- 0
  extra$y1 <- extra$y1 * 4 + 6
  extra$y2 <- extra$y2 * 3 - 5
  extra$z <- extra$z * 5
  list(data = data, padded = rbind(data, extra))
}

r2_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(condition) {
    if (inherits(condition, c("latents_boundary", "latents_thin_profiles",
                              "latents_small_classes", "latents_unconverged",
                              "latents_extreme_coefficients"))) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("restricted cross-level pseudo likelihood ignores zero-weight groups (CL-1)", {
  sets <- r2_padded()
  fit <- function(data) r2_quiet(multilpa(
    data, c("y1", "y2"), "g", 2, 2, family = "restricted_cross_level",
    weights = "w", n_starts = 3, seed = 1, tol = 1e-12))
  plain <- fit(sets$data)
  padded <- fit(sets$padded)
  expect_equal(padded$log_likelihood, plain$log_likelihood, tolerance = 1e-10)
  # The reported likelihood is sum_j w_j log L_j at the reported estimates.
  data <- sets$data
  index <- match(data$g, unique(data$g))
  weights <- .latents_sampling_weights(data, "w", index, 24L)
  x <- as.matrix(data[c("y1", "y2")])
  row_density <- vapply(1:2, function(k) {
    plain$profile_prevalences[k] *
      stats::dnorm(x[, 1], plain$profile_means[k, 1], sqrt(plain$profile_variances[k, 1])) *
      stats::dnorm(x[, 2], plain$profile_means[k, 2], sqrt(plain$profile_variances[k, 2]))
  }, numeric(nrow(x)))
  group_means <- rowsum(x, index) / as.vector(table(index))
  group_density <- vapply(1:2, function(h) {
    plain$group_probabilities[h] *
      stats::dnorm(group_means[, 1], plain$group_means[h, 1],
                   sqrt(plain$group_mean_variances[h, 1])) *
      stats::dnorm(group_means[, 2], plain$group_means[h, 2],
                   sqrt(plain$group_mean_variances[h, 2]))
  }, numeric(24L))
  direct <- sum(weights * (log(rowSums(group_density)) +
                             as.vector(rowsum(log(rowSums(row_density)), index))))
  expect_equal(plain$log_likelihood, direct, tolerance = 1e-10)
})

test_that("the equal between-group variance divides by the weighted total (CL-2)", {
  sets <- r2_padded()
  fit <- function(data) r2_quiet(multilpa(
    data, c("y1", "y2"), "g", 2, 2, family = "full_cross_level",
    between_variance = "equal", weights = "w", n_starts = 3, seed = 1, tol = 1e-12))
  plain <- fit(sets$data)
  padded <- fit(sets$padded)
  expect_equal(padded$group_mean_variances, plain$group_mean_variances, tolerance = 1e-8)
  expect_equal(padded$log_likelihood, plain$log_likelihood, tolerance = 1e-10)
})

test_that("cross-level criteria and nobs count units of positive weight (CL-3)", {
  sets <- r2_padded()
  fit <- function(data) r2_quiet(multilpa(
    data, c("y1", "y2"), "g", 2, 2, family = "full_cross_level",
    weights = "w", n_starts = 3, seed = 1, tol = 1e-12))
  plain <- fit(sets$data)
  padded <- fit(sets$padded)
  expect_identical(nobs(padded), 24L)
  expect_identical(attr(logLik(padded), "nobs"), 24L)
  expect_equal(padded$bic, plain$bic, tolerance = 1e-10)
  expect_equal(padded$bic_individual, plain$bic_individual, tolerance = 1e-10)
})

test_that("a zero-weight group with a degenerate person scale does not refuse a covariate fit (ZW-3)", {
  data <- r2_two_level(30L, 5L, 4L)
  data$w <- withr::with_seed(2, sample(1:3, 30L, TRUE))[data$g]
  data$w[data$g == 7] <- 0
  data$y1[data$g == 7] <- 0.5
  # Person-centred data carry no group-level differences for group classes to
  # hold, so one group class keeps the information regular.
  arguments <- list(vars = c("y1", "y2"), id = "g", n_profiles = 2,
                    n_group_classes = 1, weights = "w", n_starts = 2, seed = 1,
                    centering = "person", scaling = "person",
                    profile_covariates = "z")
  fit <- r2_quiet(do.call(multilpa, c(list(data = data), arguments)))
  dropped <- r2_quiet(do.call(multilpa, c(list(data = subset(data, w > 0)), arguments)))
  expect_equal(fit$log_likelihood, dropped$log_likelihood, tolerance = 1e-8)
  # Inference rebuilds the designs with the fit's weights, so it no longer
  # refuses the data (`latents_bad_data`). A person-centred fit defaults to
  # the robust sandwich (INF2-2, round 2), which the zero-weight group must
  # leave unchanged.
  with_zero <- parameter_inference(fit)
  without <- parameter_inference(dropped)
  expect_identical(attr(with_zero, "vcov_type"), "robust")
  expect_equal(with_zero$estimate, without$estimate, tolerance = 1e-6)
  expect_equal(with_zero$standard_error, without$standard_error, tolerance = 1e-6)
})

test_that("grand centring and scaling ignore zero-weight units (ZW-4)", {
  sets <- r2_padded()
  fit <- function(data, model, covariates = character()) r2_quiet(multilpa(
    data, c("y1", "y2"), "g", 2, 2, model = model, weights = "w", n_starts = 3,
    seed = 5, tol = 1e-12, centering = "grand", scaling = "grand",
    profile_covariates = covariates))
  for (model in c("EII", "VVI")) {
    plain <- fit(sets$data, model)
    padded <- fit(sets$padded, model)
    expect_equal(padded$centering_offsets, plain$centering_offsets, tolerance = 1e-12)
    expect_equal(padded$scaling_factors, plain$scaling_factors, tolerance = 1e-12)
    expect_equal(padded$means, plain$means, tolerance = 1e-8)
    expect_equal(padded$log_likelihood, plain$log_likelihood, tolerance = 1e-10)
  }
  plain <- fit(sets$data, "VVI", "z")
  padded <- fit(sets$padded, "VVI", "z")
  expect_equal(padded$scaling_factors, plain$scaling_factors, tolerance = 1e-12)
  expect_equal(padded$log_likelihood, plain$log_likelihood, tolerance = 1e-10)
})

test_that("drop_empty = FALSE still leaves zero-weight units out of every count (ZW-5)", {
  sets <- r2_padded()
  for (covariates in list(character(), "z")) {
    fit <- function(data) r2_quiet(multilpa(
      data, c("y1", "y2"), "g", 2, 2, weights = "w", n_starts = 3, seed = 1,
      tol = 1e-10, drop_empty = FALSE, profile_covariates = covariates))
    plain <- fit(sets$data)
    padded <- fit(sets$padded)
    expect_identical(nobs(padded), nobs(plain))
    expect_equal(padded$bic, plain$bic, tolerance = 1e-10)
    criteria_plain <- get_results(plain, "information_criteria")
    criteria_padded <- get_results(padded, "information_criteria")
    expect_equal(criteria_padded$value, criteria_plain$value, tolerance = 1e-8)
    expect_equal(criteria_padded$n, criteria_plain$n)
  }
})

test_that("a covariate fit's bic_individual is its criteria table's under drop_empty = FALSE (CV-1)", {
  data <- r2_two_level()
  data[c(5, 40, 41), c("y1", "y2")] <- NA
  fit <- r2_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2, profile_covariates = "z",
                           missing = "fiml", n_starts = 2, seed = 1,
                           drop_empty = FALSE))
  criteria <- .multilpa_information_criteria(fit, format = "long")
  individual <- subset(criteria, criterion == "bic" & convention == "individuals")
  expect_identical(individual$n, 141L)
  expect_equal(fit$bic_individual, individual$value, tolerance = 1e-12)
})

test_that("weighted ICL charges the weighted entropy, that of the replicated data (EN-1)", {
  data <- r2_two_level(30L, 5L, 6L)
  raw <- withr::with_seed(3, sample(1:4, 30L, TRUE))
  data$w <- raw[data$g]
  replicated <- do.call(rbind, lapply(seq_len(30L), function(group) {
    block <- subset(data, g == group)
    do.call(rbind, lapply(seq_len(raw[group]), function(copy) {
      block$g <- group * 100L + copy
      block
    }))
  }))
  weighted <- r2_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2, weights = "w",
                                n_starts = 3, seed = 2, tol = 1e-13, max_iter = 20000))
  copies <- r2_quiet(multilpa(replicated, c("y1", "y2"), "g", 2, 2, n_starts = 3,
                              seed = 2, tol = 1e-13, max_iter = 20000))
  scale <- 30 / sum(raw)
  expect_equal(copies$log_likelihood * scale, weighted$log_likelihood, tolerance = 1e-7)
  criteria <- .multilpa_information_criteria(weighted, format = "long")
  q <- weighted$n_parameters
  # Both conventions charge the complete-data entropy, group classes and
  # profiles together (selection round 2, S-1).
  group_weights <- unname(weighted$sampling_weights)
  group_entropy <- -rowSums(ifelse(weighted$group_posteriors > 0,
                                   weighted$group_posteriors *
                                     log(weighted$group_posteriors), 0))
  expected <- sum(group_weights * group_entropy) +
    sum(group_weights[weighted$group_index] * weighted$classification_entropy_rows)
  for (level in c("individuals", "groups")) {
    row <- subset(criteria, criterion == "icl" & convention == level)
    charge <- (row$value + 2 * weighted$log_likelihood - q * log(row$n)) / 2
    # Exactly the weighted sum of the fit's own per-unit entropies ...
    expect_equal(charge, expected, tolerance = 1e-10)
    # ... which is the replicated data's entropy on the weights' scale, to the
    # precision at which the two EM runs stop on a flat group-class surface.
    expect_equal(charge, .multilpa_complete_entropy(copies) * scale, tolerance = 1e-4)
  }
})

test_that("start 1 is the same whatever the seed, with several group classes (ST-1)", {
  data <- r2_two_level(40L, 6L, 12L)
  scores <- vapply(1:4, function(seed) {
    r2_quiet(multilpa(data, c("y1", "y2"), "g", 3, 3, n_starts = 1, max_iter = 0,
                      seed = seed, search = "kmeans"))$log_likelihood
  }, numeric(1))
  expect_identical(length(unique(scores)), 1L)
})

test_that("a covariate emEM pool begins with the kmeans search's starts (ST-2)", {
  data <- r2_two_level(40L, 5L, 14L)
  fit <- function(search, n_screen = NULL) r2_quiet(multilpa(
    data, c("y1", "y2"), "g", 2, 2, model = "VVI", profile_covariates = "z",
    search = search, n_screen = n_screen, n_starts = 4, seed = 3))
  kmeans <- get_results(fit("kmeans"), "starts")
  emem <- get_results(fit("emEM", 4), "starts")
  expect_identical(emem$log_likelihood, kmeans$log_likelihood)
  expect_identical(emem$iterations, kmeans$iterations)
})

test_that("covariate and cross-level fits count replicated maxima as multilpa() does (RP-1)", {
  data <- r2_two_level()
  recount <- function(fit, starts) {
    final <- if (is.null(starts$stage)) TRUE else starts$stage == "final"
    sum(final & is.finite(starts$log_likelihood) &
          abs(starts$log_likelihood - fit$log_likelihood) <=
          1e-6 * (1 + abs(fit$log_likelihood)))
  }
  covariate <- r2_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                 profile_covariates = "z", n_starts = 4, seed = 1))
  expect_true(is.numeric(covariate$n_best_replicated))
  expect_identical(covariate$n_best_replicated,
                   recount(covariate, get_results(covariate, "starts")))
  full <- r2_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2, family = "full_cross_level",
                            n_starts = 4, seed = 1))
  expect_identical(full$n_best_replicated,
                   recount(full, subset(full$starts, !is.na(log_likelihood))))
  restricted <- r2_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                  family = "restricted_cross_level",
                                  n_starts = 4, seed = 1))
  expect_identical(restricted$n_best_replicated,
                   min(restricted$starts$n_best_replicated))
})

test_that("the weights line counts the units of positive weight (DOC-1)", {
  data <- r2_two_level(30L, 5L, 4L)
  data$w <- withr::with_seed(2, sample(1:3, 30L, TRUE))[data$g]
  data$w[data$g <= 4] <- 0
  fit <- r2_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2, weights = "w",
                           n_starts = 2, seed = 1))
  printed <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "of 26 units (4 more of weight zero)", fixed = TRUE)
})
