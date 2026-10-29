# Properties of direct effects of covariates on continuous indicators
# (`direct_effects = "y3 ~ z"`, y = mean_k + b z + e): invariances, unit
# changes, weights, prediction, the refit paths and the error classes. The
# comparisons with Latent GOLD, Mplus and an independent implementation are in
# equivalence/test-features-direct-effects.R. Data are drawn in base R here
# (simulab is not a declared dependency of the package's tests).
# The comparisons with independent references are in tests/equivalence/test-features-direct-effects.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

de_unit_data <- function(seed = 7, n = 400, b = c(0.8, -0.6)) {
  withr::with_seed(seed, {
    class <- sample.int(2L, n, replace = TRUE, prob = c(0.55, 0.45))
    z <- round(stats::rnorm(n, 0.5, 1.5), 4)
    w <- round(stats::rnorm(n, -1, 1), 4)
    y <- vapply(1:4, \(j) stats::rnorm(n, c(0, c(1.8, 1.5, 1.6, 1.2)[[j]])[class],
                                       c(1, c(0.8, 1.2, 0.9, 1.1)[[j]])[class]), numeric(n))
    data <- data.frame(id = seq_len(n), z = z, w = w, round(y, 4))
    names(data)[4:7] <- paste0("y", 1:4)
    data$y3 <- data$y3 + b[[1L]] * data$z
    data$y1 <- data$y1 + b[[2L]] * data$z
    data
  })
}

de_unit_fit <- function(data, effects = "y3 ~ z", model = "VVI", n_starts = 5, ...) {
  quietly(lpa(data, paste0("y", 1:4), 2, model = model, n_starts = n_starts, seed = 1,
              tol = 1e-12, max_iter = 5000, direct_effects = effects, ...))
}

# Classes of `other` in the order of `fit`'s, from the posteriors of the same rows.
de_unit_order <- function(fit, other) {
  if (mean(abs(fit$subject_posteriors[, 1] - other$subject_posteriors[, 1])) < 0.5) 1:2 else 2:1
}

test_that("permuting the rows leaves the fit unchanged", {
  data <- de_unit_data()
  fit <- de_unit_fit(data, c("y3 ~ z", "y1 ~ z"))
  permutation <- withr::with_seed(11, sample.int(nrow(data)))
  permuted <- de_unit_fit(data[permutation, ], c("y3 ~ z", "y1 ~ z"))
  expect_equal(permuted$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  expect_equal(permuted$direct_effects, fit$direct_effects, tolerance = 1e-7)
  order <- de_unit_order(fit, list(subject_posteriors =
                                     permuted$subject_posteriors[order(permutation), ]))
  expect_equal(unname(permuted$means[order, ]), unname(fit$means), tolerance = 1e-7)
  expect_equal(unname(permuted$variances[order, ]), unname(fit$variances), tolerance = 1e-7)
  expect_equal(unname(permuted$subject_posteriors[order(permutation), order]),
               unname(fit$subject_posteriors), tolerance = 1e-7)
})

test_that("column order, term order and unused columns do not change the fit", {
  data <- de_unit_data()
  fit <- de_unit_fit(data, "y3 ~ z + w")
  shuffled <- data[c("w", "y4", "y2", "z", "y1", "id", "y3")]
  shuffled$unused <- rev(seq_len(nrow(data)))
  other <- de_unit_fit(shuffled, "y3 ~ w + z")
  expect_equal(other$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  expect_equal(other$direct_effects[c("y3:z", "y3:w")], fit$direct_effects[c("y3:z", "y3:w")],
               tolerance = 1e-7)
  expect_equal(other$n_parameters, fit$n_parameters)
  # Several indicators given in another order are the same model.
  both <- de_unit_fit(data, c("y3 ~ z", "y1 ~ z"))
  swapped <- de_unit_fit(data, c("y1 ~ z", "y3 ~ z"))
  expect_equal(swapped$log_likelihood, both$log_likelihood, tolerance = 1e-10)
  expect_equal(swapped$direct_effects[names(both$direct_effects)], both$direct_effects,
               tolerance = 1e-7)
})

test_that("relabelling the groups of a two-level fit does not change it", {
  skip_on_cran()
  data <- transform(de_unit_data(n = 300), g = rep(seq_len(30), each = 10))
  relabelled <- transform(data, g = sprintf("school-%03d", 100L - g))
  vars <- paste0("y", 1:4)
  first <- quietly(multilpa(data, vars, "g", 2, 1, n_starts = 3, seed = 1, tol = 1e-12,
                            direct_effects = "y3 ~ z"))
  second <- quietly(multilpa(relabelled, vars, "g", 2, 1, n_starts = 3, seed = 1,
                             tol = 1e-12, direct_effects = "y3 ~ z"))
  expect_equal(second$log_likelihood, first$log_likelihood, tolerance = 1e-10)
  expect_equal(second$direct_effects, first$direct_effects, tolerance = 1e-7)
})

test_that("shifting a covariate moves the means by -b c and nothing else", {
  data <- de_unit_data()
  effects <- c("y3 ~ z", "y1 ~ z")
  base <- de_unit_fit(data, effects)
  shift <- -3.7
  moved <- de_unit_fit(transform(data, z = z + shift), effects)
  order <- de_unit_order(base, moved)
  expect_equal(moved$log_likelihood, base$log_likelihood, tolerance = 1e-10)
  expect_equal(moved$direct_effects, base$direct_effects, tolerance = 1e-7)
  expect_equal(unname(moved$variances[order, ]), unname(base$variances), tolerance = 1e-7)
  expected <- base$means
  expected[, "y3"] <- expected[, "y3"] - base$direct_effects[["y3:z"]] * shift
  expected[, "y1"] <- expected[, "y1"] - base$direct_effects[["y1:z"]] * shift
  expect_equal(unname(moved$means[order, ]), unname(expected), tolerance = 1e-7)
  # The effects' standard errors exist and do not depend on where the
  # covariate is centred.
  base_errors <- get_results(base, "direct_effects")$standard_error
  expect_true(all(is.finite(base_errors) & base_errors > 0))
  expect_equal(get_results(moved, "direct_effects")$standard_error, base_errors,
               tolerance = 1e-5)
})

test_that("scaling a covariate scales its effect by the inverse", {
  data <- de_unit_data()
  base <- de_unit_fit(data, "y3 ~ z + w")
  invisible(lapply(c(-2.5, 1e3, 1e-3), \(factor) {
    scaled <- de_unit_fit(transform(data, z = z * factor), "y3 ~ z + w")
    expect_equal(scaled$log_likelihood, base$log_likelihood, tolerance = 1e-10)
    expect_equal(scaled$direct_effects[["y3:z"]] * factor, base$direct_effects[["y3:z"]],
                 tolerance = 1e-7)
    expect_equal(scaled$direct_effects[["y3:w"]], base$direct_effects[["y3:w"]],
                 tolerance = 1e-7)
    expect_equal(get_results(scaled, "direct_effects")$standard_error * c(abs(factor), 1),
                 get_results(base, "direct_effects")$standard_error, tolerance = 1e-5)
  }))
})

test_that("scaling an indicator scales its mean, effect and standard deviation", {
  data <- de_unit_data()
  base <- de_unit_fit(data, "y3 ~ z")
  factor <- 40
  scaled <- de_unit_fit(transform(data, y3 = y3 * factor), "y3 ~ z")
  order <- de_unit_order(base, scaled)
  # The density of y3 changes by the Jacobian, 1 / factor per row.
  expect_equal(scaled$log_likelihood, base$log_likelihood - nrow(data) * log(factor),
               tolerance = 1e-10)
  # EM stops on a relative likelihood change of 1e-12, which pins the
  # likelihood to that order but the parameters only to about its square root
  # (the likelihood is flat to second order at the maximum): 1e-5.
  expect_equal(scaled$direct_effects[["y3:z"]], factor * base$direct_effects[["y3:z"]],
               tolerance = 1e-5)
  expect_equal(unname(scaled$means[order, "y3"]), factor * unname(base$means[, "y3"]),
               tolerance = 1e-5)
  expect_equal(unname(sqrt(scaled$variances[order, "y3"])),
               factor * unname(sqrt(base$variances[, "y3"])), tolerance = 1e-5)
  expect_equal(unname(scaled$means[order, "y1"]), unname(base$means[, "y1"]),
               tolerance = 1e-5)
})

test_that("integer weights fit as the rows repeated", {
  data <- de_unit_data(n = 250)
  data$copies <- withr::with_seed(5, sample(1:3, nrow(data), replace = TRUE))
  weighted <- de_unit_fit(data, "y3 ~ z", weights = "copies")
  repeated <- de_unit_fit(data[rep(seq_len(nrow(data)), data$copies), ], "y3 ~ z")
  order <- de_unit_order(weighted, list(subject_posteriors =
    repeated$subject_posteriors[!duplicated(rep(seq_len(nrow(data)), data$copies)), ]))
  expect_equal(weighted$direct_effects, repeated$direct_effects, tolerance = 1e-7)
  expect_equal(unname(weighted$means[order, ]), unname(repeated$means), tolerance = 1e-7)
  expect_equal(unname(weighted$variances[order, ]), unname(repeated$variances),
               tolerance = 1e-7)
  # The weights are scaled to sum to the number of rows (Mplus's convention).
  expect_equal(weighted$log_likelihood * sum(data$copies) / nrow(data),
               repeated$log_likelihood, tolerance = 1e-9)
})

test_that("predict() on the training data reproduces the posteriors", {
  data <- de_unit_data()
  fit <- de_unit_fit(data, c("y3 ~ z + w", "y1 ~ z"))
  predicted <- predict(fit, newdata = data, type = "posterior")
  expect_equal(unname(as.matrix(predicted[c("probability_profile_1", "probability_profile_2")])),
               unname(fit$subject_posteriors), tolerance = 1e-10)
  expect_identical(predict(fit, newdata = data, type = "class")$profile,
                   unname(fit$subject_profiles))
  # The covariates are read from newdata: moving them moves the posteriors.
  moved <- predict(fit, newdata = transform(data, z = z + 2), type = "posterior")
  expect_gt(max(abs(moved$probability_profile_1 - predicted$probability_profile_1)), 0.1)
  expect_error(predict(fit, newdata = transform(data, w = as.character(w))),
               class = "latents_bad_data")
})

test_that("bootstrap, sensitivity and refits carry the direct effects", {
  skip_on_cran()
  data <- de_unit_data(n = 300)
  fit <- de_unit_fit(data, c("y3 ~ z", "y1 ~ z"), n_starts = 3)
  arguments <- .multilpa_refit_arguments(fit)
  expect_setequal(arguments$direct_effects, c("y3 ~ z", "y1 ~ z"))
  seen <- new.env()
  seen$effects <- list()
  suppressMessages(trace("multilpa", where = asNamespace("latents"), print = FALSE,
                         exit = bquote(assign("effects", c(get("effects", .(seen)),
                                                           list(returnValue()$direct_effects)),
                                              envir = .(seen)))))
  on.exit(suppressMessages(untrace("multilpa", where = asNamespace("latents"))), add = TRUE)
  booted <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3, n_starts = 2,
                                        seed = 2))
  stable <- quietly(sensitivity(fit, seeds = 1:2, n_starts = 2))
  expect_gte(length(seen$effects), 5L)
  expect_true(all(vapply(seen$effects, \(effects) setequal(names(effects), c("y3:z", "y1:z")),
                         logical(1))))
  expect_true(all(is.finite(subset(booted, parameter == "direct_effect")$standard_error)))
})

test_that("direct effects the model cannot hold raise their documented classes", {
  data <- transform(de_unit_data(n = 120), g = rep(seq_len(12), each = 10),
                    t = rep(1:2, 60), constant = 3, gappy = replace(z, 5, NA),
                    text = as.character(round(z)), infinite = replace(z, 2, Inf),
                    item = factor(rep(c("a", "b"), 60)))
  vars <- paste0("y", 1:4)
  fit_with <- function(effects, ..., indicators = vars) {
    lpa(data, indicators, 2, n_starts = 1, seed = 1, direct_effects = effects, ...)
  }
  # Covariates: absent, constant, missing, not numeric, not finite.
  invisible(lapply(c("absent", "constant", "gappy", "text", "infinite", "item"), \(covariate) {
    expect_error(fit_with(sprintf("y3 ~ %s", covariate)), class = "latents_bad_data")
  }))
  # Terms and targets.
  invisible(lapply(list("y3", "~ z", "y3 ~", "y3 ~ + z", "y3 ~ z ~ w", c("y3 ~ z", "y3 ~ z"),
                        "y3 ~ y2", "nope ~ z"), \(effects) {
    expect_error(fit_with(effects), class = "latents_bad_argument")
  }))
  expect_error(fit_with("item ~ z", indicators = c(vars, "item"), categorical = "item"),
               class = "latents_bad_argument")
  # Settings direct effects are not combined with.
  expect_error(fit_with("y3 ~ z", model = "EEE"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", model = "VVV"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", local_dependence = "y1:y2"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", noise = TRUE), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", prior = prior_control()), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", centering = "grand"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", scaling = "grand"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", fixed = "means"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", start = list(means = matrix(0, 2, 4),
                                               variances = matrix(1, 2, 4))),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars, NULL, 2, 1, profile_covariates = "w",
                        direct_effects = "y3 ~ z", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars, "g", 2, 2, group_covariates = "w",
                        direct_effects = "y3 ~ z", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars, "g", 2, 1, time = "t", direct_effects = "y3 ~ z",
                        n_starts = 1, seed = 1), class = "latents_bad_argument")
  plain <- quietly(lpa(data, vars, 2, n_starts = 1, seed = 1))
  expect_error(get_results(plain, "direct_effects"), class = "latents_no_direct_effects")
})

test_that("every malformed `direct_effects` value is refused as a bad argument", {
  # The documented contract: a term that is not "indicator ~ covariate" raises
  # latents_bad_argument. These three used not to be refused that way
  # (FINDINGS.md in equivalence/features/direct-effects, findings 2 and 3): a
  # trailing `+` was dropped by strsplit() and the term fitted as "y3 ~ z";
  # NA and a number stopped with an unclassed stopifnot() error.
  data <- de_unit_data(n = 120)
  invisible(lapply(list("y3 ~ z +", NA_character_, 1), \(effects) {
    outcome <- tryCatch(quietly(lpa(data, paste0("y", 1:4), 2, n_starts = 1, seed = 1,
                                    direct_effects = effects)), error = identity)
    expect_s3_class(outcome, "latents_bad_argument")
  }))
})
