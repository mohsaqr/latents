# Covariances for chosen pairs of continuous indicators
# (`local_dependence = "y1:y2"`): invariances, refits that must keep the
# pattern, and the documented error classes. The comparisons with Latent
# GOLD, Mplus, mclust and an independent implementation are in
# equivalence/test-features-pair-covariance.R.

pair_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  class <- sample.int(2L, n, replace = TRUE, prob = c(0.55, 0.45))
  rho <- c(0.6, -0.5)[class]
  first <- stats::rnorm(n)
  second <- rho * first + sqrt(1 - rho^2) * stats::rnorm(n)
  third <- stats::rnorm(n)
  fourth <- -0.4 * third + sqrt(1 - 0.16) * stats::rnorm(n)
  data.frame(y1 = c(0, 2.2)[class] + first, y2 = c(0, 1.8)[class] + 1.2 * second,
             y3 = c(0, 1.5)[class] + third, y4 = c(0, -1.4)[class] + 0.8 * fourth,
             y5 = stats::rnorm(n, c(0, 1)[class]))
}

pair_fit <- function(data, vars = paste0("y", 1:4), pairs = c("y1:y2", "y3:y4"), ...) {
  quietly(lpa(data, vars, 2, local_dependence = pairs, n_starts = 4, seed = 1,
              tol = 1e-10, max_iter = 5000, ...))
}

# Fitted classes in the order of a reference posterior matrix.
pair_order <- function(posteriors, reference) {
  k <- ncol(posteriors)
  orders <- if (k == 2L) list(1:2, 2:1) else stop("two classes only")
  gaps <- vapply(orders, \(o) max(abs(posteriors[, o] - reference)), numeric(1))
  orders[[which.min(gaps)]]
}

# The pattern a fit must keep: zero outside the pairs, symmetric, every
# eigenvalue at or above the floor, one shared matrix under "equal".
pair_pattern_holds <- function(fit, pairs) {
  vars <- colnames(fit$means)
  free <- diag(length(vars)) > 0
  dimnames(free) <- list(vars, vars)
  ends <- strsplit(pairs, ":", fixed = TRUE)
  free[do.call(rbind, ends)] <- TRUE
  free[do.call(rbind, lapply(ends, rev))] <- TRUE
  matrices <- lapply(seq_len(fit$n_profiles), \(j) fit$covariances[, , j])
  all(vapply(matrices, function(s) {
    # Symmetric to rounding here; exact symmetry is asserted by the
    # equivalence suite (equivalence/features/pair-covariance/FINDINGS.md).
    # The M-step clamps an eigenvalue to the floor exactly, but eigen()
    # recomputes it from the rebuilt matrix to within a few ulps.
    all(s[!free] == 0) && max(abs(s - t(s))) <= 4 * .Machine$double.eps * max(abs(s)) &&
      min(eigen(s, symmetric = TRUE, only.values = TRUE)$values) >=
        fit$min_variance * (1 - 8 * .Machine$double.eps)
  }, logical(1))) &&
    (!identical(fit$variance_model, "equal") ||
       all(vapply(matrices, identical, logical(1), matrices[[1L]])))
}

test_that("the fit keeps the pattern and counts one parameter per class and pair", {
  data <- pair_data()
  fit <- pair_fit(data)
  expect_true(pair_pattern_holds(fit, c("y1:y2", "y3:y4")))
  diagonal <- quietly(lpa(data, paste0("y", 1:4), 2, model = "VVI", n_starts = 4, seed = 1))
  expect_equal(fit$n_parameters, diagonal$n_parameters + 2 * 2)
  equal <- pair_fit(data, variance_model = "equal")
  expect_true(pair_pattern_holds(equal, c("y1:y2", "y3:y4")))
  shared <- quietly(lpa(data, paste0("y", 1:4), 2, model = "EEI", n_starts = 4, seed = 1))
  expect_equal(equal$n_parameters, shared$n_parameters + 2)
  # Freeing covariances can only raise the maximum of the nested diagonal model.
  expect_gte(fit$log_likelihood, diagonal$log_likelihood - 1e-8)
  expect_identical(get_results(fit, "model")$covariance_structure, "VVI + y1:y2, y3:y4")
  expect_identical(get_results(equal, "model")$covariance_structure, "EEI + y1:y2, y3:y4")
})

test_that("the order of rows, pairs, pair ends and indicators does not matter", {
  data <- pair_data()
  fit <- pair_fit(data)
  set.seed(3)
  permutation <- sample.int(nrow(data))
  shuffled <- pair_fit(data[permutation, ])
  expect_equal(shuffled$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  order <- pair_order(shuffled$subject_posteriors, fit$subject_posteriors[permutation, ])
  expect_equal(unname(shuffled$subject_posteriors[, order]),
               unname(fit$subject_posteriors[permutation, ]), tolerance = 1e-6)
  expect_equal(unname(shuffled$covariances[, , order]), unname(fit$covariances),
               tolerance = 1e-6)
  reversed <- pair_fit(data, pairs = c("y4:y3", "y2:y1"))
  expect_equal(reversed$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  expect_identical(reversed$n_parameters, fit$n_parameters)
  expect_identical(reversed$covariance_pattern, fit$covariance_pattern)
  order <- pair_order(reversed$subject_posteriors, fit$subject_posteriors)
  expect_equal(unname(reversed$covariances[, , order]), unname(fit$covariances),
               tolerance = 1e-6)
  columns <- pair_fit(data, vars = c("y4", "y2", "y3", "y1"))
  expect_equal(columns$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  order <- pair_order(columns$subject_posteriors, fit$subject_posteriors)
  vars <- paste0("y", 1:4)
  expect_equal(unname(columns$covariances[vars, vars, order]), unname(fit$covariances),
               tolerance = 1e-6)
  expect_equal(unname(columns$means[order, vars]), unname(fit$means), tolerance = 1e-6)
})

test_that("an affine change of an indicator's units moves its estimates and the likelihood by the Jacobian", {
  data <- pair_data()
  fit <- pair_fit(data)
  scale <- -250
  shift <- 40
  moved <- transform(data, y1 = shift + scale * y1)
  refit <- pair_fit(moved)
  # log |dy1'/dy1| per row.
  expect_equal(refit$log_likelihood, fit$log_likelihood - nrow(data) * log(abs(scale)),
               tolerance = 1e-9)
  order <- pair_order(refit$subject_posteriors, fit$subject_posteriors)
  expect_equal(unname(refit$subject_posteriors[, order]), unname(fit$subject_posteriors),
               tolerance = 1e-6)
  expect_equal(unname(refit$means[order, "y1"]), shift + scale * unname(fit$means[, "y1"]),
               tolerance = 1e-6)
  expect_equal(unname(refit$covariances["y1", "y1", order]),
               scale^2 * unname(fit$covariances["y1", "y1", ]), tolerance = 1e-6)
  expect_equal(unname(refit$covariances["y1", "y2", order]),
               scale * unname(fit$covariances["y1", "y2", ]), tolerance = 1e-6)
  expect_equal(unname(refit$covariances["y3", "y4", order]),
               unname(fit$covariances["y3", "y4", ]), tolerance = 1e-6)
})

test_that("integer weights equal duplicated rows, up to the documented scaling", {
  data <- pair_data(n = 200)
  set.seed(5)
  data$w <- sample(1:3, nrow(data), replace = TRUE)
  weighted <- pair_fit(data, weights = "w")
  duplicated <- pair_fit(data[rep(seq_len(nrow(data)), data$w), ])
  # Weights are scaled to sum to the number of rows, so the pseudo log
  # likelihood is the duplicated one times n / sum(w).
  expect_equal(weighted$log_likelihood * sum(data$w) / nrow(data),
               duplicated$log_likelihood, tolerance = 1e-9)
  order <- pair_order(weighted$subject_posteriors,
                      duplicated$subject_posteriors[!duplicated(rep(seq_len(nrow(data)),
                                                                    data$w)), ])
  expect_equal(unname(weighted$covariances[, , order]), unname(duplicated$covariances),
               tolerance = 1e-6)
  expect_equal(unname(weighted$means[order, ]), unname(duplicated$means), tolerance = 1e-6)
})

test_that("predict() on the training rows reproduces the posteriors", {
  data <- pair_data()
  fit <- pair_fit(data)
  predicted <- predict(fit, newdata = data, type = "posterior")
  expect_equal(unname(as.matrix(predicted[c("probability_profile_1", "probability_profile_2")])),
               unname(fit$subject_posteriors), tolerance = 1e-10)
  densities <- predict(fit, newdata = data, type = "density")
  expect_equal(sum(densities$log_density), fit$log_likelihood, tolerance = 1e-10)
})

test_that("refits keep the pattern: refit arguments, enumeration, comparison, bootstrap", {
  skip_on_cran()
  data <- pair_data(n = 250)
  fit <- pair_fit(data, pairs = "y1:y2", vars = paste0("y", 1:3))
  arguments <- .multilpa_refit_arguments(fit)
  expect_identical(arguments$local_dependence, "y1:y2")
  expect_identical(arguments$variance_model, "varying")
  expect_null(arguments$covariance_model)
  measurement <- .multilpa_refit_measurement(fit)
  expect_identical(measurement$local_dependence, "y1:y2")
  equal <- pair_fit(data, pairs = "y1:y2", vars = paste0("y", 1:3), variance_model = "equal")
  expect_identical(.multilpa_refit_arguments(equal)$variance_model, "equal")
  expect_identical(.latents_pattern_label(equal), "EEI + y1:y2")
  compared <- compare_models(varying = fit, equal = equal)
  expect_identical(compared$structure, c("VVI + y1:y2", "EEI + y1:y2"))
  expect_equal(compared$n_parameters, c(fit$n_parameters, equal$n_parameters))
  grid <- as.data.frame(quietly(enumerate_lpa(data, paste0("y", 1:3), n_profiles = 1:2,
                                              local_dependence = "y1:y2", n_starts = 4,
                                              seed = 1, tol = 1e-10, max_iter = 5000)))
  cell <- subset(grid, n_profiles == 2 & model == "VVI + y1:y2")
  expect_equal(cell$log_likelihood, fit$log_likelihood, tolerance = 1e-8)
  expect_equal(cell$n_parameters, fit$n_parameters)
  # A bootstrap replicate refits the same pattern: with the same rows and the
  # refit arguments it reproduces the fit's maximum and its count.
  refit <- quietly(do.call(multilpa, c(list(data = .multilpa_resolve_data(fit, NULL),
                                            n_starts = 4, seed = 1, tol = 1e-10,
                                            max_iter = 5000), arguments)))
  expect_identical(refit$covariance_pattern, fit$covariance_pattern)
  expect_equal(refit$n_parameters, fit$n_parameters)
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-8)
  # Every refit the resampling verbs make carries the pairs and the variance
  # model, never a full covariance (recorded through the namespace binding).
  original <- multilpa
  seen <- list()
  local_mocked_bindings(multilpa = function(...) {
    arguments <- list(...)
    seen[[length(seen) + 1L]] <<- list(pairs = arguments$local_dependence,
                                       variance = arguments$variance_model,
                                       covariance = arguments$covariance_model,
                                       model = arguments$model)
    original(...)
  })
  boot <- quietly(parameter_inference(fit, method = "bootstrap", iter = 5L, seed = 1))
  expect_true(all(is.finite(subset(boot, parameter == "covariance" &
                                     term == "y2:y1")$standard_error)))
  expect_false(any(subset(boot, parameter == "covariance")$term %in% c("y3:y1", "y3:y2")))
  null <- quietly(lpa(data, paste0("y", 1:3), 1, local_dependence = "y1:y2", n_starts = 1,
                      seed = 1))
  seen <- list()
  lrt <- quietly(bootstrap_lrt(null, fit, iter = 3L, n_starts = 2L, seed = 1))
  expect_identical(lrt$n_valid, 3L)
  # Two refits (null and alternative) per replicate, each with the pattern.
  expect_length(seen, 6L)
  expect_true(all(vapply(seen, \(call) identical(call$pairs, "y1:y2") &&
                           identical(call$variance, "varying") && is.null(call$covariance) &&
                           is.null(call$model), logical(1))))
})

test_that("a variance floor binds on a block's eigenvalue and keeps the pattern", {
  set.seed(21)
  n <- 300
  class <- sample.int(2L, n, replace = TRUE)
  first <- stats::rnorm(n)
  tight <- class == 1L
  # In the first class y1 and y2 are almost the same variable at a small
  # scale; the pair's smaller eigenvalue is about 1e-6.
  data <- data.frame(y1 = ifelse(tight, 0.05 * first, 1 + stats::rnorm(n)),
                     y2 = ifelse(tight, 0.05 * first + 0.001 * stats::rnorm(n),
                                 1 + stats::rnorm(n)),
                     y3 = stats::rnorm(n, c(0, 2)[class]))
  fit <- quietly(lpa(data, c("y1", "y2", "y3"), 2, local_dependence = "y1:y2",
                     min_variance = 1e-3, n_starts = 4, seed = 1, tol = 1e-10))
  expect_true(fit$boundary)
  expect_true(pair_pattern_holds(fit, "y1:y2"))
  smallest <- min(vapply(1:2, \(j) min(eigen(fit$covariances[1:2, 1:2, j], symmetric = TRUE,
                                             only.values = TRUE)$values), numeric(1)))
  expect_equal(smallest, 1e-3, tolerance = 1e-8)
  expect_error(parameter_inference(fit), class = "latents_boundary_fit")
})

test_that("pairs the model cannot hold are refused by their documented class", {
  data <- transform(pair_data(n = 120), item = factor(rep(c("a", "b"), 60)),
                    g = rep(seq_len(12), each = 10), t = rep(1:2, 60))
  vars <- paste0("y", 1:4)
  refused <- function(...) {
    expect_error(lpa(data, vars, 2, n_starts = 1, seed = 1, ...),
                 class = "latents_bad_argument")
  }
  refused(local_dependence = c("y1:y2", "y2:y3"))
  refused(local_dependence = "y1:y1")
  refused(local_dependence = "y1-y2")
  refused(local_dependence = "y1:")
  refused(local_dependence = NA_character_)
  refused(local_dependence = "y1:y9")
  refused(local_dependence = "y1:y2", covariance_model = "full")
  refused(local_dependence = "y1:y2", covariance_model = "diagonal")
  refused(local_dependence = "y1:y2", model = "VVI")
  refused(local_dependence = "y1:y2", volume = "equal")
  refused(local_dependence = "y1:y2", prior = TRUE)
  refused(local_dependence = "y1:y2", noise = TRUE)
  refused(local_dependence = "y1:y2", direct_effects = "y3 ~ y5")
  expect_error(lpa(data, c(vars, "item"), 2, categorical = "item",
                   local_dependence = "y1:item", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars, "g", 2, 2, profile_covariates = "y5",
                        local_dependence = "y1:y2", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars, "g", 2, 1, time = "t", local_dependence = "y1:y2",
                        n_starts = 1, seed = 1), class = "latents_bad_argument")
  fit <- pair_fit(pair_data(n = 120), pairs = "y1:y2")
  expect_error(lpa(data, vars, 2, local_dependence = "y1:y2", fixed = "variances",
                   start = fit, n_starts = 1, seed = 1), class = "latents_bad_argument")
  # A pair item as a covariate is refused by the same rule, before the
  # indicator-reuse check sees it.
  expect_error(multilpa(data, vars, "g", 2, 2, profile_covariates = "y1",
                        local_dependence = "y1:y2", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
})

test_that("a defect that leaves a covariance outside the pairs is detected", {
  data <- pair_data(n = 200)
  # Sabotage in memory: the M-step forgets the pattern.
  original <- .multilpa_maximize_moments
  local_mocked_bindings(.multilpa_maximize_moments = function(..., covariance_pattern = NULL) {
    original(..., covariance_pattern = NULL)
  })
  broken <- pair_fit(data)
  expect_false(pair_pattern_holds(broken, c("y1:y2", "y3:y4")))
})
