# Regression checks from the independent package-wide audit. Known answers
# are raw binomial/Gaussian formulas or invariances, rather than saved outputs.

test_that("Gaussian dispersion identifies labels when means are held equal", {
  withr::local_seed(882)
  z <- stats::rnorm(1000)
  data <- list(data.frame(y = z * c(rep(1, 400), rep(5, 600))),
               data.frame(y = z * c(rep(1, 600), rep(5, 400))))
  start <- starting_values(quietly(lpa(data[[1]], "y", 2, n_starts = 1, seed = 3)))
  start$means[,] <- 0
  start$variances[,] <- c(1, 25)
  start$profile_probabilities[,] <- .5
  arguments <- list(vars = "y", id = NULL, n_profiles = 2, n_group_classes = 1,
                    start = start, fixed = "means", n_starts = 1, max_iter = 5000,
                    tol = 1e-12, acceleration = "none", search = "kmeans", seed = 3)
  fits <- lapply(data, function(rows) quietly(do.call(multilpa,
    c(list(data = rows), arguments))))
  pooled <- quietly(do.call(pool_imputations, c(list(imputed = data), arguments)))
  expected <- rowMeans(vapply(fits, function(fit) sort(as.numeric(fit$variances)),
                              numeric(2)))
  variances <- subset(pooled$estimates, parameter == "variance")
  expect_equal(variances$estimate, expected, tolerance = 1e-8)
  expect_lt(variances$estimate[1], 2)
  expect_gt(variances$estimate[2], 20)
  expect_identical(pooled$fits$profile_order, c("1,2", "1,2"))
})

test_that("covariance orientation also identifies otherwise equal Gaussian profiles", {
  block <- list(means = matrix(0, 2, 2), variances = matrix(1, 2, 2),
                covariances = array(c(1, .8, .8, 1, 1, -.8, -.8, 1), c(2, 2, 2)))
  moved <- block
  moved$covariances <- block$covariances[, , 2:1, drop = FALSE]
  reference <- .latents_measurement_signature(list(block), list(block))
  candidate <- .latents_measurement_signature(list(moved), list(block))
  expect_identical(.multilpa_match_order(reference, candidate), 2:1)
  expect_equal(candidate, reference[2:1, , drop = FALSE])
})

test_that("categorical labels with delimiters retain their own errors and intervals", {
  data <- data.frame(a = rep(c("b:c", "b:d"), c(60, 40)), check.names = FALSE)
  data[["a:b"]] <- rep(c("c", "d"), c(20, 80))
  fit <- lca(data, c("a:b", "a"), 1, n_starts = 1, seed = 1)
  inference <- parameter_inference(fit)
  responses <- get_results(fit, "responses")
  expected <- sqrt(responses$probability * (1 - responses$probability) / nrow(data))
  expect_equal(responses$probability_standard_error, expected, tolerance = 1e-8)
  expect_false(anyDuplicated(names(coef(fit))) > 0)
  expect_false(anyDuplicated(rownames(vcov(fit))) > 0)
  expect_equal(unname(confint(fit, parm = 3)),
               matrix(unlist(inference[3, c("conf_low", "conf_high")]), 1),
               tolerance = 1e-12)
  expect_equal(unname(confint(fit)),
               unname(as.matrix(inference[c("conf_low", "conf_high")])),
               tolerance = 1e-12)
  # Escaping the escape character prevents literal encoded-looking names
  # from colliding with a real delimiter in another component.
  expect_identical(.latents_parameter_component(c("a:b", "a%3Ab", "c,d", "x[y]")),
                   c("a%3Ab", "a%253Ab", "c%2Cd", "x%5By%5D"))
})

test_that("a univariate full-covariance start retains one row per profile", {
  withr::local_seed(51)
  data <- data.frame(y = c(stats::rnorm(120, -3), stats::rnorm(120, 3)))
  fit <- quietly(lpa(data, "y", 2, model = "VVV", n_starts = 1, seed = 1))
  start <- starting_values(fit)
  table <- get_results(start, "covariances")
  expect_identical(table$profile, 1:2)
  expect_identical(table$indicator, c(1L, 1L))
  expect_identical(table$indicator_2, c(1L, 1L))
  expect_equal(table$covariance, as.numeric(fit$variances))
})

test_that("BIC and nobs use the fit's informative-group convention", {
  withr::local_seed(84)
  data <- data.frame(y = c(stats::rnorm(60), rep(NA_real_, 20)))
  invisible(lapply(c(TRUE, FALSE), function(drop) {
    fit <- quietly(lpa(data, "y", 1, missing = "fiml", drop_empty = drop,
                        n_starts = 1, seed = 1))
    expect_equal(nobs(fit), if (drop) 60 else 80)
    expect_equal(attr(logLik(fit), "nobs"), nobs(fit))
    expect_equal(stats::BIC(fit), fit$bic, tolerance = 1e-12)
  }))
  # Both LTA engines use the same group-count convention, including saved
  # objects whose information count differs from their original row count.
  invisible(lapply(c("multilpa_transitions", "multilpa_lta"), function(kind) {
    object <- structure(list(log_likelihood = -10, n_parameters = 2,
      n_groups = 5L, n_informative_groups = 3L, drop_empty = TRUE), class = kind)
    expect_identical(nobs(object), 3L)
    expect_identical(attr(logLik(object), "nobs"), 3L)
    expect_equal(stats::BIC(object), 20 + 2 * log(3))
    object$drop_empty <- FALSE
    expect_identical(nobs(object), 5L)
  }))
})

test_that("large indicator origins do not relax data identity checks", {
  withr::local_seed(19)
  data <- data.frame(y = 1e10 + stats::rnorm(100, sd = 10))
  fit <- quietly(lpa(data, "y", 1, n_starts = 1, seed = 1))
  changed <- data
  changed$y[1] <- changed$y[1] + 1
  expect_error(parameter_inference(fit, changed), class = "latents_bad_inference_data")
  expect_error(get_results(fit, "assignments", data = changed),
                class = "latents_bad_inference_data")
  other <- quietly(lpa(changed, "y", 1, n_starts = 1, seed = 1))
  expect_error(compare_models(original = fit, changed = other),
                class = "latents_incomparable_models")
  expect_false(.multilpa_same_values(data$y, changed$y))
  expect_true(.multilpa_same_values(c(1L, NA_integer_), c(1, NA_real_)))
  expect_false(.multilpa_same_values(c("0.3", "0.3"), c(.3, .1 + .2)))
  expect_false(.multilpa_same_values(c(.3, .1 + .2), c("0.3", "0.3")))
  expect_true(.multilpa_same_values(c("0.3", "0.4"), c(.3, .4)))
  expect_true(.multilpa_same_values(factor(c("a", "b")), c("a", "b")))
  expect_false(.multilpa_alignment_values(c(.3, .3), c(.3, .1 + .2)))
  original <- c(pi, exp(1), 1 / 3)
  connection <- textConnection("encoded", "w", local = TRUE)
  utils::write.table(data.frame(value = original), connection,
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(connection)
  serialized <- as.numeric(encoded)
  expect_true(.multilpa_alignment_values(serialized, original))
  expect_false(.multilpa_same_values(serialized, original))
  transformed <- quietly(lpa(data, "y", 1, centering = "grand", scaling = "grand",
                              n_starts = 1, seed = 1))
  expect_silent(.multilpa_check_alignment(transformed, data))
  expect_silent(.multilpa_check_alignment(transformed, .multilpa_model_frame(transformed)))
  expect_s3_class(parameter_inference(transformed), "data.frame")
  expect_error(.multilpa_check_alignment(transformed, changed),
                class = "latents_bad_inference_data")
})

test_that("assignment alignment recognizes R's own CSV encoding at rounding ties", {
  withr::local_seed(193)
  data <- data.frame(id = rep(seq_len(20), each = 4), t = rep(1:4, 20),
                     a = stats::rnorm(80), b = stats::rnorm(80))
  fit <- multilpa(data, c("a", "b"), "id", 1, 1, n_starts = 1, seed = 1)
  connection <- textConnection("csv", "w", local = TRUE)
  utils::write.csv(data, connection, row.names = FALSE)
  close(connection)
  read_back <- utils::read.csv(text = paste(csv, collapse = "\n"))
  expect_identical(get_results(fit, "assignments", data = read_back)$profile,
                   fit$subject_profiles)
  # The CSV copy reproduces the fit, so the Wald path accepts it too and
  # evaluates the stored values: the same errors as without `data` (review
  # 2026-10-26 round 2, INF-3). A changed value is still refused.
  expect_equal(parameter_inference(fit, read_back)$standard_error,
               parameter_inference(fit)$standard_error)
  read_back$a[1] <- read_back$a[1] + 1e-6
  expect_error(parameter_inference(fit, read_back), class = "latents_bad_inference_data")
})

test_that("numeric ordinal alignment uses native values rather than rounded labels", {
  data <- data.frame(o = rep(c(1, 1e15, 1e15 + 1), 30))
  fit <- lca(data, "o", 1, ordinal = "o", n_starts = 1, seed = 1)
  expect_silent(.multilpa_check_alignment(fit, data))
  expect_identical(get_results(fit, "assignments", data = data)$profile,
                   fit$subject_profiles)
  collapsed <- data
  collapsed$o[collapsed$o == 1e15 + 1] <- 1e15
  expect_error(.multilpa_check_alignment(fit, collapsed),
                class = "latents_bad_inference_data")
})

test_that("additive inference is invariant to indicator units and origin", {
  withr::local_seed(123)
  data <- data.frame(g = rep(seq_len(100), each = 5))
  data$y <- rep(stats::rnorm(100, sd = 2), each = 5) + stats::rnorm(500)
  fit <- multilpa(data, "y", "g", n_group_classes = 1, family = "additive",
                  n_starts = 1, seed = 1, tol = 1e-12)
  baseline <- parameter_inference(fit)
  invisible(lapply(list(c(multiplier = 1e6, shift = 0),
                       c(multiplier = 1, shift = 1e8)), function(change) {
    moved <- data
    moved$y <- moved$y * change[["multiplier"]] + change[["shift"]]
    other <- multilpa(moved, "y", "g", n_group_classes = 1, family = "additive",
                      n_starts = 1, seed = 1, tol = 1e-12)
    inference <- parameter_inference(other)
    factor <- ifelse(baseline$parameter == "mean", change[["multiplier"]],
                     ifelse(baseline$parameter == "variance", change[["multiplier"]]^2, 1))
    expect_equal(inference$standard_error / factor, baseline$standard_error,
                 tolerance = 1e-7)
    expect_equal(attr(inference, "condition_ratio"), attr(baseline, "condition_ratio"),
                 tolerance = 1e-7)
  }))
  # Analytic normal random-intercept information in the balanced one-class case.
  n <- 5
  expected_mean_se <- sqrt((fit$between_variances[1] + fit$within_variances[1] / n) / 100)
  expect_equal(baseline$standard_error[baseline$parameter == "mean"],
               expected_mean_se, tolerance = 1e-9)
})

test_that("additive supplied data and weights must reproduce sufficient statistics", {
  withr::local_seed(993)
  data <- data.frame(g = rep(seq_len(60), each = 5), w = rep(c(1, 2), each = 150))
  data$y <- 1e10 + rep(stats::rnorm(60, sd = 2), each = 5) + stats::rnorm(300)
  fit <- multilpa(data, "y", "g", n_group_classes = 1, family = "additive",
                  weights = "w", n_starts = 1, seed = 1)
  changed <- data
  changed$y[changed$g == 1] <- changed$y[changed$g == 1] + 1
  expect_error(parameter_inference(fit, changed), class = "latents_bad_inference_data")
  changed <- data
  changed$w[changed$g == 1] <- 3
  expect_error(parameter_inference(fit, changed), class = "latents_bad_inference_data")
  rescaled <- data
  rescaled$w <- rescaled$w * 1e-12
  expect_s3_class(parameter_inference(fit, rescaled), "data.frame")
})

test_that("comparison data retain equivalent normalized sampling-weight units", {
  withr::local_seed(552)
  data <- data.frame(y = stats::rnorm(150), z = stats::rnorm(150),
                     w = exp(stats::rnorm(150)))
  rescaled <- data
  rescaled$w <- rescaled$w * 1e-12
  profile <- lapply(list(data, rescaled), function(rows) quietly(
    lpa(rows, "y", 1, weights = "w", n_starts = 1, seed = 1)))
  regression <- lapply(list(data, rescaled), function(rows) {
    mixture_regression(y ~ z, rows, 1, weights = "w", n_starts = 0)
  })
  expect_s3_class(compare_models(profile[[1]], profile[[2]]), "data.frame")
  expect_s3_class(compare_models(regression[[1]], regression[[2]]), "data.frame")
})

test_that("comparison identity proves centered LTA representations by the forward map", {
  raw <- matrix(c(.1, 1 / 3, pi), ncol = 1L)
  centers <- colMeans(raw)
  centered <- sweep(raw, 2L, centers, "-")
  original <- list(indicators = raw, indicator_centers = NULL,
                   centered_indicators = NULL)
  general <- list(indicators = sweep(centered, 2L, centers, "+"),
                  indicator_centers = centers, centered_indicators = centered)
  expect_true(.latents_same_comparison_data(original, general))
  expect_true(.latents_same_comparison_data(general, original))
  general$centered_indicators[1, 1] <- general$centered_indicators[1, 1] + 1e-6
  general$indicators[1, 1] <- general$indicators[1, 1] + 1e-6
  expect_false(.latents_same_comparison_data(original, general))
})

test_that("large common count densities retain the occasion's profile prior", {
  object <- list(indicator_data = matrix(NA_real_, 1, 1), means = matrix(0, 2, 1),
    variances = matrix(1, 2, 1), categorical_data = NULL,
    extra_data = list(ordinal = NULL, count = matrix(1e16, 1, 1)),
    count_means = matrix(1.5, 2, 1), subject_posteriors = matrix(c(.8, .2), 1),
    informative = TRUE)
  expect_equal(unname(.lta_occasion_posteriors(object)), matrix(c(.8, .2), 1),
               tolerance = 1e-12)
  general <- object
  class(general) <- "multilpa_lta"
  general$n_profiles <- 2L
  general$x <- object$indicator_data
  general$continuous <- "y"
  general$centers <- 0
  general$occasion_of_row <- 1L
  general$measurement <- list(object[c("means", "variances", "count_means")])
  general$expectation <- list(subject_posteriors = object$subject_posteriors)
  expect_equal(unname(.lta_occasion_posteriors(general)), matrix(c(.8, .2), 1),
               tolerance = 1e-12)
})

test_that("empty persons do not enlarge step-three robust degrees of freedom", {
  withr::local_seed(451)
  n <- 60
  z <- stats::rnorm(n)
  paths <- matrix(sample(c(-3, 3), n * 4, TRUE), n, 4)
  data <- data.frame(id = rep(seq_len(n), each = 4), time = rep(1:4, n),
    z = rep(z, each = 4), y = as.vector(t(paths)) + stats::rnorm(n * 4, sd = .5))
  empty <- data.frame(id = rep(61:100, each = 4), time = rep(1:4, 40), z = 0, y = NA_real_)
  fit_one <- function(rows) lta(rows, "y", "id", 2, time = "time", occasions = "grid",
    missing = "fiml", n_starts = 1, seed = 1, tol = 1e-11)
  base <- fit_one(data)
  padded <- fit_one(rbind(data, empty))
  expected <- r3step(base, data, initial_covariates = "z")
  result <- r3step(padded, rbind(data, empty), initial_covariates = "z")
  expect_equal(base$log_likelihood, padded$log_likelihood, tolerance = 1e-10)
  expect_equal(result$estimate, expected$estimate, tolerance = 1e-6)
  expect_equal(result$standard_error, expected$standard_error, tolerance = 1e-6)
})

test_that("LTA three-step covariate standardization ignores the choice of origin", {
  values <- seq(-2, 2, length.out = 100)
  basis <- .lta_three_step_basis(cbind(1, values))
  translated <- .lta_three_step_basis(cbind(1, values + 1e8))
  expect_equal(unname(cbind(1, values) %*% basis),
               unname(cbind(1, values + 1e8) %*% translated), tolerance = 1e-7)
  expect_error(.lta_three_step_basis(cbind(1, rep(1e8, 100))),
                class = "latents_bad_covariate")
})

test_that("LTA BLRT requires identical variance and probability constraints", {
  base <- list(arguments = list(vars = "y", id = "id", time = "time",
    n_profiles = 2L, n_group_classes = 1L, variance_model = "varying",
    categorical = character(), occasions = "grid", transitions = "homogeneous",
    transition_covariates = character(), initial_covariates = character(),
    measurement = "invariant", order = 1L, mover_stayer = FALSE,
    missing = "fiml", covariance_model = "diagonal", model = "VVI",
    min_variance = 1e-6, min_probability = 1e-10))
  class(base) <- c("multilpa_lta", "multilpa")
  alternative <- base
  alternative$arguments$n_profiles <- 3L
  expect_true(.lta_nested(base, alternative))
  alternative$arguments$min_variance <- .1
  expect_false(.lta_nested(base, alternative))
  alternative$arguments$min_variance <- base$arguments$min_variance
  alternative$arguments$min_probability <- .01
  expect_false(.lta_nested(base, alternative))
})

test_that("an all-missing LTA occasion has no fabricated step-three assignment", {
  withr::local_seed(12)
  n <- 200L
  w <- stats::rnorm(n)
  z <- matrix(stats::rnorm(n * 4), n)
  first <- ifelse(stats::runif(n) < stats::plogis(.3 + .8 * w), 1L, 2L)
  paths <- Reduce(function(previous, time) {
    logit <- ifelse(previous == 1L, -1.2 + .8 * z[, time], -1 - .6 * z[, time])
    ifelse(stats::runif(n) < stats::plogis(logit), 3L - previous, previous)
  }, 2:4, init = first, accumulate = TRUE)
  state <- as.vector(t(do.call(cbind, paths)))
  data <- data.frame(id = rep(seq_len(n), each = 4), time = rep(1:4, n),
                     w = rep(w, each = 4))
  data$y1 <- stats::rnorm(nrow(data), c(-1.2, 1.2)[state], 1)
  data$y2 <- stats::rnorm(nrow(data), c(-1, 1)[state], 1)
  missing <- data$time == 2 & data$id %% 3 == 0
  data[missing, c("y1", "y2")] <- NA_real_
  removed <- data[!missing, ]
  fit_one <- function(rows) lta(rows, c("y1", "y2"), "id", 2, time = "time",
    occasions = "grid", missing = "fiml", n_starts = 1, seed = 1, tol = 1e-11)
  present <- fit_one(data)
  absent <- fit_one(removed)
  inputs <- .lta_three_step_inputs(present, data, "w", character())
  expect_equal(sum(is.na(inputs$modal)), sum(missing))
  expect_true(all(is.na(inputs$modal[missing])))
  expect_equal(present$log_likelihood, absent$log_likelihood, tolerance = 1e-10)
  result <- r3step(present, data, initial_covariates = "w")
  reference <- r3step(absent, removed, initial_covariates = "w")
  expect_equal(attr(result, "classification_errors"),
               attr(reference, "classification_errors"), tolerance = 1e-6)
  expect_equal(result$estimate, reference$estimate, tolerance = 1e-6)
  expect_equal(result$standard_error, reference$standard_error, tolerance = 1e-6)
  translated <- data
  translated$w <- translated$w + 1e8
  moved <- r3step(present, translated, initial_covariates = "w")
  slope <- result$term == "w"
  expect_equal(moved$estimate[slope], result$estimate[slope], tolerance = 1e-6)
  expect_equal(moved$standard_error[slope], result$standard_error[slope], tolerance = 1e-6)
  data$sampling <- rep(seq_len(n), each = 4)
  weighted <- lta(data, c("y1", "y2"), "id", 2, time = "time",
    occasions = "grid", missing = "fiml", weights = "sampling", n_starts = 1,
    seed = 1, tol = 1e-10)
  expect_silent(.lta_check_data(weighted, data))
  changed <- data
  changed$sampling[changed$id == 1] <- 2 * changed$sampling[changed$id == 1]
  expect_error(.lta_check_data(weighted, changed), class = "latents_bad_inference_data")
  expect_false(.lta_nested(present, weighted))
})
