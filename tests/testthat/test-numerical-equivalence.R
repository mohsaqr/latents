source(testthat::test_path("helper-independent-likelihood.R"), local = TRUE)

test_that("exhaustive reference agrees with an analytically solvable model", {
  x <- matrix(c(-0.4, 0.2, 0.8), ncol = 1L)
  parameters <- list(means = matrix(c(0, 0), 2L, 1L),
                     variances = matrix(c(1, 1), 2L, 1L),
                     profile_probabilities = rbind(c(0.8, 0.2), c(0.3, 0.7)),
                     group_probabilities = c(0.4, 0.6))
  reference <- enumerate_latent_assignments(x, c("b", "a", "b"), parameters)
  expect_equal(reference$log_likelihood, sum(stats::dnorm(x, log = TRUE)),
               tolerance = 1e-12)
  expect_equal(reference$group_posteriors, matrix(c(0.4, 0.6), 2L, 2L,
                                                byrow = TRUE))
  expect_equal(reference$subject_posteriors, matrix(0.5, 3L, 2L))
  expect_error(enumerate_latent_assignments(x, "one", parameters))
  expect_error(enumerated_em_update(x, rep("one", 3L), parameters,
                                    variance_model = "invalid"))
  expect_error(optim_multilpa_reference(cbind(x, x), rep("one", 3L), parameters))
  updated <- enumerated_em_update(x, c("b", "a", "b"), parameters)
  expect_equal(dim(updated$variances), c(2L, 1L))
  expect_equal(updated$variances, matrix(mean((x - mean(x))^2), 2L, 1L))
})

test_that("fitted likelihood and posteriors match exhaustive latent assignments", {
  synthetic <- data.frame(group = c("b", "a", "b", "a"),
                          x = c(-1.8, -0.4, 1.2, 2.1),
                          y = c(0.2, 0.8, 2.8, 3.4))
  start <- list(means = rbind(c(-1, 0.5), c(1.5, 3)),
                variances = rbind(c(0.8, 0.7), c(0.9, 0.6)),
                profile_probabilities = rbind(c(0.8, 0.2), c(0.3, 0.7)),
                group_probabilities = c(0.45, 0.55))
  expect_warning(
    expect_warning(
      fit <- multilpa(synthetic, c("x", "y"), "group", 2L,
                        n_group_classes = 2L, start = start,
                        n_starts = 1L, max_iter = 1L),
      "did not converge"),
    "effective membership")
  reference <- enumerate_latent_assignments(as.matrix(synthetic[c("x", "y")]),
                               synthetic$group, fit)
  expect_equal(fit$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  expect_equal(unname(fit$subject_posteriors), reference$subject_posteriors,
               tolerance = 1e-10)
  expect_equal(unname(fit$group_posteriors), reference$group_posteriors,
               tolerance = 1e-10)
  expect_equal(rowSums(fit$subject_posteriors), rep(1, nrow(synthetic)),
               ignore_attr = TRUE)
})

test_that("one EM update matches an independent exhaustive E-step and M-step", {
  synthetic <- data.frame(group = rep(c("c", "a", "b"), each = 2L),
                          x = c(-2, -1, -0.5, 1.5, 0.8, 2.3),
                          y = c(0.1, 0.7, 1.2, 2.9, 2.4, 3.5))
  invisible(lapply(c("varying", "equal"), function(variance_model) {
    stopifnot(is.character(variance_model), length(variance_model) == 1L)
    start <- list(means = rbind(c(-1, 0.5), c(1.5, 3)),
                  variances = rbind(c(0.8, 0.7), c(0.8, 0.7)),
                  profile_probabilities = rbind(c(0.8, 0.2), c(0.3, 0.7)),
                  group_probabilities = c(0.45, 0.55))
    reference <- enumerated_em_update(as.matrix(synthetic[c("x", "y")]),
                                       synthetic$group, start, variance_model)
    expect_warning(
      fit <- multilpa(synthetic, c("x", "y"), "group", 2L,
                        n_group_classes = 2L, start = start,
                        variance_model = variance_model,
                        n_starts = 1L, max_iter = 1L),
      "did not converge")
    invisible(lapply(names(reference), function(parameter) {
      stopifnot(is.character(parameter), length(parameter) == 1L)
      expect_equal(unname(fit[[parameter]]), unname(reference[[parameter]]),
                   tolerance = 1e-9, info = paste(variance_model, parameter))
    }))
  }))
})

test_that("EM agrees with independent direct numerical maximum likelihood", {
  set.seed(481L)
  n_groups <- 40L
  group <- rep(seq_len(n_groups), each = 5L)
  group_class <- rep(c(1L, 2L), each = n_groups / 2L)
  profile <- 1L + (stats::runif(length(group)) > c(0.86, 0.16)[group_class[group]])
  synthetic <- data.frame(group = group,
                          x = stats::rnorm(length(group), c(-2.5, 2.5)[profile],
                                           c(0.6, 0.8)[profile]))
  start <- list(means = matrix(c(-2.5, 2.5), 2L, 1L),
                variances = matrix(c(0.36, 0.64), 2L, 1L),
                profile_probabilities = rbind(c(0.86, 0.14), c(0.16, 0.84)),
                group_probabilities = c(0.5, 0.5))
  reference <- optim_multilpa_reference(as.matrix(synthetic["x"]), group, start)
  fit <- multilpa(synthetic, "x", "group", 2L, n_group_classes = 2L,
                    start = start, n_starts = 1L, max_iter = 3000L, tol = 1e-12)
  expect_equal(reference$optim$convergence, 0L)
  expect_true(fit$converged)
  expect_equal(fit$log_likelihood, -reference$optim$value, tolerance = 1e-7)
  invisible(lapply(names(reference$parameters), function(parameter) {
    stopifnot(is.character(parameter), length(parameter) == 1L)
    expect_equal(unname(fit[[parameter]]), unname(reference$parameters[[parameter]]),
                 tolerance = 2e-4, info = parameter)
  }))
})

test_that("single group class reproduces mclust diagonal Gaussian mixtures", {
  skip_if_not_installed("mclust")
  # Mclust evaluates an unqualified mclustBIC call in its caller's environment.
  mclustBIC <- mclust::mclustBIC
  set.seed(192L)
  synthetic <- data.frame(group = rep(seq_len(30L), each = 4L),
                          x = c(stats::rnorm(60L, -3, 0.6), stats::rnorm(60L, 3, 1)),
                          y = c(stats::rnorm(60L, 0, 0.8), stats::rnorm(60L, 4, 0.5)))
  x <- as.matrix(synthetic[c("x", "y")])
  invisible(lapply(c(varying = "VVI", equal = "EEI"), function(model) {
    stopifnot(is.character(model), length(model) == 1L)
    reference <- mclust::Mclust(x, G = 2L, modelNames = model, verbose = FALSE)
    means <- t(reference$parameters$mean)
    variances <- t(vapply(seq_len(2L), function(profile) {
      stopifnot(length(profile) == 1L)
      diag(reference$parameters$variance$sigma[, , profile])
    }, numeric(ncol(x))))
    start <- list(means = means, variances = variances,
                  profile_probabilities = matrix(reference$parameters$pro, 1L),
                  group_probabilities = 1)
    fit <- multilpa(synthetic, c("x", "y"), "group", 2L,
                      n_group_classes = 1L, start = start, n_starts = 1L,
                      variance_model = if (model == "VVI") "varying" else "equal",
                      tol = 1e-12)
    expect_equal(fit$log_likelihood, as.numeric(reference$loglik), tolerance = 1e-8)
    expect_equal(unname(fit$means), unname(means), tolerance = 1e-6)
    expect_equal(unname(fit$variances), unname(variances), tolerance = 1e-6)
    expect_equal(unname(fit$subject_posteriors), unname(reference$z),
                 tolerance = 1e-6)
  }))
})
