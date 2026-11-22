# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# The uniform noise component (`multilpa(noise = TRUE)`). The fourteen-structure
# comparison with mclust lives in equivalence/test-mclust-noise.R.

noise_frame <- function(seed = 3L) {
  set.seed(seed)
  x <- rbind(matrix(stats::rnorm(135L * 3L, 2), 135L, 3L),
             matrix(stats::rnorm(135L * 3L, -2), 135L, 3L),
             matrix(stats::runif(30L * 3L, -10, 10), 30L, 3L))
  colnames(x) <- c("a", "b", "c")
  as.data.frame(x)
}

noise_fit <- function(data, ...) {
  quietly(multilpa(data, c("a", "b", "c"), id = NULL, n_profiles = 2L,
                   noise = TRUE, n_starts = 3L, seed = 1L, tol = 1e-10, ...))
}

test_that("the hypervolume is mclust's", {
  data <- as.matrix(noise_frame())
  # The axis-aligned box, by hand, is an upper bound on it (a unit check).
  box <- sum(log(apply(data, 2L, function(column) diff(range(column)))))
  expect_lte(latents:::.multilpa_log_hypervolume(data), box + 1e-12)
  skip_without_reference("mclust")
  expect_reference(exp(latents:::.multilpa_log_hypervolume(data)), mclust::hypvol(data),
                   reference = "mclust", covers = "multilpa()",
                   quantity = "hypervolume of the data", tolerance = 1e-12, measure = "relative",
                   note = "Package side: the internal hypervolume (.multilpa_log_hypervolume) that sets the uniform density of multilpa(noise = TRUE); reference: mclust::hypvol().")
})

test_that("a noise fit reproduces mclust's EM from the same start", {
  skip_without_reference("mclust")
  skip_on_cran()
  data <- noise_frame()
  x <- as.matrix(data)
  set.seed(5)
  flagged <- c(rep(FALSE, 270L), rep(TRUE, 30L))
  clusters <- stats::kmeans(x[!flagged, ], 2L)$cluster
  z <- matrix(0, nrow(x), 3L)
  z[cbind(which(!flagged), clusters)] <- 1
  z[flagged, 3L] <- 1
  inverse_volume <- mclust::hypvol(x, reciprocal = TRUE)
  runs <- lapply(c("VVI", "VEV"), function(code) {
    run <- function(itmax, tol) {
      do.call(mclust::me, list(modelName = code, data = x, z = z, Vinv = inverse_volume,
                               control = mclust::emControl(itmax = c(itmax, 1e6),
                                                           tol = c(tol, 1e-14))),
              envir = asNamespace("mclust"))
    }
    first <- run(1, 1e-5)$parameters
    reference <- run(1e6, 1e-13)
    start <- list(means = t(first$mean), profile_probabilities = matrix(first$pro[1:2], 1L),
                  group_probabilities = 1, noise_probability = first$pro[3L])
    if (endsWith(code, "I")) {
      start$variances <- t(vapply(1:2, function(k) diag(first$variance$sigma[, , k]),
                                  numeric(3L)))
    } else start$covariances <- first$variance$sigma
    fit <- quietly(do.call(multilpa, c(
      list(data = data, vars = c("a", "b", "c"), id = NULL, n_profiles = 2L,
           n_starts = 1L, start = start, tol = 1e-13, max_iter = 20000L,
           noise = TRUE, acceleration = "none"),
      latents:::.multilpa_structure_arguments(code))))
    cells <- expand.grid(variable = c("a", "b", "c"), profile = 1:2, stringsAsFactors = FALSE)
    list(actual = data.frame(code = code, log_likelihood = fit$log_likelihood,
                             noise_probability = fit$noise_probability,
                             n_parameters = fit$n_parameters),
         expected = data.frame(code = code, log_likelihood = reference$loglik,
                               noise_probability = reference$parameters$pro[3L],
                               n_parameters = mclust::nMclustParams(code, 3L, 2L, noise = TRUE)),
         actual_means = data.frame(code = code, cells, mean = as.vector(t(fit$means))),
         expected_means = data.frame(code = code, cells,
                                     mean = as.vector(reference$parameters$mean)))
  })
  stack <- function(part) do.call(rbind, lapply(runs, `[[`, part))
  actual <- stack("actual")
  expected <- stack("expected")
  expect_reference(actual[c("code", "log_likelihood")], expected[c("code", "log_likelihood")],
                   reference = "mclust", covers = "multilpa()",
                   quantity = "noise-model log likelihood from mclust's start", tolerance = 1e-10,
                   measure = "relative", key = "code")
  # Both EMs stop on a relative log likelihood change of 1e-13, which fixes
  # the estimates only to about sqrt(1e-13) = 3e-7 on a flat maximum: on
  # Linux under R 4.6.2 the two stopped 1e-7 apart in the noise proportion
  # with log likelihoods equal to 1e-10. 1e-6 is the stopping rule's bound.
  expect_reference(actual[c("code", "noise_probability")], expected[c("code", "noise_probability")],
                   reference = "mclust", covers = "multilpa()",
                   quantity = "noise proportion from mclust's start", tolerance = 1e-6,
                   measure = "relative", key = "code")
  # Same stopping rule, same bound.
  expect_reference(stack("actual_means"), stack("expected_means"),
                   reference = "mclust", covers = "multilpa()",
                   quantity = "profile means of the noise model from mclust's start", tolerance = 1e-6,
                   measure = "relative", key = c("code", "variable", "profile"))
  # A count: exact.
  expect_reference(actual[c("code", "n_parameters")], expected[c("code", "n_parameters")],
                   reference = "mclust", covers = "multilpa()",
                   quantity = "number of free parameters of the noise model", tolerance = 0,
                   key = "code")
})

# Moved from tests/testthat/test-noise.R: comparisons with independent references.

test_that("the likelihood is the mixture with a uniform component", {
  data <- noise_frame()
  fit <- noise_fit(data)
  x <- as.matrix(data)
  volume <- fit$hypervolume
  densities <- cbind(vapply(seq_len(2L), function(k) {
    fit$profile_probabilities[1L, k] * exp(rowSums(stats::dnorm(
      x, rep(fit$means[k, ], each = nrow(x)),
      rep(sqrt(fit$variances[k, ]), each = nrow(x)), log = TRUE)))
  }, numeric(nrow(x))), fit$noise_probability / volume)
  # Gaussian densities from stats::dnorm() and the uniform density 1 / volume.
  expect_reference(fit$log_likelihood, sum(log(rowSums(densities))),
                   reference = "independent: Gaussian mixture with a uniform component, stats::dnorm()",
                   covers = "multilpa()", quantity = "log likelihood with a noise component",
                   tolerance = 1e-10, measure = "relative")
  # Posteriors over the Gaussian profiles and the noise sum to one, and so do
  # the proportions.
  expect_equal(unname(rowSums(fit$subject_posteriors)) + unname(fit$noise_posteriors),
               rep(1, nrow(x)), tolerance = 1e-12)
  expect_equal(sum(fit$profile_probabilities) + fit$noise_probability, 1,
               tolerance = 1e-12)
  # One proportion and the hypervolume: mclust's count.
  plain <- quietly(multilpa(data, c("a", "b", "c"), id = NULL, n_profiles = 2L,
                            n_starts = 1L, seed = 1L))
  expect_equal(fit$n_parameters, plain$n_parameters + 2)
  # The scattered points are the ones the component takes.
  expect_gte(mean(fit$subject_profiles[271:300] == 0L), 0.9)
})
