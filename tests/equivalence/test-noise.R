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
  # The axis-aligned box, by hand, is an upper bound on it.
  box <- sum(log(apply(data, 2L, function(column) diff(range(column)))))
  expect_lte(latents:::.multilpa_log_hypervolume(data), box + 1e-12)
  skip_if_not_installed("mclust")
  expect_equal(exp(latents:::.multilpa_log_hypervolume(data)),
               mclust::hypvol(data), tolerance = 1e-12)
})

test_that("a noise fit reproduces mclust's EM from the same start", {
  skip_if_not_installed("mclust")
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
  invisible(lapply(c("VVI", "VEV"), function(code) {
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
    expect_equal(fit$log_likelihood, reference$loglik, tolerance = 1e-10, info = code)
    expect_equal(fit$noise_probability, reference$parameters$pro[3L],
                 tolerance = 1e-8, info = code)
    expect_equal(t(fit$means), reference$parameters$mean, tolerance = 1e-8,
                 ignore_attr = TRUE, info = code)
    expect_equal(fit$n_parameters, mclust::nMclustParams(code, 3L, 2L, noise = TRUE),
                 info = code)
  }))
})
