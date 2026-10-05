# Maximum a posteriori estimation under mclust's conjugate prior
# (`multilpa(prior = prior_control())`). The full ten-structure comparison
# with mclust lives in equivalence/test-mclust-prior.R.

prior_frame <- function(seed = 3L) {
  set.seed(seed)
  n <- 150L
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  x[1:50, ] <- x[1:50, ] %*% diag(c(3, 1, 0.5)) + 2
  x[51:100, ] <- x[51:100, ] %*% matrix(c(1, .5, 0, 0, 2, .3, 0, 0, 1), 3L) - 1
  colnames(x) <- c("a", "b", "c")
  as.data.frame(x)
}

prior_fit <- function(data, code, prior = prior_control(), ...) {
  quietly(do.call(multilpa, c(
    list(data = data, vars = c("a", "b", "c"), id = NULL, n_profiles = 3L,
         n_starts = 2L, seed = 1L, tol = 1e-12, max_iter = 5000L,
         prior = prior),
    latents:::.multilpa_structure_arguments(code), list(...))))
}

test_that("the MAP M-step matches mclust's EM step from the same posteriors", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  data <- as.matrix(prior_frame())
  set.seed(9)
  z <- matrix(stats::runif(nrow(data) * 3L), nrow(data), 3L)
  z <- z / rowSums(z)
  invisible(lapply(c("EII", "EEI", "EVI", "VEV", "EEV", "VVV"), function(code) {
    # One iteration of mclust's own EM is its M-step from `z`.
    reference <- do.call(mclust::me, list(
      modelName = code, data = data, z = z, prior = mclust::priorControl(),
      control = mclust::emControl(itmax = c(1, 1e6), tol = c(1e-5, 1e-14))),
      envir = asNamespace("mclust"))
    hyper <- latents:::.multilpa_prior_parameters(prior_control(), data, 3L, code)
    ours <- latents:::.multilpa_prior_maximize(data, z, code, hyper, 1e-12,
                                               tol = 1e-14)
    sigma <- ours$covariances %||% array(unlist(lapply(seq_len(3L), function(k) {
      diag(ours$variances[k, ])
    })), c(3L, 3L, 3L))
    expect_equal(sigma, reference$parameters$variance$sigma, tolerance = 1e-10,
                 ignore_attr = TRUE, info = code)
    expect_equal(t(ours$means), reference$parameters$mean, tolerance = 1e-10,
                 ignore_attr = TRUE, info = code)
  }))
})
