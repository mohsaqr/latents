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

# One row per cell of a d x d x K covariance array or a d x K mean matrix.
prior_cells <- function(code, values) {
  cells <- expand.grid(lapply(dim(values), seq_len))
  names(cells) <- c("row", "column", "profile")[seq_along(dim(values))]
  data.frame(code = code, cells, value = as.vector(values))
}

test_that("the MAP M-step matches mclust's EM step from the same posteriors", {
  skip_without_reference("mclust")
  skip_on_cran()
  data <- as.matrix(prior_frame())
  set.seed(9)
  z <- matrix(stats::runif(nrow(data) * 3L), nrow(data), 3L)
  z <- z / rowSums(z)
  steps <- lapply(c("EII", "EEI", "EVI", "VEV", "EEV", "VVV"), function(code) {
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
    list(actual_sigma = prior_cells(code, array(sigma, c(3L, 3L, 3L))),
         expected_sigma = prior_cells(code, reference$parameters$variance$sigma),
         actual_mean = prior_cells(code, t(ours$means)),
         expected_mean = prior_cells(code, reference$parameters$mean))
  })
  stack <- function(part) do.call(rbind, lapply(steps, `[[`, part))
  note <- "Package side: the internal MAP M-step (.multilpa_prior_maximize) of multilpa(prior = prior_control()), from fixed posteriors; reference: one iteration of mclust::me() with mclust::priorControl()."
  expect_reference(stack("actual_sigma"), stack("expected_sigma"), reference = "mclust",
                   covers = "multilpa()", quantity = "MAP covariance matrices, one M-step",
                   tolerance = 1e-10, measure = "relative",
                   key = c("code", "row", "column", "profile"), note = note)
  expect_reference(stack("actual_mean"), stack("expected_mean"), reference = "mclust",
                   covers = "multilpa()", quantity = "MAP means, one M-step",
                   tolerance = 1e-10, measure = "relative",
                   key = c("code", "row", "column"), note = note)
})
