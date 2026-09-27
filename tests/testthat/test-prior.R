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

test_that("prior_control() carries mclust's defaults", {
  prior <- prior_control()
  expect_s3_class(prior, "latents_prior")
  expect_identical(prior$shrinkage, 0.01)
  expect_null(prior$mean)
  expect_null(prior$dof)
  expect_null(prior$scale)
  expect_output(print(prior), "Fraley")
  # The resolved hyperparameters are defaultPrior()'s formulas.
  data <- as.matrix(prior_frame())
  diagonal <- latents:::.multilpa_prior_parameters(prior, data, 3L, "VVI")
  expect_equal(diagonal$scale, (1 / 3)^(2 / 3) * sum(apply(data, 2L, stats::var)) / 3)
  expect_identical(diagonal$dof, 5)
  full <- latents:::.multilpa_prior_parameters(prior, data, 3L, "VVV")
  expect_equal(full$scale, unname((1 / 3)^(2 / 3) * stats::var(data)))
  expect_equal(full$mean, unname(colMeans(data)))
})

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

test_that("a prior the package does not define is refused by class", {
  data <- prior_frame()
  expect_error(prior_fit(data, "VEE"), class = "latents_unsupported_prior")
  expect_error(prior_fit(data, "EVV"), class = "latents_unsupported_prior")
  with_missing <- data
  with_missing$a[3] <- NA
  expect_error(prior_fit(with_missing, "VVV", missing = "fiml"),
               class = "latents_unsupported_prior")
  expect_error(prior_control(shrinkage = -1), class = "latents_bad_argument")
  expect_error(prior_control(shrinkage = 0, mean = c(0, 0, 0)),
               class = "latents_bad_argument")
  expect_error(prior_fit(data, "VVV", prior = list(shrinkage = 0.01)),
               class = "latents_bad_argument")
  # Wald errors invert the likelihood's curvature at its maximum; a posterior
  # mode is not one.
  fit <- prior_fit(data, "VVI")
  expect_error(parameter_inference(fit), class = "latents_unsupported_inference")
})

test_that("a MAP fit is translation equivariant and reports its convention", {
  data <- prior_frame()
  shifted <- data
  shifted[c("a", "b", "c")] <- sweep(as.matrix(data), 2L, c(10, -5, 3), "+")
  original <- prior_fit(data, "VEV")
  moved <- prior_fit(shifted, "VEV")
  # The default prior mean is the data's mean, so it moves with the data.
  expect_equal(unname(moved$means), unname(sweep(original$means, 2L, c(10, -5, 3), "+")),
               tolerance = 1e-8)
  expect_equal(unname(moved$covariances), unname(original$covariances),
               tolerance = 1e-8)
  expect_equal(moved$log_likelihood, original$log_likelihood, tolerance = 1e-10)
  # The criteria are the unpenalized likelihood's, as in mclust.
  expect_equal(original$bic, -2 * original$log_likelihood +
                 log(nrow(data)) * original$n_parameters)
  expect_identical(original$prior, prior_control())
  expect_equal(original$prior_parameters$mean, unname(colMeans(data)))
  expect_output(print(original), "conjugate prior")
})

test_that("shrinkage pulls the means towards the prior mean", {
  data <- prior_frame()
  loose <- prior_fit(data, "VVI", prior = prior_control(shrinkage = 0.01))
  tight <- prior_fit(data, "VVI", prior = prior_control(shrinkage = 500))
  centre <- colMeans(as.matrix(data))
  spread <- function(fit) sum(sweep(fit$means, 2L, centre, "-")^2)
  expect_lt(spread(tight), spread(loose))
})
