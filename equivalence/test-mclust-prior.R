# multilpa(prior = prior_control()) against mclust's conjugate-prior EM,
# `mclust::me(..., prior = priorControl())`, which is the engine `Mclust()`
# runs. Both start from the same posterior matrix, so the comparison is of two
# EM iterations, not of two optimisers' luck.

prior_reference_data <- function() {
  set.seed(3)
  n <- 300L
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  x[1:100, ] <- x[1:100, ] %*% diag(c(3, 1, 0.5)) + 2
  x[101:200, ] <- x[101:200, ] %*% matrix(c(1, .5, 0, 0, 2, .3, 0, 0, 1), 3L) - 1
  colnames(x) <- c("a", "b", "c")
  x
}

as_sigma <- function(fit_or_step, d, k) {
  if (!is.null(fit_or_step$covariances)) return(unname(fit_or_step$covariances))
  array(unlist(lapply(seq_len(k), function(profile) {
    diag(unname(fit_or_step$variances[profile, ]), d)
  })), c(d, d, k))
}

test_that("each MAP M-step is mclust's EM step from the same posteriors", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  x <- prior_reference_data()
  set.seed(9)
  z <- matrix(stats::runif(nrow(x) * 3L), nrow(x), 3L)
  z <- z / rowSums(z)
  worst <- vapply(latents:::.multilpa_prior_structures(), function(code) {
    reference <- do.call(mclust::me, list(
      modelName = code, data = x, z = z, prior = mclust::priorControl(),
      control = mclust::emControl(itmax = c(1, 1e6), tol = c(1e-5, 1e-14))),
      envir = asNamespace("mclust"))
    hyper <- latents:::.multilpa_prior_parameters(prior_control(), x, 3L, code)
    ours <- latents:::.multilpa_prior_maximize(x, z, code, hyper, 1e-12, tol = 1e-14)
    difference <- max(abs(as_sigma(ours, 3L, 3L) -
                            unname(reference$parameters$variance$sigma)))
    expect_lt(difference, 1e-8, label = code)
    expect_equal(t(ours$means), unname(reference$parameters$mean),
                 tolerance = 1e-10, ignore_attr = TRUE, info = code)
    difference
  }, numeric(1))
  message("max |covariance difference| of one M-step, by structure:\n",
          paste(sprintf("%s %.2e", names(worst), worst), collapse = "\n"))
})

test_that("the whole MAP EM reproduces mclust::me from the same start", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  x <- prior_reference_data()
  frame <- as.data.frame(x)
  set.seed(4)
  z <- mclust::unmap(stats::kmeans(scale(x), 3L, nstart = 5L)$cluster)
  z <- 0.8 * z + 0.2 / 3
  rows <- lapply(latents:::.multilpa_prior_structures(), function(code) {
    reference <- do.call(mclust::me, list(
      modelName = code, data = x, z = z, prior = mclust::priorControl(),
      control = mclust::emControl(tol = c(1e-13, 1e-14))),
      envir = asNamespace("mclust"))
    # mclust's `me()` starts with an M-step from `z`; so does this start.
    hyper <- latents:::.multilpa_prior_parameters(prior_control(), x, 3L, code)
    first <- latents:::.multilpa_prior_maximize(x, z, code, hyper, 1e-12, tol = 1e-14)
    start <- list(means = first$means,
                  profile_probabilities = matrix(colMeans(z), 1L),
                  group_probabilities = 1)
    if (is.null(first$covariances)) start$variances <- first$variances else
      start$covariances <- first$covariances
    fit <- suppressWarnings(do.call(multilpa, c(
      list(data = frame, vars = c("a", "b", "c"), id = NULL, n_profiles = 3L,
           n_starts = 1L, start = start, tol = 1e-13, max_iter = 20000L,
           prior = prior_control()),
      latents:::.multilpa_structure_arguments(code))))
    expect_true(fit$converged, label = code)
    expect_equal(fit$log_likelihood, reference$loglik, tolerance = 1e-10, info = code)
    # mclust's BIC is 2 logL - df log n on the unpenalized likelihood; ours is
    # its negative, with the same parameter count.
    expect_equal(fit$bic, -mclust::bic(code, reference$loglik, nrow(x), 3L, 3L),
                 tolerance = 1e-10, info = code)
    data.frame(
      code = code,
      loglik = abs(fit$log_likelihood - reference$loglik),
      means = max(abs(t(fit$means) - reference$parameters$mean)),
      covariances = max(abs(as_sigma(fit, 3L, 3L) -
                              unname(reference$parameters$variance$sigma))),
      proportions = max(abs(fit$profile_probabilities - reference$parameters$pro)),
      posteriors = max(abs(unname(fit$subject_posteriors) - reference$z)))
  })
  table <- do.call(rbind, rows)
  expect_true(all(unlist(table[, -1L]) < 1e-8))
  message(paste(utils::capture.output(print(table, digits = 3)), collapse = "\n"))
})
