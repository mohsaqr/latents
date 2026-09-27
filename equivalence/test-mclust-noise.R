# multilpa(noise = TRUE) against mclust's noise component: `me(..., Vinv = )`,
# the EM `Mclust(..., initialization = list(noise = ))` runs. Both start from
# the same posterior matrix (Gaussian columns, then noise), so what is compared
# is two EM iterations, not two optimisers' luck.

noise_reference_data <- function() {
  set.seed(3)
  x <- rbind(matrix(stats::rnorm(135L * 3L, 2), 135L, 3L),
             matrix(stats::rnorm(135L * 3L, -2), 135L, 3L) %*%
               matrix(c(1, .5, 0, 0, 1, .3, 0, 0, 1), 3L),
             matrix(stats::runif(30L * 3L, -10, 10), 30L, 3L))
  colnames(x) <- c("a", "b", "c")
  x
}

noise_compare <- function(x, z, code, prior = FALSE) {
  inverse_volume <- mclust::hypvol(x, reciprocal = TRUE)
  run <- function(itmax, tol) {
    arguments <- list(modelName = code, data = x, z = z, Vinv = inverse_volume,
                      control = mclust::emControl(itmax = c(itmax, 1e6),
                                                  tol = c(tol, 1e-14)))
    if (prior) arguments$prior <- mclust::priorControl()
    do.call(mclust::me, arguments, envir = asNamespace("mclust"))
  }
  # One iteration of mclust's EM is its M-step from `z`; our start is that.
  first <- run(1, 1e-5)$parameters
  reference <- run(1e6, 1e-13)
  k <- ncol(z) - 1L
  start <- list(means = t(first$mean), profile_probabilities = matrix(first$pro[seq_len(k)], 1L),
                group_probabilities = 1, noise_probability = first$pro[k + 1L])
  if (endsWith(code, "I")) {
    start$variances <- t(vapply(seq_len(k), function(profile) {
      diag(first$variance$sigma[, , profile])
    }, numeric(ncol(x))))
  } else start$covariances <- first$variance$sigma
  fit <- suppressWarnings(do.call(multilpa, c(
    list(data = as.data.frame(x), vars = colnames(x), id = NULL, n_profiles = k,
         n_starts = 1L, start = start, tol = 1e-13, max_iter = 20000L,
         noise = TRUE, prior = if (prior) prior_control()),
    latents:::.multilpa_structure_arguments(code))))
  sigma <- if (is.null(fit$covariances)) {
    array(unlist(lapply(seq_len(k), function(profile) diag(fit$variances[profile, ]))),
          c(ncol(x), ncol(x), k))
  } else unname(fit$covariances)
  data.frame(
    code = code, prior = prior, converged = fit$converged,
    loglik = fit$log_likelihood - reference$loglik,
    means = max(abs(t(fit$means) - reference$parameters$mean)),
    covariances = max(abs(sigma - unname(reference$parameters$variance$sigma))),
    noise = abs(fit$noise_probability - reference$parameters$pro[k + 1L]),
    posteriors = max(abs(cbind(unname(fit$subject_posteriors), fit$noise_posteriors) -
                           reference$z)),
    parameters = fit$n_parameters - mclust::nMclustParams(code, ncol(x), k, noise = TRUE),
    bic = fit$bic + mclust::bic(code, reference$loglik, nrow(x), ncol(x), k,
                                noise = TRUE))
}

noise_start_z <- function(x) {
  set.seed(5)
  flagged <- c(rep(FALSE, 270L), rep(TRUE, 30L))
  flagged[sample(270L, 10L)] <- TRUE
  clusters <- stats::kmeans(x[!flagged, ], 2L)$cluster
  z <- matrix(0, nrow(x), 3L)
  z[cbind(which(!flagged), clusters)] <- 1
  z[flagged, 3L] <- 1
  z
}

test_that("the hypervolume is mclust's hypvol()", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  x <- noise_reference_data()
  expect_equal(exp(latents:::.multilpa_log_hypervolume(x)), mclust::hypvol(x),
               tolerance = 1e-12)
})

test_that("a noise fit reproduces mclust's EM for all fourteen structures", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  x <- noise_reference_data()
  z <- noise_start_z(x)
  table <- do.call(rbind, lapply(latents:::.multilpa_structures(), noise_compare,
                                 x = x, z = z))
  message(paste(utils::capture.output(print(table, digits = 3)), collapse = "\n"))
  expect_true(all(table$converged))
  expect_true(all(table$parameters == 0))
  # EVE and VVE maximise their shared orientation by an iteration in both
  # packages; there the requirement is a likelihood no worse than mclust's.
  exact <- !table$code %in% c("EVE", "VVE")
  expect_true(all(abs(table$loglik[exact]) < 1e-8))
  expect_true(all(unlist(table[exact, c("means", "covariances", "noise",
                                        "posteriors")]) < 1e-6))
  expect_true(all(abs(table$bic[exact]) < 1e-7))
  expect_true(all(table$loglik[!exact] > -1e-8))
})

test_that("noise and prior together reproduce mclust's EM", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  x <- noise_reference_data()
  z <- noise_start_z(x)
  table <- do.call(rbind, lapply(latents:::.multilpa_prior_structures(),
                                 noise_compare, x = x, z = z, prior = TRUE))
  message(paste(utils::capture.output(print(table, digits = 3)), collapse = "\n"))
  expect_true(all(abs(table$loglik) < 1e-8))
  expect_true(all(unlist(table[, c("means", "covariances", "noise", "posteriors")]) < 1e-6))
})

test_that("from Mclust's own noise initialization the fits agree", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  mclustBIC <- mclust::mclustBIC
  x <- noise_reference_data()
  flagged <- c(rep(FALSE, 270L), rep(TRUE, 30L))
  reference <- mclust::Mclust(x, G = 2L, modelNames = "VVV", verbose = FALSE,
                              initialization = list(noise = flagged))
  # Mclust's own final posteriors, as the common start of both EMs.
  row <- noise_compare(x, reference$z, "VVV")
  message(paste(utils::capture.output(print(row, digits = 3)), collapse = "\n"))
  expect_lt(abs(row$loglik), 1e-8)
  expect_lt(row$covariances, 1e-6)
  expect_equal(reference$hypvol, exp(latents:::.multilpa_log_hypervolume(x)),
               tolerance = 1e-12)
})
