# Which side holds the higher likelihood under VVE, and does either move the other?
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(mclust))
source("equivalence/argument-sweep/helpers.R")
source("equivalence/argument-sweep/datasets.R")
datasets <- sweep_continuous()
rows <- lapply(sweep_grid(dataset = names(datasets), k = 2:3, model = c("VVE", "EVE")), function(o) {
  data <- datasets[[o$dataset]]
  x <- as.matrix(data)
  fit <- lpa(data, names(data), n_profiles = o$k, model = o$model, seed = 11L, n_starts = 10L,
             tol = 1e-10, max_iter = 20000L)
  tight <- lpa(data, names(data), n_profiles = o$k, model = o$model, seed = 11L, n_starts = 10L,
               tol = 1e-13, max_iter = 200000L, acceleration = "none")
  continued <- suppressWarnings(mclust::me(x, modelName = o$model, z = sweep_posteriors(fit),
    control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(20000L, 20000L))))
  rebuilt <- sweep_mixture(data, sweep_parameters(fit))$log_likelihood
  # mclust's own density at mclust's continued parameters, and one E-step of ours from it.
  z <- continued$z
  back <- lpa(data, names(data), n_profiles = o$k, model = o$model, seed = 11L, n_starts = 1L,
              tol = 1e-13, max_iter = 200000L,
              start = list(means = t(continued$parameters$mean),
                           variances = t(apply(continued$parameters$variance$sigma, 3L, diag)),
                           covariances = unname(continued$parameters$variance$sigma),
                           profile_probabilities = matrix(continued$parameters$pro, 1L),
                           group_probabilities = 1))
  data.frame(dataset = o$dataset, k = o$k, model = o$model,
             latents = fit$log_likelihood, latents_iterations = fit$iterations,
             latents_converged = fit$converged, latents_tight = tight$log_likelihood,
             tight_iterations = tight$iterations, rebuilt = rebuilt,
             mclust_continued = continued$loglik,
             mclust_minus_latents = continued$loglik - fit$log_likelihood,
             latents_from_mclust = back$log_likelihood,
             from_mclust_minus_mclust = back$log_likelihood - continued$loglik)
})
options(width = 250, digits = 12)
print(do.call(rbind, rows))
