# Is the likelihood at the posterior mode monotone in the prior's strength?
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(mclust))
source("equivalence/argument-sweep/helpers.R")
source("equivalence/argument-sweep/datasets.R")
RNGkind("L'Ecuyer-CMRG")
datasets <- sweep_continuous()
settings <- list(list(dataset = "srl", k = 2L, model = "VVI"),
                 list(dataset = "blobs", k = 3L, model = "EEI"))
options(digits = 11, width = 200)
print(do.call(rbind, lapply(settings, function(s) {
  data <- datasets[[s$dataset]]
  fit <- function(...) suppressWarnings(lpa(data, names(data), n_profiles = s$k, model = s$model,
                                            seed = 11L, n_starts = 10L, tol = 1e-10,
                                            max_iter = 20000L, ...))
  base <- fit()
  shrinkages <- c(1, 0.1, 0.01, 1e-4, 1e-8)
  # mclust's own posterior mode under the same prior, from the package's maximum likelihood fit.
  rows <- lapply(shrinkages, function(shrinkage) {
    ours <- fit(prior = prior_control(shrinkage = shrinkage))
    theirs <- suppressWarnings(mclust::me(
      as.matrix(data), modelName = s$model, z = sweep_posteriors(ours),
      prior = mclust::priorControl(shrinkage = shrinkage),
      control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(20000L, 20000L))))
    data.frame(dataset = s$dataset, model = s$model, shrinkage = shrinkage,
               maximum_likelihood = base$log_likelihood, at_mode = ours$log_likelihood,
               gap = base$log_likelihood - ours$log_likelihood,
               mclust_continued = theirs$loglik, iterations = ours$iterations,
               replicated = ours$n_best_replicated)
  })
  do.call(rbind, rows)
})))
