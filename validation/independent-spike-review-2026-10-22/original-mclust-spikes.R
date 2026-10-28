# mclust on the benchmark's heteroscedastic datasets (same seeds as pilot.R),
# spike measure as in equivalence/search/spike-measure.R.
suppressPackageStartupMessages({ library(simulab); library(mclust); library(latents) })
source("/Users/mohammedsaqr/Documents/Github/latents/equivalence/search/spike-measure.R")
S <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/c17a2edf-450c-4631-b62e-18ced1b5c86e/scratchpad"
conditions <- list(
  lpa_vvi_k3 = list(index = 1L, model = "VVI", n = 300, fit_k = 3:4,
                    means = rbind(c(0, 0, 0, 0), c(1.5, 1.5, 0, 0), c(0, 1.5, 1.5, 1.5)),
                    sds = rbind(rep(1, 4), rep(0.7, 4), rep(1.3, 4))),
  lpa_vvv_k2 = list(index = 2L, model = "VVV", n = 300, fit_k = 2:3,
                    means = rbind(rep(0, 5), c(1.2, 1.2, 1.2, 0, 0)),
                    sds = rbind(rep(1, 5), rep(0.8, 5))))
one <- function(name, replicate, k, prior) {
  spec <- conditions[[name]]
  seed <- 20261021L + 1000L * spec$index + replicate
  d <- simulate_lpa(n = spec$n, means = spec$means, sds = spec$sds, seed = seed)
  d <- d[grep("^indicator_", names(d))]
  overall <- vapply(d, function(v) mean((v - mean(v))^2), numeric(1))
  fit <- Mclust(d, G = k, modelNames = spec$model, verbose = FALSE,
                prior = if (prior) priorControl())
  if (is.null(fit)) return(data.frame(condition = name, replicate = replicate, k = k,
    prior = prior, mclust_ll = NA_real_, ratio = NA_real_, min_class = NA_real_))
  data.frame(condition = name, replicate = replicate, k = k, prior = prior,
             mclust_ll = fit$loglik,
             ratio = spike_measure(list(covariances = fit$parameters$variance$sigma), overall),
             min_class = min(colSums(fit$z)))
}
grid <- expand.grid(replicate = 1:20, k = 1:2, condition = names(conditions), prior = c(FALSE, TRUE),
                    stringsAsFactors = FALSE)
rows <- do.call(rbind, Map(function(cn, r, ki, p) one(cn, r, conditions[[cn]]$fit_k[ki], p),
                           grid$condition, grid$replicate, grid$k, grid$prior))
# Data check: latents kmeans_10 on replicate 1 reproduces the pilot's likelihood.
check <- vapply(names(conditions), function(name) {
  spec <- conditions[[name]]; seed <- 20261021L + 1000L * spec$index + 1L
  d <- simulate_lpa(n = spec$n, means = spec$means, sds = spec$sds, seed = seed)
  d <- d[grep("^indicator_", names(d))]
  f <- suppressWarnings(lpa(d, names(d), spec$fit_k[1], model = spec$model, seed = seed, search = "kmeans"))
  pilot <- readRDS(file.path(S, "pilot", sprintf("%s_001.rds", name)))$rows
  f$log_likelihood - subset(pilot, fit == "kmeans_10" & k == spec$fit_k[1])$log_likelihood
}, numeric(1))
cat("data check (latents now minus pilot, replicate 1):\n"); print(check)
saveRDS(rows, file.path(S, "mclust", "rows.rds"))
