# Do the higher maxima recover the generating classes better?
# heavy_02 (low vs high) and the seven correct-k VVV misses (kmeans_10 vs the
# heavy emEM reference search of the benchmark).
suppressMessages({pkgload::load_all(".", quiet = TRUE); library(mclust); library(simulab)})
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
source(file.path("equivalence", "search", "spike-measure.R"))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
fits <- readRDS("/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/heavy/fits.rds")
cl <- \(fit) max.col(fit$subject_posteriors)
cat("heavy_02 ARI with the generating classes: low", round(adjustedRandIndex(cl(fits$low), case$truth$membership), 4),
    " high", round(adjustedRandIndex(cl(fits$high), case$truth$membership), 4), "\n\n")
spec <- list(n = 300, means = rbind(rep(0, 5), c(1.2, 1.2, 1.2, 0, 0)), sds = rbind(rep(1, 5), rep(0.8, 5)))
replicates <- c(2, 6, 9, 10, 12, 16, 19)
rows <- lapply(replicates, function(replicate) {
  seed <- 20261021L + 1000L * 2L + replicate
  d <- simulate_lpa(n = spec$n, means = spec$means, sds = spec$sds, seed = seed)
  vars <- grep("^indicator_", names(d), value = TRUE)
  truth <- d$profile
  overall <- vapply(d[vars], \(v) mean((v - mean(v))^2), numeric(1))
  quiet <- \(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
  km <- quiet(lpa(d, vars, 2, model = "VVV", seed = seed, search = "kmeans"))
  ref <- quiet(lpa(d, vars, 2, model = "VVV", seed = seed, search = "emEM", n_starts = 50, n_screen = 1000, boundary = "exclude"))
  data.frame(replicate, kmeans_ll = km$log_likelihood, reference_ll = ref$log_likelihood,
             kmeans_ARI = adjustedRandIndex(cl(km), truth), reference_ARI = adjustedRandIndex(cl(ref), truth),
             kmeans_min_class = min(colSums(km$subject_posteriors)), reference_min_class = min(colSums(ref$subject_posteriors)),
             kmeans_ratio = spike_measure(km, overall), reference_ratio = spike_measure(ref, overall))
})
print(do.call(rbind, rows), digits = 4, row.names = FALSE)
