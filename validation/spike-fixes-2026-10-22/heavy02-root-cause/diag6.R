# On the archived heteroscedastic benchmark (lpa_vvi_k3, lpa_vvv_k2; 40 datasets,
# 80 fits), one start from mclust's model-based hierarchical clustering (VVV
# criterion, SVD, mclust's default), run by latents' own EM. Compared with the
# recorded kmeans_10 and emem_10 results and the recorded best non-spike maximum.
suppressMessages({pkgload::load_all(".", quiet = TRUE); library(mclust); library(simulab)})
source(file.path("equivalence", "search", "spike-measure.R"))
recorded <- read.csv(file.path("validation", "review-fit-search-2026-10-22", "recomputed.csv"))
recorded <- subset(recorded, condition %in% c("lpa_vvi_k3", "lpa_vvv_k2"))
specs <- list(
  lpa_vvi_k3 = list(model = "VVI", n = 300, means = rbind(c(0, 0, 0, 0), c(1.5, 1.5, 0, 0), c(0, 1.5, 1.5, 1.5)),
                    sds = rbind(rep(1, 4), rep(0.7, 4), rep(1.3, 4))),
  lpa_vvv_k2 = list(model = "VVV", n = 300, means = rbind(rep(0, 5), c(1.2, 1.2, 1.2, 0, 0)),
                    sds = rbind(rep(1, 5), rep(0.8, 5))))
start_from <- function(y, partition, full) {
  k <- max(partition); p <- ncol(y)
  overall <- stats::cov(y) * (nrow(y) - 1) / nrow(y)
  spread <- lapply(seq_len(k), function(g) {
    rows <- partition == g
    if (sum(rows) <= p + 1L) return(overall)
    stats::cov(y[rows, , drop = FALSE]) * (sum(rows) - 1) / sum(rows)
  })
  means <- t(vapply(seq_len(k), \(g) colMeans(y[partition == g, , drop = FALSE]), numeric(p)))
  start <- list(means = means, variances = t(vapply(spread, diag, numeric(p))),
                profile_probabilities = matrix(tabulate(partition, k) / nrow(y), 1L),
                group_probabilities = 1)
  if (full) start$covariances <- array(unlist(spread), c(p, p, k))
  structure(start, class = "multilpa_start")
}
keys <- unique(subset(recorded, select = c(condition, replicate, seed, k, best)))
started <- proc.time()[["elapsed"]]
hc_rows <- lapply(seq_len(nrow(keys)), function(i) {
  key <- keys[i, ]
  spec <- specs[[key$condition]]
  d <- simulate_lpa(n = spec$n, means = spec$means, sds = spec$sds, seed = key$seed)
  vars <- grep("^indicator_", names(d), value = TRUE)
  y <- as.matrix(subset(d, select = vars))
  partition <- as.vector(hclass(hc(y, modelName = "VVV", use = "SVD"), key$k))
  fit <- withCallingHandlers(lpa(d, vars, key$k, model = spec$model, n_starts = 1, seed = key$seed,
                                 start = start_from(y, partition, spec$model == "VVV")),
                             warning = function(w) invokeRestart("muffleWarning"))
  overall <- apply(y, 2, \(v) mean((v - mean(v))^2))
  data.frame(key, hc_ll = fit$log_likelihood, hc_ratio = spike_measure(fit, overall),
             hc_min_class = min(colSums(fit$subject_posteriors)), hc_start_smallest = min(tabulate(partition, key$k)))
})
hc_rows <- do.call(rbind, hc_rows)
cat("seconds:", round(proc.time()[["elapsed"]] - started, 1), "\n")
spike <- function(r) !is.na(r) & r < .01
wide <- reshape(subset(recorded, select = c(condition, replicate, k, fit, log_likelihood, min_relative_variance)),
                idvar = c("condition", "replicate", "k"), timevar = "fit", direction = "wide")
all_rows <- merge(hc_rows, wide, by = c("condition", "replicate", "k"))
# A non-spike hc solution above the recorded reference raises the reference.
all_rows$reference <- pmax(all_rows$best, ifelse(spike(all_rows$hc_ratio), -Inf, all_rows$hc_ll))
status <- function(ll, ratio, reference) ifelse(spike(ratio), "spike",
  ifelse(reference - ll <= 1e-6 * (1 + abs(reference)), "best", "miss"))
pool_takes_hc <- all_rows$hc_ll > all_rows$log_likelihood.kmeans_10 + 1e-6
all_rows <- within(all_rows, {
  kmeans_10 <- status(log_likelihood.kmeans_10, min_relative_variance.kmeans_10, reference)
  emem_10 <- status(log_likelihood.emem_10, min_relative_variance.emem_10, reference)
  hc_only <- status(hc_ll, hc_ratio, reference)
  kmeans_10_plus_hc <- status(pmax(log_likelihood.kmeans_10, hc_ll),
                              ifelse(pool_takes_hc, hc_ratio, min_relative_variance.kmeans_10), reference)
})
long <- reshape(subset(all_rows, select = c(condition, replicate, k, kmeans_10, emem_10, hc_only, kmeans_10_plus_hc)),
                direction = "long", varying = c("kmeans_10", "emem_10", "hc_only", "kmeans_10_plus_hc"),
                v.names = "status", timevar = "search", times = c("kmeans_10", "emem_10", "hc_only", "kmeans_10_plus_hc"),
                idvar = c("condition", "replicate", "k"))
print(with(long, table(search, status)))
print(with(long, table(paste(condition, k), search, status)))
cat("\nreference raised by a non-spike hc maximum in", sum(all_rows$reference > all_rows$best + 1e-6), "fits\n")
print(subset(all_rows, hc_only == "spike" | (kmeans_10_plus_hc != kmeans_10),
             select = c(condition, replicate, k, best, hc_ll, hc_ratio, hc_min_class, log_likelihood.kmeans_10,
                        kmeans_10, kmeans_10_plus_hc)), row.names = FALSE)
saveRDS(all_rows, "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/heavy/hc-benchmark.rds")
