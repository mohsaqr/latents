# With a well-posed objective (the existing Fraley-Raftery prior as a stand-in),
# do kmeans and emEM still disagree, and what do they recover?
suppressMessages({pkgload::load_all(".", quiet = TRUE); library(mclust); library(simulab)})
source(file.path("equivalence", "search", "spike-measure.R"))
cl <- \(fit) max.col(fit$subject_posteriors)
quiet <- \(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
spec <- list(n = 300, means = rbind(rep(0, 5), c(1.2, 1.2, 1.2, 0, 0)), sds = rbind(rep(1, 5), rep(0.8, 5)))
rows <- lapply(c(2, 6, 9, 10, 12, 16, 19), function(replicate) {
  seed <- 20261021L + 1000L * 2L + replicate
  d <- simulate_lpa(n = spec$n, means = spec$means, sds = spec$sds, seed = seed)
  vars <- grep("^indicator_", names(d), value = TRUE)
  overall <- vapply(d[vars], \(v) mean((v - mean(v))^2), numeric(1))
  one <- function(search, prior) {
    fit <- quiet(lpa(d, vars, 2, model = "VVV", seed = seed, search = search, prior = prior))
    c(ARI = adjustedRandIndex(cl(fit), d$profile), min_class = min(colSums(fit$subject_posteriors)),
      ratio = spike_measure(fit, overall))
  }
  out <- rbind(kmeans = one("kmeans", NULL), emEM = one("emEM", NULL),
               kmeans_prior = one("kmeans", prior_control()), emEM_prior = one("emEM", prior_control()))
  data.frame(replicate, fit = rownames(out), out, row.names = NULL)
})
all_rows <- do.call(rbind, rows)
print(all_rows, digits = 3, row.names = FALSE)
print(aggregate(cbind(ARI, min_class) ~ fit, all_rows, median))
