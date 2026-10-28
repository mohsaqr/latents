# heavy_02 (direct effect removed at the high solution's b): which starting
# partitions lead latents' own EM to which maximum, and how they differ.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(mclust))
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
vars <- case$spec$vars
adjusted <- transform(case$data, y3 = y3 - 0.8063057 * z)
y <- as.matrix(subset(adjusted, select = vars))
standardized <- scale(y)
set.seed(1)
partitions <- list(
  kmeans = stats::kmeans(standardized, 2, iter.max = 100, algorithm = "Lloyd")$cluster,
  ward = .multilpa_ward_assignments(standardized, 2L),
  hc_VVV_vars = as.vector(hclass(hc(y, modelName = "VVV", use = "VARS"), 2)),
  hc_VVV_svd = as.vector(hclass(hc(y, modelName = "VVV", use = "SVD"), 2)),
  hc_EII_vars = as.vector(hclass(hc(y, modelName = "EII", use = "VARS"), 2)),
  hc_VII_vars = as.vector(hclass(hc(y, modelName = "VII", use = "VARS"), 2)))
start_from <- function(partition) {
  groups <- sort(unique(partition))
  means <- t(vapply(groups, \(g) colMeans(y[partition == g, , drop = FALSE]), numeric(ncol(y))))
  variances <- t(vapply(groups, \(g) apply(y[partition == g, , drop = FALSE], 2, \(v) mean((v - mean(v))^2)),
                        numeric(ncol(y))))
  structure(list(means = means, variances = variances,
                 profile_probabilities = matrix(tabulate(partition) / length(partition), 1L),
                 group_probabilities = 1), class = "multilpa_start")
}
reached <- vapply(partitions, function(partition) {
  suppressWarnings(lpa(adjusted, vars, 2, model = "VVI", start = start_from(partition),
                       tol = 1e-10, max_iter = 5000))$log_likelihood
}, numeric(1))
sizes <- vapply(partitions, \(p) paste(tabulate(p), collapse = "/"), character(1))
agreement_with_kmeans <- vapply(partitions, \(p) adjustedRandIndex(p, partitions$kmeans), numeric(1))
print(data.frame(start = names(partitions), sizes, ARI_vs_kmeans = round(agreement_with_kmeans, 3),
                 log_likelihood = round(reached, 4)), row.names = FALSE)
high <- suppressWarnings(lpa(adjusted, vars, 2, model = "VVI", n_starts = 10, seed = 1, search = "emEM"))
low <- suppressWarnings(lpa(adjusted, vars, 2, model = "VVI", n_starts = 10, seed = 1, search = "kmeans"))
cl <- \(fit) max.col(fit$subject_posteriors)
cat("\nfinal classifications: ARI(low, high) =", round(adjustedRandIndex(cl(low), cl(high)), 3),
    "; differing rows =", sum(cl(low) != cl(high)) , "or", sum(cl(low) == cl(high)), "\n")
cat("ARI of each start with final low / high:\n")
print(round(t(vapply(partitions, \(p) c(low = adjustedRandIndex(p, cl(low)), high = adjustedRandIndex(p, cl(high))),
                     numeric(2))), 3))
