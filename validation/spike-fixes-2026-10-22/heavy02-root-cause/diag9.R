# The actual heavy_02 model (direct effect y3 ~ z estimated, data not adjusted):
# does a model-based hierarchical start on the raw indicators reach -3098.487?
suppressMessages({pkgload::load_all(".", quiet = TRUE); library(mclust)})
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
vars <- case$spec$vars
y <- as.matrix(subset(case$data, select = vars))
start_from <- function(partition) {
  groups <- sort(unique(partition))
  structure(list(
    means = t(vapply(groups, \(g) colMeans(y[partition == g, , drop = FALSE]), numeric(ncol(y)))),
    variances = t(vapply(groups, \(g) apply(y[partition == g, , drop = FALSE], 2, \(v) mean((v - mean(v))^2)),
                         numeric(ncol(y)))),
    profile_probabilities = matrix(tabulate(partition) / length(partition), 1L),
    group_probabilities = 1), class = "multilpa_start")
}
partitions <- list(
  ward = .multilpa_ward_assignments(scale(y), 2L),
  hc_VII_vars = as.vector(hclass(hc(y, modelName = "VII", use = "VARS"), 2)),
  hc_VVV_vars = as.vector(hclass(hc(y, modelName = "VVV", use = "VARS"), 2)),
  hc_VVV_svd = as.vector(hclass(hc(y, modelName = "VVV", use = "SVD"), 2)),
  hc_VII_std = as.vector(hclass(hc(y, modelName = "VII", use = "STD"), 2)),
  hc_VVV_std = as.vector(hclass(hc(y, modelName = "VVV", use = "STD"), 2)),
  hc_VII_svd = as.vector(hclass(hc(y, modelName = "VII", use = "SVD"), 2)),
  hc_EII_std = as.vector(hclass(hc(y, modelName = "EII", use = "STD"), 2)))
# Replace the first (hierarchical) start's partition in a scratch session only.
original <- .multilpa_ward_assignments
reached <- vapply(partitions, function(partition) {
  assignInNamespace(".multilpa_ward_assignments", function(input, n_clusters) partition, "latents")
  on.exit(assignInNamespace(".multilpa_ward_assignments", original, "latents"), add = TRUE)
  fit <- withCallingHandlers(
    lpa(case$data, vars, 2, model = "VVI", n_starts = 1, seed = 1,
        tol = 1e-10, max_iter = 5000, direct_effects = "y3 ~ z"),
    warning = function(w) invokeRestart("muffleWarning"))
  stopifnot(identical(get_results(fit, "starts")$start, 1L))
  fit$log_likelihood
}, numeric(1))
print(round(reached, 4))
