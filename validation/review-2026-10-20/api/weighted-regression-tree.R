library(latents, lib.loc = "validation/review-2026-10-20/rlib")
set.seed(318)
d <- data.frame(y = c(rnorm(160, -4, .3), rnorm(40, 4, .3)),
                w = c(rep(1, 160), rep(16, 40)))
grid <- suppressWarnings(enumerate_regressions(y ~ 1, d, n_classes = 1:2,
  weights = "w", n_starts = 5, seed = 27))
print(as.data.frame(grid))
print(subset(plot_views(grid), type == "tree"))
tree <- get_results(grid, "tree")
print(tree)
fit <- get_results(grid, "model", n_classes = 2)
p <- fit$expectation$tau
w <- fit$spec$row_weights
cat("Weighted fitted shares:", colSums(p * w) / sum(w), "\n")
cat("Unweighted fitted shares:", colMeans(p), "\n")
cat("Tree shares:", tree$to_share, "\n")
cat("Tree discrepancy:", max(abs(tree$to_share - colSums(p * w) / sum(w))), "\n")
stopifnot(max(abs(tree$to_share - colSums(p * w) / sum(w))) > .5)
