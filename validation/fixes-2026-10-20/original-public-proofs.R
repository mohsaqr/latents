library(latents)
set.seed(17)
data <- data.frame(g = rep(1:5, each = 30), y = stats::rnorm(150))
padded <- rbind(data, data.frame(g = 6:100, y = NA_real_))
fit <- function(data) multilpa(data, "y", "g", n_profiles = 1,
  n_group_classes = 1, missing = "fiml", n_starts = 1, seed = 3)
original <- parameter_inference(fit(data), vcov_type = "robust")
expanded <- parameter_inference(fit(padded), vcov_type = "robust")
cat("Original review's five-group example after repair:\n")
print(data.frame(parameter = original$parameter[1:2],
  original_se = original$standard_error[1:2],
  padded_se = expanded$standard_error[1:2]))
stopifnot(isTRUE(all.equal(original$standard_error, expanded$standard_error,
  tolerance = 1e-10)))

set.seed(318)
data <- data.frame(y = c(stats::rnorm(160, -4, .3), stats::rnorm(40, 4, .3)),
                   w = c(rep(1, 160), rep(16, 40)))
grid <- enumerate_regressions(y ~ 1, data, n_classes = 1:2,
  weights = "w", n_starts = 5, seed = 27)
tree <- get_results(grid, "tree")
fitted <- get_results(grid, "model", n_classes = 2)
expected <- colSums(fitted$expectation$tau * fitted$spec$row_weights) /
  sum(fitted$spec$row_weights)
cat("\nOriginal review's weighted tree after repair:\n")
print(data.frame(class = 1:2, weighted_fit = expected,
  tree_share = tree$to_share, tree_mass = tree$mass))
stopifnot(max(abs(tree$to_share - expected)) < 1e-12,
          max(abs(tree$mass - expected)) < 1e-12)

data <- expand.grid(a = 0:1, b = 0:1, c = 0:1)
data[["a:b"]] <- seq_len(nrow(data))
collision <- tryCatch(lca(data, names(data), 1, local_dependence = "a:b",
  n_starts = 1, seed = 1), error = identity)
cat("\nOriginal pair-name collision is refused:\n")
print(collision)
stopifnot(inherits(collision, "latents_bad_data"))
ambiguous <- data.frame(a = rep(c(.3, .1 + .2), 100),
  b = rep(c(0, 0, 1, 1), 50), c = rep(0:1, each = 100))
labels <- tryCatch(lca(ambiguous, names(ambiguous), 1,
  local_dependence = "a:b", n_starts = 1), error = identity)
cat("\nOriginal ambiguous native categories are refused:\n")
print(labels)
stopifnot(inherits(labels, "latents_bad_categorical"))
