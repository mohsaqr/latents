library(latents, lib.loc = "validation/review-2026-10-20/rlib")
set.seed(1020)
data <- data.frame(g = rep(seq_len(20), each = 5), t = rep(0:4, 20))
data$y <- rep(stats::rnorm(20), each = 5) + 0.4 * data$t + stats::rnorm(100, sd = .7)
data$weight <- rep(rep(c(1, 2), 10), each = 5)
added <- data[data$g <= 5, ]
added$g <- added$g + 20
added$weight <- 0
padded <- rbind(data, added)
fitters <- list(
  additive = function(d) multilpa(d, "y", "g", n_group_classes = 1,
    family = "additive", weights = "weight", n_starts = 1, seed = 1, tol = 1e-10),
  regression = function(d) mixture_regression(y ~ t, d, n_classes = 1, id = "g",
    weights = "weight", n_starts = 1, seed = 1, tol = 1e-10),
  growth = function(d) mixture_regression(y ~ t, d, n_classes = 1, id = "g",
    class_level = "group", random = ~ 1, weights = "weight", n_starts = 1,
    seed = 1, tol = 1e-10))
invisible(lapply(names(fitters), function(name) {
  a <- fitters[[name]](data)
  b <- fitters[[name]](padded)
  cat("\n", name, "\n")
  str(a, max.level = 1)
  covariance_a <- vcov(a, vcov_type = "robust")
  covariance_b <- vcov(b, vcov_type = "robust")
  print(cbind(diag(covariance_a), diag(covariance_b),
    ratio = diag(covariance_b) / diag(covariance_a)))
  cat("Expected old CR1 ratio:", (25 / 24) / (20 / 19), "\n")
}))
