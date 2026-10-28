.libPaths(c("/private/tmp/latents-review-20261022/rlib", .libPaths()))
library(latents)
set.seed(20261022)
cases <- lapply(seq_len(120), function(i) {
  k <- sample(2:4, 1)
  d <- sample(2:4, 1)
  weights <- 10^runif(k, -2, 2)
  scatter <- matrix(10^runif(k * d, -7, 4), k, d) * weights
  if (i %% 4 == 0) scatter[sample.int(length(scatter), 1)] <- 0
  list(scatter = scatter, weights = weights,
       structure = sample(c("EVI", "VEI", "EEI"), 1), min_variance = 1)
})
results <- lapply(seq_along(cases), function(i) {
  value <- tryCatch(do.call(latents:::.multilpa_bounded_variances, cases[[i]]),
                    error = identity)
  list(index = i, value = value, case = cases[[i]])
})
saveRDS(results, "/private/tmp/latents-review-20261022/probe-covariance.rds")
failed <- Filter(function(result) inherits(result$value, "error"), results)
cat("CASES:", length(results), "FAILURES:", length(failed), "\n")
invisible(lapply(failed, function(result) {
  print(result$index)
  print(result$case)
  print(result$value)
}))
