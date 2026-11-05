# Supplementary null calibration: common orientation is unidentified at sphericity.
# Strongly separated profiles keep the class-number issue out of this check.
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
iterations <- if (length(args)) as.integer(args[[1L]]) else 500L
set.seed(91471)
quiet <- function(expr) withCallingHandlers(expr,
  warning = function(w) invokeRestart("muffleWarning"),
  message = function(m) invokeRestart("muffleMessage"))
rows <- lapply(seq_len(iterations), function(i) {
  state <- rep(1:2, each = 500L)
  data <- data.frame(a = rnorm(1000L, c(-6, 6)[state]), b = rnorm(1000L))
  fit <- function(model) quiet(lpa(data, c("a", "b"), 2, model = model,
                                   n_starts = 1, seed = 1, tol = 1e-10))
  null <- fit("EII")
  alternative <- fit("EVE")
  tested <- as.data.frame(compare_models(null, alternative, test = "lrt"))
  data.frame(iteration = i, statistic = tested$lr_statistic[2L],
             df = tested$lr_df[2L], p_value = tested$lr_p_value[2L],
             boundary = null$boundary || alternative$boundary)
})
results <- do.call(rbind, rows)
utils::write.csv(results, "validation/review-2026-11-01/lrt-orientation.csv", row.names = FALSE)
cat("Datasets:", nrow(results), "\n")
cat("Boundary fits:", sum(results$boundary), "\n")
cat("Degrees of freedom:", unique(results$df), "\n")
cat("Mean statistic:", mean(results$statistic), "\n")
cat("5% rejection rate:", mean(results$p_value < 0.05), "\n")
print(stats::ks.test(results$p_value, "punif"))
