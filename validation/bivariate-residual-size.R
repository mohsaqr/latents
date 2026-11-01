# Size of the bivariate-residual tests behind the local-independence check of
# diagnostics(), under local independence. Run from the package root:
#   Rscript validation/bivariate-residual-size.R        # 200 per cell, 2 cores
#   Rscript validation/bivariate-residual-size.R 5      # a quick pilot
#
# Eight cells: lca() with five binary items and lpa() with four continuous
# indicators (VVI); K = 2 and 3 classes; n = 300 and 1000. Data are drawn from
# a latent class or profile model with local independence (well separated
# classes), refitted with the true K, and get_results(fit, "residuals") read.
# Reported per cell: the share of single tests with p < 0.05 (should be about
# 0.05), and the share of data sets with a flagged pair: any BH-adjusted
# p < 0.05, the criterion diagnostics() states (should be at most about 0.05:
# under the global null the false discovery rate is the family-wise rate),
# and any unadjusted p < 0.05, what diagnostics() computed up to 1.1.9. Declared before
# running: a test is taken as calibrated in a cell when its per-test rate lies
# within 0.025-0.075 (200 data sets of 10-30 tests each); the flag rate is
# reported against 0.05. If the tests are not calibrated, the check reports
# the residuals without a decision (as Latent GOLD does), which is the owner's
# decision to make.
pkgload::load_all(".", quiet = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(arguments) > 0L) as.integer(arguments[[1L]]) else 200L
cells <- expand.grid(kind = c("lca", "lpa"), classes = 2:3, n = c(300L, 1000L),
                     stringsAsFactors = FALSE)

replicate_cell <- function(replication, cell) {
  k <- cell$classes
  truth <- sample.int(k, cell$n, replace = TRUE)
  quiet <- function(expr) {
    withCallingHandlers(expr, warning = function(condition) invokeRestart("muffleWarning"))
  }
  fit <- if (identical(cell$kind, "lca")) {
    probability <- rbind(c(0.85, 0.8, 0.75, 0.8, 0.7), c(0.2, 0.25, 0.3, 0.2, 0.25),
                         c(0.8, 0.2, 0.75, 0.25, 0.5))[seq_len(k), , drop = FALSE]
    data <- as.data.frame(lapply(stats::setNames(seq_len(5L), paste0("u", seq_len(5L))),
                                 function(item) {
      factor(stats::rbinom(cell$n, 1L, probability[truth, item]), levels = 0:1)
    }))
    quiet(lca(data, paste0("u", seq_len(5L)), n_classes = k, n_starts = 5, seed = 1))
  } else {
    centre <- rbind(c(0, 0, 0, 0), c(2.5, 2, 2.5, 2), c(0, 2.5, -2.5, 1))[seq_len(k), , drop = FALSE]
    data <- as.data.frame(lapply(stats::setNames(seq_len(4L), paste0("y", seq_len(4L))),
                                 function(item) stats::rnorm(cell$n, centre[truth, item])))
    quiet(lpa(data, paste0("y", seq_len(4L)), n_profiles = k, model = "VVI", n_starts = 5,
              seed = 1))
  }
  residuals <- get_results(fit, "residuals", adjust = "BH")
  data.frame(kind = cell$kind, classes = k, n = cell$n, replication = replication,
             tests = nrow(residuals), rejected = sum(residuals$p_value < 0.05, na.rm = TRUE),
             flagged_unadjusted = any(residuals$p_value < 0.05, na.rm = TRUE),
             flagged = any(residuals$p_adjusted < 0.05, na.rm = TRUE),
             mean_statistic = mean(residuals$statistic, na.rm = TRUE))
}

started <- proc.time()[["elapsed"]]
cluster <- parallel::makePSOCKcluster(2L)
parallel::clusterSetRNGStream(cluster, 20261103)
invisible(parallel::clusterEvalQ(cluster, {
  Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  pkgload::load_all(".", quiet = TRUE)
}))
parallel::clusterExport(cluster, "replicate_cell")
results <- do.call(rbind, lapply(seq_len(nrow(cells)), function(index) {
  do.call(rbind, parallel::parLapply(cluster, seq_len(datasets), replicate_cell,
                                     cell = cells[index, ]))
}))
parallel::stopCluster(cluster)
summary_table <- do.call(rbind, lapply(split(results, with(results, list(kind, classes, n)),
                                             drop = TRUE), function(part) {
  with(part, data.frame(kind = kind[1L], classes = classes[1L], n = n[1L],
                        tests_per_dataset = mean(tests), per_test_rate = sum(rejected) / sum(tests),
                        flag_rate_unadjusted = mean(flagged_unadjusted),
                        flag_rate = mean(flagged), mean_statistic = mean(mean_statistic),
                        calibrated = abs(sum(rejected) / sum(tests) - 0.05) <= 0.025))
}))
summary_table <- summary_table[with(summary_table, order(kind, classes, n)), ]
print(summary_table, digits = 3, row.names = FALSE)
cat(sprintf("%.1f min\n", (proc.time()[["elapsed"]] - started) / 60))
if (datasets >= 100L) {
  utils::write.csv(summary_table, "validation/bivariate-residual-size.csv", row.names = FALSE)
}
