# Predeclared simulation study of multilpa(family = "additive").
# Design: validation/ADDITIVE_MODEL_DESIGN.md, "Simulation design and release
# decisions". The cell registry is written to additive-simulation-registry.csv
# BEFORE any replicate is run, and never edited afterwards.
#
# Run from the project root:
#   Rscript validation/additive-simulation.R pilot   # 100 datasets per cell
#   Rscript validation/additive-simulation.R full    # 1000 datasets per cell
# Outputs (validation/): additive-simulation-<stage>-replicates.csv,
# additive-simulation-<stage>-summary.csv, additive-simulation-<stage>.txt.

stage <- commandArgs(trailingOnly = TRUE)[1L] %||% "pilot"
stopifnot("stage must be `pilot` or `full`" = stage %in% c("pilot", "full"))
n_datasets <- if (identical(stage, "pilot")) 100L else 1000L
source(file.path("validation", "additive-simulation-setup.R"))

# ---- Run ----------------------------------------------------------------------------
started <- Sys.time()
jobs <- expand.grid(cell = names(cells), replicate = seq_len(n_datasets),
                    stringsAsFactors = FALSE)
replicates <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  # A script error is a recorded failure of that replicate, not a lost one.
  tryCatch(run_one(jobs$cell[k], jobs$replicate[k]), error = function(error) {
    failure_row(data.frame(cell = jobs$cell[k], replicate = jobs$replicate[k],
                           seed = NA_real_, stringsAsFactors = FALSE),
                "script_error", conditionMessage(error))
  })
}, mc.cores = max(1L, parallel::detectCores() - 1L), mc.set.seed = TRUE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
utils::write.csv(replicates, file.path(out_dir, sprintf(
  "additive-simulation-%s-replicates.csv", stage)), row.names = FALSE)

# ---- Summaries ------------------------------------------------------------------------
fitted_rows <- subset(replicates, fitted)
keys <- unique(fitted_rows[, c("cell", "parameter_id")])
summary_table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  rows <- subset(fitted_rows, cell == keys$cell[k] & parameter_id == keys$parameter_id[k])
  attempted <- n_datasets
  available <- is.finite(rows$se_observed)
  available_robust <- is.finite(rows$se_robust)
  cover <- function(low, high, ok) {
    if (!any(ok)) return(c(NA_real_, NA_real_))
    hit <- rows$truth[ok] >= low[ok] & rows$truth[ok] <= high[ok]
    p <- mean(hit)
    c(p, sqrt(p * (1 - p) / sum(ok)))
  }
  observed_cover <- cover(rows$low_observed, rows$high_observed, available)
  robust_cover <- cover(rows$low_robust, rows$high_robust, available_robust)
  empirical_sd <- stats::sd(rows$estimate)
  data.frame(
    cell = keys$cell[k], parameter_id = keys$parameter_id[k],
    truth = rows$truth[1L], attempted = attempted, fitted = nrow(rows),
    fit_and_interval = sum(available), availability = sum(available) / attempted,
    bias = mean(rows$estimate) - rows$truth[1L],
    bias_mcse = empirical_sd / sqrt(nrow(rows)),
    empirical_sd = empirical_sd,
    bias_over_sd = (mean(rows$estimate) - rows$truth[1L]) / empirical_sd,
    mean_se_observed = mean(rows$se_observed[available]),
    se_ratio_observed = mean(rows$se_observed[available]) /
      stats::sd(rows$estimate[available]),
    coverage_observed = observed_cover[1L], coverage_observed_mcse = observed_cover[2L],
    mean_se_robust = mean(rows$se_robust[available_robust]),
    se_ratio_robust = mean(rows$se_robust[available_robust]) /
      stats::sd(rows$estimate[available_robust]),
    coverage_robust = robust_cover[1L], coverage_robust_mcse = robust_cover[2L],
    width_observed = mean(rows$high_observed[available] - rows$low_observed[available]),
    stringsAsFactors = FALSE)
}))
summary_table$meets_release <- with(summary_table,
  availability >= 0.98 & abs(bias_over_sd) <= 0.1 &
    se_ratio_observed >= 0.9 & se_ratio_observed <= 1.1 &
    coverage_observed >= 0.925 & coverage_observed <= 0.975)
utils::write.csv(summary_table, file.path(out_dir, sprintf(
  "additive-simulation-%s-summary.csv", stage)), row.names = FALSE)

cell_status <- do.call(rbind, lapply(names(cells), function(name) {
  rows <- unique(subset(replicates, cell == name)[, c("replicate", "fitted",
                                                      "converged", "boundary",
                                                      "inference")])
  data.frame(cell = name, attempted = n_datasets, fitted = sum(rows$fitted),
             converged = sum(rows$converged %in% TRUE),
             boundary = sum(rows$boundary %in% TRUE),
             interval_available = sum(rows$inference == "available"),
             interval_withheld = sum(rows$inference == "withheld"),
             fit_failed = sum(rows$inference == "fit_failed"),
             script_error = sum(rows$inference == "script_error"))
}))
sink(file.path(out_dir, sprintf("additive-simulation-%s.txt", stage)))
cat(sprintf("Additive family simulation, stage %s: %d datasets per cell, %.1f minutes\n",
            stage, n_datasets, elapsed))
cat(sprintf("%s; latents %s\n\n", R.version.string, utils::packageVersion("latents")))
print(cell_status, row.names = FALSE)
cat("\n")
print(summary_table, digits = 3, row.names = FALSE)
sink()
cat(sprintf("done: %d rows, %.1f minutes\n", nrow(replicates), elapsed))
