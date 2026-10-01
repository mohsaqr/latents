# Evaluates the weak-class diagnostic exactly as declared in
# additive-diagnostic-protocol.md. Refits every registry dataset with its
# stored seed, computes the diagnostic statistics, joins them to the stored
# coverage records, chooses thresholds on replicates 1-500 and evaluates on
# 501-1000. Run from the project root:
#   Rscript validation/additive-diagnostic.R
# Outputs: additive-diagnostic-datasets.csv, additive-diagnostic.txt.

source(file.path("validation", "additive-simulation-setup.R"))
n_datasets <- 1000L

diagnose_one <- function(cell_name, replicate) {
  cell <- cells[[cell_name]]
  seed <- registry$base_seed[registry$cell == cell_name] + replicate
  data <- draw(cell, seed)
  fit <- withCallingHandlers(
    multilpa(data, c("y1", "y2"), "group", n_group_classes = cell$H,
             family = "additive", between_variance = cell$fit, seed = seed),
    warning = function(w) invokeRestart("muffleWarning"))
  inference <- tryCatch(.additive_inference(fit, "observed"),
    latents_boundary_fit = function(e) NULL,
    latents_no_converge = function(e) NULL,
    latents_singular_information = function(e) NULL)
  classification <- get_results(fit, "classification")
  avepp <- min(subset(classification, assigned == group_class)$mean_posterior,
               na.rm = FALSE)
  effective <- if (is.null(inference)) NA_real_ else
    min(.additive_effective_groups(fit, inference)$effective_groups)
  data.frame(cell = cell_name, replicate = replicate,
             effective_groups = effective,
             entropy = get_results(fit, "fit")$entropy,
             avepp = avepp, stringsAsFactors = FALSE)
}

started <- Sys.time()
# Datasets whose intervals were withheld are excluded by the protocol, so only
# those with intervals in the stored run are refitted.
stored <- utils::read.csv(file.path(out_dir, "additive-simulation-full-replicates.csv"),
                          stringsAsFactors = FALSE)
jobs <- unique(subset(stored, fitted & inference == "available",
                      c("cell", "replicate")))
diagnostics <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  diagnose_one(jobs$cell[k], jobs$replicate[k])
}, mc.cores = max(1L, parallel::detectCores() - 1L), mc.set.seed = TRUE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

# Coverage of class means and weights per dataset, from the stored records.
class_rows <- subset(stored, fitted & inference == "available" &
                       grepl("^(mean|weight)\\.", parameter_id))
class_rows$covered <- class_rows$truth >= class_rows$low_observed &
  class_rows$truth <= class_rows$high_observed
per_dataset <- stats::aggregate(covered ~ cell + replicate, data = class_rows,
                                FUN = function(v) c(n = length(v), hit = sum(v)))
per_dataset <- data.frame(cell = per_dataset$cell, replicate = per_dataset$replicate,
                          n_class_parameters = per_dataset$covered[, "n"],
                          n_covered = per_dataset$covered[, "hit"])
datasets <- merge(diagnostics, per_dataset, by = c("cell", "replicate"))
stopifnot("every dataset with intervals has a diagnostic" =
            nrow(datasets) == nrow(per_dataset),
          "the diagnostic exists wherever intervals were available" =
            all(is.finite(datasets$effective_groups)))
datasets$half <- ifelse(datasets$replicate <= 500L, "training", "evaluation")
utils::write.csv(datasets, file.path(out_dir, "additive-diagnostic-datasets.csv"),
                 row.names = FALSE)

reference <- c("regular_varying", "regular_shared", "unequal_sizes",
               "singleton_admixture", "three_classes")
choose <- function(statistic, grid) {
  training <- subset(datasets, half == "training" & cell %in% reference)
  rates <- vapply(grid, function(t) mean(training[[statistic]] < t), numeric(1))
  eligible <- grid[rates <= 0.01]
  list(threshold = if (length(eligible) > 0L) max(eligible) else NA_real_,
       training_rates = data.frame(threshold = grid, reference_flag_rate = rates))
}
rules <- list(
  effective_groups = choose("effective_groups", c(5, 10, 15, 20, 25, 30, 40, 50, 75, 100)),
  entropy = choose("entropy", seq(0.5, 0.9, by = 0.05)),
  avepp = choose("avepp", seq(0.6, 0.95, by = 0.05)))

evaluate <- function(statistic) {
  threshold <- rules[[statistic]]$threshold
  held <- subset(datasets, half == "evaluation")
  held$flag <- held[[statistic]] < threshold
  by_cell <- do.call(rbind, lapply(split(held, held$cell), function(d) {
    data.frame(cell = d$cell[1L], datasets = nrow(d), flag_rate = mean(d$flag))
  }))
  coverage <- do.call(rbind, lapply(c(FALSE, TRUE), function(f) {
    d <- subset(held, flag == f)
    p <- sum(d$n_covered) / sum(d$n_class_parameters)
    data.frame(flagged = f, datasets = nrow(d), parameters = sum(d$n_class_parameters),
               coverage = p, mcse = sqrt(p * (1 - p) / sum(d$n_class_parameters)))
  }))
  reference_rate <- mean(subset(held, cell %in% reference)$flag)
  target_rate <- mean(subset(held, cell == "rare_weak")$flag)
  criteria <- c(reference_false_alarm_at_most_2pct = reference_rate <= 0.02,
                rare_weak_detection_at_least_80pct = target_rate >= 0.80,
                unflagged_coverage_closer_to_95 =
                  abs(coverage$coverage[1L] - 0.95) < abs(coverage$coverage[2L] - 0.95))
  list(threshold = threshold, reference_rate = reference_rate,
       target_rate = target_rate, by_cell = by_cell, coverage = coverage,
       criteria = criteria)
}
results <- lapply(stats::setNames(names(rules), names(rules)), evaluate)

sink(file.path(out_dir, "additive-diagnostic.txt"))
cat(sprintf("Weak-class diagnostic evaluation (protocol: additive-diagnostic-protocol.md)\n%s; %.1f minutes; %d datasets with intervals\n\n",
            R.version.string, elapsed, nrow(datasets)))
invisible(lapply(names(results), function(name) {
  r <- results[[name]]
  cat(sprintf("== %s: threshold %s (flag when below), chosen on replicates 1-500\n",
              name, format(r$threshold)))
  cat("Training flag rates in reference cells:\n")
  print(rules[[name]]$training_rates, row.names = FALSE, digits = 3)
  cat(sprintf("\nHeld out (501-1000): reference flag rate %.4f; rare_weak flag rate %.4f\n",
              r$reference_rate, r$target_rate))
  print(r$by_cell, row.names = FALSE, digits = 3)
  cat("\nCoverage of class means and weights by flag (held out):\n")
  print(r$coverage, row.names = FALSE, digits = 4)
  cat("\nDeclared criteria:\n")
  print(r$criteria)
  cat("\n")
}))
sink()
cat(sprintf("done in %.1f minutes\n", elapsed))
