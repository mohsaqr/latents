# Rscript equivalence/search/compare-norule.R <old_pilot_dir> <new_pilot_dir>
# The emEM default without the admissibility rule against the 2026-10-21
# pilot: share of fits that reach the best admissible maximum known (over
# every reported solution of both runs), share that report a spike (smallest
# relative variance below 0.01, equivalence/search/spike-measure.R), time.
args <- commandArgs(trailingOnly = TRUE)
read_rows <- function(dir, tag) {
  files <- list.files(dir, "_[0-9]{3}[.]rds$", full.names = TRUE)
  rows <- do.call(rbind, lapply(files, function(f) readRDS(f)$rows))
  rows$run <- tag
  rows
}
old <- read_rows(args[1], "with rule")
new <- read_rows(args[2], "without rule")
rows <- rbind(old, new)
key <- paste(rows$condition, rows$replicate, rows$k)
spike <- !is.na(rows$min_relative_variance) & rows$min_relative_variance < 0.01
ok <- !spike & is.finite(rows$log_likelihood)
best <- tapply(rows$log_likelihood[ok], key[ok], max)
rows$best <- unname(best[key])
tolerance <- 1e-6 * (1 + abs(rows$best))
rows$status <- ifelse(spike, "spike", ifelse(rows$best - rows$log_likelihood <= tolerance, "best", "miss"))
summarise <- function(g) data.frame(run = g$run[1], fit = g$fit[1], fits = nrow(g),
  best = round(mean(g$status == "best"), 3), spike = round(mean(g$status == "spike"), 3),
  median_seconds = round(median(g$seconds), 2))
pick <- subset(rows, fit %in% c("kmeans_10", "emem_10"))
out <- do.call(rbind, lapply(split(pick, list(pick$run, pick$fit), drop = TRUE), summarise))
print(out, row.names = FALSE)
cat("\nSpikes by condition and k (without rule, emem_10):\n")
nr <- subset(rows, run == "without rule" & fit == "emem_10")
print(aggregate(cbind(spike = status == "spike", best = status == "best") ~ condition + k,
                data = nr, FUN = mean), row.names = FALSE)
