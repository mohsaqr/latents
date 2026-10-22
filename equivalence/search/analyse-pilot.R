# Rscript equivalence/search/analyse-pilot.R <pilot_dir>
# Per condition, number of classes and search: how often the reported
# solution is the best admissible maximum known for the dataset, how often it
# is a spike, how far below it falls, and what it costs.
args <- commandArgs(trailingOnly = TRUE)
files <- list.files(args[1], "_[0-9]{3}[.]rds$", full.names = TRUE)
runs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(runs, `[[`, "rows"))
finals <- do.call(rbind, lapply(runs, `[[`, "finals"))
bound <- 0.01
tolerance <- function(value) 1e-6 * (1 + abs(value))
spike <- function(ratio) !is.na(ratio) & ratio < bound

# The best admissible maximum known for a dataset and k: over every reported
# solution and every final emEM start that is not a spike.
key <- function(d) paste(d$condition, d$replicate, d$k)
candidates <- rbind(
  data.frame(key = key(rows), log_likelihood = rows$log_likelihood,
             ratio = rows$min_relative_variance),
  data.frame(key = key(finals), log_likelihood = finals$log_likelihood,
             ratio = finals$min_relative_variance))
admissible <- subset(candidates, !spike(ratio) & is.finite(log_likelihood))
best <- tapply(admissible$log_likelihood, admissible$key, max)
rows$best <- unname(best[key(rows)])
rows$gap <- rows$best - rows$log_likelihood
rows$status <- ifelse(is.na(rows$log_likelihood), "error",
                ifelse(spike(rows$min_relative_variance), "spike",
                 ifelse(rows$gap <= tolerance(rows$best), "best", "miss")))

summarise <- function(g) {
  n <- nrow(g)
  misses <- subset(g, status == "miss")
  data.frame(condition = g$condition[1], k = g$k[1], fit = g$fit[1], datasets = n,
             best = mean(g$status == "best"), miss = mean(g$status == "miss"),
             spike = mean(g$status == "spike"),
             median_gap_when_missed = if (nrow(misses)) median(misses$gap) else NA_real_,
             max_gap = if (nrow(misses)) max(misses$gap) else 0,
             median_seconds = median(g$seconds))
}
table <- do.call(rbind, lapply(split(rows, list(rows$condition, rows$k, rows$fit), drop = TRUE),
                               summarise))
order_fit <- c("kmeans_10", "kmeans_20", "kmeans_100", "emem_10", "emem_10_pool100",
               "emem_20", "reference_emem")
# Rows in condition, k and search order.
sorted_by_fit <- function(frame, keys = c("condition", "k")) {
  present <- intersect(keys, names(frame))
  ordering <- do.call(order, c(unname(as.list(frame[present])),
                               list(match(frame$fit, order_fit))))
  frame[ordering, , drop = FALSE]
}
table <- sorted_by_fit(table)
print(table, row.names = FALSE, digits = 3)

cat("\nOverall by search (all conditions and k):\n")
overall <- do.call(rbind, lapply(split(rows, rows$fit), function(g) {
  data.frame(fit = g$fit[1], fits = nrow(g), best = mean(g$status == "best"),
             miss = mean(g$status == "miss"), spike = mean(g$status == "spike"),
             median_seconds = median(g$seconds), total_seconds = sum(g$seconds))
}))
print(sorted_by_fit(overall, keys = character()), row.names = FALSE, digits = 3)

# Sensitivity to the admissibility bound: re-select each emEM fit's best final
# start under a stricter bound of 0.05.
cat("\nemEM finals re-selected at bound 0.05 (share whose choice changes):\n")
changed <- do.call(rbind, lapply(split(finals, list(key(finals), finals$fit), drop = TRUE),
                                 function(g) {
  pick <- function(b) {
    ok <- is.na(g$min_relative_variance) | g$min_relative_variance >= b
    if (!any(ok)) ok <- rep(TRUE, nrow(g))
    max(g$log_likelihood[ok])
  }
  data.frame(fit = g$fit[1], condition = g$condition[1], k = g$k[1],
             changed = abs(pick(0.05) - pick(0.01)) > tolerance(pick(0.01)))
}))
print(aggregate(changed ~ fit + condition + k, data = changed, FUN = mean), row.names = FALSE)
saveRDS(list(rows = rows, table = table, overall = overall),
        file.path(args[1], "summary.rds"))
