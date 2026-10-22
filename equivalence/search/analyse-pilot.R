# Rscript equivalence/search/analyse-pilot.R <pilot_dir>
# Per condition, number of classes and search: how often the reported
# solution is the best above-cutoff maximum known for the dataset, how often it
# is a spike, how far below it falls, and what it costs.
args <- commandArgs(trailingOnly = TRUE)
files <- list.files(args[1], "_[0-9]{3}[.]rds$", full.names = TRUE)
runs <- lapply(files, readRDS)
# Legacy pilot files predate this metadata. Only their known condition names
# establish whether a missing metric is inapplicable or was not measured.
#' Normalize the Gaussian applicability metadata in a saved pilot table
#' @param frame A pilot rows/finalists data frame, or `NULL`.
#' @return The frame with a logical `has_gaussian` column, or `NULL`.
#' @noRd
normalize_metrics <- function(frame) {
  stopifnot(is.null(frame) || is.data.frame(frame))
  if (is.null(frame)) return(NULL)
  stopifnot("condition" %in% names(frame))
  if (!"has_gaussian" %in% names(frame)) frame$has_gaussian <- NA
  unknown <- is.na(frame$has_gaussian)
  frame$has_gaussian[unknown & frame$condition %in%
    c("lpa_vvi_k3", "lpa_vvv_k2", "lpa_eee_k3", "ml_lpa")] <- TRUE
  frame$has_gaussian[unknown & frame$condition == "lca_k3"] <- FALSE
  frame
}
rows <- do.call(rbind, lapply(runs, function(run) normalize_metrics(run$rows)))
finals <- do.call(rbind, lapply(runs, function(run) normalize_metrics(run$finals)))
if (is.null(finals)) finals <- rows[FALSE, , drop = FALSE]
bound <- 0.01
tolerance <- function(value) 1e-6 * (1 + abs(value))
#' Distinguish measured, inapplicable and unmeasured benchmark diagnostics
#' @param has_gaussian Logical applicability flag, possibly missing.
#' @param ratio Numeric variance/eigenvalue metric, possibly missing.
#' @param cutoff Finite positive descriptive benchmark cutoff.
#' @return One character status per candidate. This does not select a fit.
#' @noRd
metric_status <- function(has_gaussian, ratio, cutoff = bound) {
  stopifnot(is.logical(has_gaussian), is.numeric(ratio),
            length(has_gaussian) == length(ratio), is.numeric(cutoff),
            length(cutoff) == 1L, is.finite(cutoff), cutoff > 0)
  ifelse(!is.na(has_gaussian) & !has_gaussian, "not_applicable",
    ifelse(!is.na(has_gaussian) & has_gaussian & is.finite(ratio),
      ifelse(ratio < cutoff, "spike", "above_cutoff"), "unclassified"))
}
rows$metric_status <- metric_status(rows$has_gaussian, rows$min_relative_variance)
finals$metric_status <- metric_status(finals$has_gaussian, finals$min_relative_variance)

# The reference over reported solutions and measured final emEM starts.
# Unmeasured Gaussian candidates never establish an above-cutoff reference.
key <- function(d) paste(d$condition, d$replicate, d$k)
candidates <- rbind(
  data.frame(key = key(rows), log_likelihood = rows$log_likelihood,
             metric_status = rows$metric_status),
  data.frame(key = key(finals), log_likelihood = finals$log_likelihood,
             metric_status = finals$metric_status))
reference <- subset(candidates, metric_status %in% c("above_cutoff", "not_applicable") &
                      is.finite(log_likelihood))
best <- tapply(reference$log_likelihood, reference$key, max)
rows$best <- unname(best[key(rows)])
rows$gap <- rows$best - rows$log_likelihood
rows$status <- ifelse(!is.finite(rows$log_likelihood), "error",
                ifelse(rows$metric_status == "unclassified", "unclassified",
                 ifelse(rows$metric_status == "spike", "spike",
                 ifelse(rows$gap <= tolerance(rows$best), "best", "miss"))))
metric_counts <- as.data.frame(table(source = c(rep("reported", nrow(rows)),
                                                rep("finalist", nrow(finals))),
  metric_status = c(rows$metric_status, finals$metric_status)))
cat("Metric coverage (unclassified candidates cannot establish the reference):\n")
print(subset(metric_counts, Freq > 0), row.names = FALSE)

summarise <- function(g) {
  n <- nrow(g)
  misses <- subset(g, status == "miss")
  data.frame(condition = g$condition[1], k = g$k[1], fit = g$fit[1], datasets = n,
             best = mean(g$status == "best"), miss = mean(g$status == "miss"),
             spike = mean(g$status == "spike"),
             unclassified = mean(g$status == "unclassified"), error = mean(g$status == "error"),
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
             unclassified = mean(g$status == "unclassified"), error = mean(g$status == "error"),
             median_seconds = median(g$seconds), total_seconds = sum(g$seconds))
}))
print(sorted_by_fit(overall, keys = character()), row.names = FALSE, digits = 3)

# Benchmark-only sensitivity: compare the measured finalists above each cutoff.
# No fallback silently admits unmeasured or below-cutoff Gaussian candidates.
cat("\nemEM finals re-selected at bound 0.05 (share whose choice changes):\n")
changed <- do.call(rbind, lapply(split(finals, list(key(finals), finals$fit), drop = TRUE),
                                 function(g) {
  pick <- function(b) {
    status <- metric_status(g$has_gaussian, g$min_relative_variance, b)
    ok <- status %in% c("above_cutoff", "not_applicable") & is.finite(g$log_likelihood)
    if (!any(ok)) return(NA_real_)
    max(g$log_likelihood[ok])
  }
  data.frame(fit = g$fit[1], condition = g$condition[1], k = g$k[1],
             changed = abs(pick(0.05) - pick(0.01)) > tolerance(pick(0.01)))
}))
if (!is.null(changed)) {
  sensitivity <- do.call(rbind, lapply(split(changed,
    list(changed$fit, changed$condition, changed$k), drop = TRUE), function(g) {
      data.frame(fit = g$fit[1L], condition = g$condition[1L], k = g$k[1L],
        changed = if (all(is.na(g$changed))) NA_real_ else mean(g$changed, na.rm = TRUE),
        unclassified = sum(is.na(g$changed)), fits = nrow(g))
    }))
  print(sensitivity, row.names = FALSE)
}
saveRDS(list(rows = rows, table = table, overall = overall, metric_counts = metric_counts),
        file.path(args[1], "summary.rds"))
