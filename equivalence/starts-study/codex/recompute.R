# Independent, serial audit. Run from the repository root with BLAS threads = 1.
options(digits = 12, width = 180)
out <- file.path("equivalence", "starts-study", "codex")
folder <- file.path("equivalence", "starts-study")
source(file.path(folder, "simulate.R"))
source(file.path(folder, "analyse.R"))
source(file.path(folder, "compare-mplus.R"))
raw <- do.call(rbind, lapply(list.files(file.path(folder, "fits"), full.names = TRUE), readRDS))
str(raw)
print(head(raw)); print(summary(raw)); print(vapply(raw, class, character(1)))
stopifnot(nrow(raw) == nrow(unique(raw[c("dataset", "k", "start")])), !anyNA(raw$log_likelihood))
cells <- split(raw, paste(raw$dataset, raw$k))
stopifnot(all(vapply(cells, function(z) identical(sort(z$start), 1:500), logical(1))))
cat("DATASETS", length(unique(raw$dataset)), "CELLS", length(cells), "ROWS", nrow(raw), "\n")
print(table(unique(raw[c("dataset", "family", "true_k", "separation", "n", "rep")])$rep))
print(with(unique(raw[c("dataset", "family", "true_k", "separation", "n")]), table(family, true_k, separation, n)))
# Independent explicit slicing rather than starts_blocks' split-and-pmax pipeline.
blocks <- do.call(rbind, lapply(c(10L, 20L), function(s) do.call(rbind, lapply(cells, function(z) {
  z <- z[order(z$start), ]
  do.call(rbind, lapply(seq_len(499L %/% (s - 1L)), function(b) {
    indices <- c(1L, 1L + (b - 1L) * (s - 1L) + seq_len(s - 1L))
    ll <- max(z$log_likelihood[indices]); ref <- max(z$log_likelihood)
    data.frame(z[1, c("dataset", "family", "true_k", "separation", "n", "rep", "k", "n_parameters")],
      n_starts = s, block = b, log_likelihood = ll, reference = ref, gap = ll - ref,
      reached = ll - ref > -0.001, ward_reached = z$log_likelihood[1] - ref > -0.001,
      single_start_rate = mean(z$log_likelihood - ref > -0.001), row.names = NULL)
  }))
}))))
order_blocks <- function(x) x[order(x$dataset, x$k, x$n_starts, x$block), ]
saved_blocks <- read.csv(file.path(folder, "blocks.csv"))
a <- order_blocks(blocks); b <- order_blocks(saved_blocks); rownames(a) <- rownames(b) <- NULL
stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10)))
cat("BLOCKS verified", nrow(blocks), "\n")
selection <- do.call(rbind, lapply(split(blocks, paste(blocks$dataset, blocks$n_starts, blocks$block)), function(z) {
  bic <- -2 * z$log_likelihood + z$n_parameters * log(z$n)
  ref_bic <- -2 * z$reference + z$n_parameters * log(z$n)
  stopifnot(length(unique(z$k)) == 4L, length(unique(z$n)) == 1L)
  chosen <- z$k[order(bic, z$k)[1L]]; reference <- z$k[order(ref_bic, z$k)[1L]]
  data.frame(z[1, c("dataset", "family", "true_k", "separation", "n", "n_starts", "block")],
    chosen_k = chosen, reference_k = reference, agrees = chosen == reference,
    direction = sign(chosen - reference), reference_margin = diff(sort(ref_bic))[1L], row.names = NULL)
}))
order_selection <- function(x) x[order(x$dataset, x$n_starts, x$block), ]
a <- order_selection(selection); b <- order_selection(read.csv(file.path(folder, "selection.csv"))); rownames(a) <- rownames(b) <- NULL
stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10)))
cat("SELECTION verified", nrow(selection), "\n")
stopifnot(all(raw$n_parameters == ifelse(raw$family == "lpa", 9 * raw$k - 1, 7 * raw$k - 1)))
saved_summaries <- readRDS(file.path(folder, "summaries.rds"))
stopifnot(isTRUE(all.equal(starts_summaries(saved_blocks, read.csv(file.path(folder, "selection.csv"))), saved_summaries)))
cat("SUMMARIES verified\n")
blocks$relative_k <- blocks$k - blocks$true_k
cat("REACH COUNTS\n"); print(aggregate(cbind(hits = as.integer(reached), total = rep(1L, nrow(blocks))) ~ family + relative_k + n_starts, blocks, sum))
cat("ALL FAMILY REACH\n"); print(aggregate(cbind(hits = as.integer(reached), total = rep(1L, nrow(blocks))) ~ relative_k + n_starts, blocks, sum))
cat("TRUE K FAILING DATASETS\n"); print(aggregate(cbind(reached, gap) ~ dataset + n_starts, subset(blocks, k == true_k & !reached), function(x) c(n = length(x), min = min(x))))
cat("SELECT DIRECTIONS\n"); print(with(selection, table(family, n_starts, direction)))
cat("SELECT BY K\n"); print(aggregate(cbind(fewer = as.integer(direction < 0), more = as.integer(direction > 0), true_run = as.integer(chosen_k == true_k), true_ref = as.integer(reference_k == true_k), total = rep(1L, nrow(selection))) ~ family + true_k + n_starts, selection, sum))
changed <- subset(selection, direction != 0)
cat("CHANGED CHOICES\n"); print(aggregate(block ~ dataset + n_starts + chosen_k + reference_k + true_k, changed, length))
cell_summary <- do.call(rbind, lapply(cells, function(z) {
  ref <- max(z$log_likelihood); near <- z$log_likelihood > ref - 0.001
  data.frame(z[1, c("dataset", "family", "true_k", "separation", "n", "k", "n_parameters", "smallest_class", "seconds", "warnings")],
    reference = ref, selected_ll = z$final_log_likelihood[1], selected_difference = z$final_log_likelihood[1] - ref,
    best_hits = sum(near), best_converged = sum(near & z$converged), best_boundary = sum(near & z$boundary, na.rm = TRUE),
    failed = sum(!is.finite(z$log_likelihood)), nonconverged = sum(!z$converged), ward_reached = near[z$start == 1])
}))
cat("CELL DIAGNOSTICS\n"); print(cell_summary, row.names = FALSE)
write.csv(cell_summary, file.path(out, "cell-diagnostics.csv"), row.names = FALSE)
cat("TOLERANCE SENSITIVITY\n"); print(sapply(c(1e-4, 1e-3, 1e-2), function(t) with(blocks, tapply(gap > -t, list(relative_k, n_starts), mean))))
# Parse the independent MODEL FIT INFORMATION section as a check on the ranked list.
paths <- list.files(file.path(folder, "mplus"), pattern = "[.]out$", full.names = TRUE)
mplus <- do.call(rbind, lapply(paths, function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines))) return(NULL)
  z <- starts_mplus_read(path)
  h0 <- grep("^\\s+H0 Value\\s+", lines, value = TRUE)
  p <- grep("^Number of Free Parameters\\s+", lines, value = TRUE)
  stopifnot(length(h0) == 1L, length(p) == 1L)
  z$h0 <- as.numeric(tail(strsplit(trimws(h0), "[[:space:]]+")[[1L]], 1))
  z$n_parameters_mplus <- as.integer(tail(strsplit(trimws(p), "[[:space:]]+")[[1L]], 1))
  z$nonpd <- any(grepl("NON-POSITIVE DEFINITE", lines))
  z
}))
stopifnot(all(mplus$h0 == mplus$mplus_best))
comparison <- read.csv(file.path(folder, "mplus-comparison.csv"))
comparison <- merge(comparison, cell_summary[c("dataset", "k", "smallest_class", "best_boundary", "best_hits")], by = c("dataset", "k"))
comparison <- merge(comparison, mplus[c("dataset", "k", "h0", "n_parameters_mplus", "nonpd")], by = c("dataset", "k"))
stopifnot(all(comparison$h0 == comparison$mplus_best), all(comparison$n_parameters_mplus == ifelse(comparison$family == "lpa", 9 * comparison$k - 1, 7 * comparison$k - 1)))
cat("MPLUS", nrow(mplus), "completed outputs", nrow(comparison), "matched cells\n")
print(aggregate(latents_minus_mplus ~ family + relative_k, comparison, function(x) c(n = length(x), min = min(x), median = median(x), max = max(x))))
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(out, "comparison-diagnostics.csv"), row.names = FALSE)
cat("MPLUS REPLICATION COUNTS\n"); print(with(comparison, table(family, relative_k, mplus_replicated, mplus_not_replicated_warning)))

pkgload::load_all(".", quiet = TRUE)
# Check actual package tie selection on every saved 500-start cell.
package_gaps <- vapply(cells, function(z) z$log_likelihood[.multilpa_select_start(z$log_likelihood, z$converged)] - max(z$log_likelihood), numeric(1))
cat("PACKAGE SELECTION difference range", range(package_gaps), "\n")
cat("REFIT BEGIN\n"); flush.console()
row <- subset(starts_design(), dataset == "lpa_k4_high_n300_r1")
sim <- starts_simulate(row)
str(sim); print(head(sim$data)); print(summary(sim$data))
warnings <- character()
timing <- system.time(fit <- withCallingHandlers(lpa(sim$data, paste0("y", 1:4), n_profiles = 5L,
  model = "VVI", n_starts = 500L, seed = row$data_seed), warning = function(w) {
    warnings <<- c(warnings, paste(class(w)[1L], conditionMessage(w))); invokeRestart("muffleWarning")
  }))
cat("REFIT elapsed", timing[["elapsed"]], "seed", row$data_seed, "\n")
cat("REFIT names\n"); print(names(fit)); str(fit$variances)
counts <- tabulate(fit$subject_profiles, 5L); effective <- colSums(fit$subject_posteriors)
small <- which.min(counts)
cat("REFIT LL", fit$log_likelihood, "n_parameters", fit$n_parameters, "BIC", fit$bic, "\n")
cat("REFIT hard counts", counts, "effective", effective, "smallest index", small, "\n")
cat("REFIT smallest variances\n"); print(fit$variances[small, ])
cat("REFIT smallest means\n"); print(fit$means[small, ])
cat("REFIT relative variance\n"); print(fit$variances[small, ] / vapply(sim$data[paste0("y", 1:4)], var, numeric(1)))
cat("REFIT boundary", fit$boundary, "small_classes", fit$small_classes, "warnings", warnings, "\n")
cat("REFIT subject\n"); print(sim$data[fit$subject_profiles == small, ])
stopifnot(abs(fit$log_likelihood - subset(cell_summary, dataset == row$dataset & k == 5)$reference) < 0.001)
saveRDS(fit, file.path(out, "singleton-refit.rds"))
saveRDS(list(timing = timing, counts = counts, effective = effective, smallest = small, warnings = warnings, seed = row$data_seed), file.path(out, "refit-diagnostics.rds"))
# One small independent optimizer check, one default hierarchical start per dataset.
if (requireNamespace("mclust", quietly = TRUE)) {
  library(mclust)
  cat("MCLUST VERSION", as.character(packageVersion("mclust")), "\n")
  checks <- lapply(c("lpa_k2_low_n300_r1", "lpa_k3_high_n300_r1", "lpa_k4_high_n300_r1"), function(name) {
    row <- subset(starts_design(), dataset == name); x <- as.matrix(starts_simulate(row)$data[paste0("y", 1:4)])
    set.seed(row$data_seed)
    tm <- system.time(mc <- mclust::Mclust(x, G = row$true_k + 1L, modelNames = "VVI", verbose = FALSE))
    stopifnot(is.list(mc), is.numeric(mc$loglik))
    ref <- subset(comparison, dataset == name & k == row$true_k + 1L)
    result <- data.frame(dataset = name, k = row$true_k + 1L, seconds = tm[["elapsed"]],
      mclust_ll = mc$loglik, latents_ll = ref$reference, mplus_ll = ref$mplus_best,
      mclust_minus_latents = mc$loglik - ref$reference, mclust_minus_mplus = mc$loglik - ref$mplus_best,
      smallest_hard_class = min(tabulate(mc$classification, row$true_k + 1L)),
      min_variance = min(mc$parameters$variance$sigmasq))
    print(result); result
  })
  write.csv(do.call(rbind, checks), file.path(out, "mclust-check.csv"), row.names = FALSE)
}
cat("AUDIT COMPLETE\n")
print(sessionInfo())
