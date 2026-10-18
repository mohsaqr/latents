# Serial follow-up; no latents refitting.
options(digits = 12, width = 180)
library(mclust)
source(file.path("equivalence", "starts-study", "simulate.R"))
out <- file.path("equivalence", "starts-study", "codex")
comparison <- read.csv(file.path(out, "comparison-diagnostics.csv"))
cat("MCLUST VERSION", as.character(packageVersion("mclust")), "\n")
results <- lapply(c("lpa_k2_low_n300_r1", "lpa_k3_high_n300_r1", "lpa_k4_high_n300_r1"), function(name) {
  row <- subset(starts_design(), dataset == name)
  x <- as.matrix(starts_simulate(row)$data[paste0("y", 1:4)])
  set.seed(row$data_seed)
  timing <- system.time(mc <- Mclust(x, G = row$true_k + 1L, modelNames = "VVI", verbose = FALSE))
  str(mc, max.level = 1L)
  stopifnot(is.numeric(mc$loglik), is.finite(mc$loglik), mc$G == row$true_k + 1L)
  ref <- subset(comparison, dataset == name & k == row$true_k + 1L)
  result <- data.frame(dataset = name, k = mc$G, seconds = timing[["elapsed"]], mclust_ll = mc$loglik,
    latents_ll = ref$reference, mplus_ll = ref$mplus_best,
    mclust_minus_latents = mc$loglik - ref$reference, mclust_minus_mplus = mc$loglik - ref$mplus_best,
    smallest_hard_class = min(tabulate(mc$classification, mc$G)),
    smallest_effective_class = min(colSums(mc$z)), min_variance = min(mc$parameters$variance$sigmasq))
  print(result); result
})
write.csv(do.call(rbind, results), file.path(out, "mclust-check.csv"), row.names = FALSE)

# Boundary-only sensitivity: not a complete degeneracy screen.
source(file.path("equivalence", "starts-study", "analyse.R"))
raw <- starts_read()
nonboundary <- raw
nonboundary$log_likelihood[is.na(raw$boundary) | raw$boundary] <- -Inf
blocks <- rbind(starts_blocks(nonboundary, 10L), starts_blocks(nonboundary, 20L))
selection <- starts_selection(blocks)
cat("NONBOUNDARY REFERENCE CHOICES\n")
print(unique(selection[c("dataset", "family", "true_k", "reference_k")]), row.names = FALSE)
cat("NONBOUNDARY SELECTION DIRECTIONS\n")
print(with(selection, table(family, n_starts, direction)))
cat("NONBOUNDARY LPA4\n")
print(aggregate(cbind(fewer = direction < 0, more = direction > 0, true = chosen_k == true_k) ~ n_starts,
  subset(selection, family == "lpa" & true_k == 4), mean))
write.csv(selection, file.path(out, "nonboundary-selection.csv"), row.names = FALSE)
cat("FOLLOW-UP COMPLETE\n")
