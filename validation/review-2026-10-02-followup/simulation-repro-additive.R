# Re-run sampled replicates of the additive simulation; compare to stored. Run from the repository root.
source(file.path("validation", "additive-simulation-setup.R"))
reps <- c(1L, 2L, 3L, 500L)
jobs <- expand.grid(cell = names(cells), replicate = reps, stringsAsFactors = FALSE)
new <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), \(k)
  run_one(jobs$cell[k], jobs$replicate[k]), mc.cores = 4L))
old <- read.csv("validation/additive-simulation-full-replicates.csv", stringsAsFactors = FALSE)
m <- merge(old, new, by = c("cell", "replicate", "parameter_id"), suffixes = c(".old", ".new"))
d <- abs(m$estimate.old - m$estimate.new)
ds <- abs(m$se_observed.old - m$se_observed.new)
cat("rows", nrow(m), "of", nrow(new), "| max |d est|", signif(max(d, na.rm = TRUE), 3),
    "| max |d se|", signif(max(ds, na.rm = TRUE), 3),
    "| inference changed:", sum(m$inference.old != m$inference.new), "\n")
print(aggregate(list(max_d_est = d), list(cell = m$cell), max))
