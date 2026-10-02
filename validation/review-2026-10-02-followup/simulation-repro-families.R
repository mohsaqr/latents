# Re-run sampled replicates of the families simulation; compare to stored. Run from the repository root.
reps <- c(1L, 2L, 3L, 500L)
compare <- function(env, stored_file, key = "parameter_id") {
  jobs <- expand.grid(cell = names(env$cells), replicate = reps, stringsAsFactors = FALSE)
  new <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), \(k)
    env$run_one(jobs$cell[k], jobs$replicate[k]), mc.cores = 4L))
  old <- read.csv(stored_file, stringsAsFactors = FALSE)
  m <- merge(old, new, by = c("cell", "replicate", key), suffixes = c(".old", ".new"))
  d <- abs(m$estimate.old - m$estimate.new)
  ds <- abs(m$se_observed.old - m$se_observed.new)
  cat(stored_file, ": rows", nrow(m), "of", nrow(new),
      "| max |d est|", signif(max(d, na.rm = TRUE), 3),
      "| max |d se|", signif(max(ds, na.rm = TRUE), 3),
      "| NA est new:", sum(is.na(m$estimate.new)), "old:", sum(is.na(m$estimate.old)), "\n")
  print(aggregate(list(max_d_est = d), list(cell = m$cell), max))
}
src <- readLines("validation/families-simulation.R")
src <- src[seq_len(grep("^started <- Sys.time", src) - 1L)]
fam <- new.env(); eval(parse(text = src), envir = fam)
compare(fam, "validation/families-simulation-replicates.csv")
