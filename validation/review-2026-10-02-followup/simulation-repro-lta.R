# Re-run a sample of archived LTA simulation replicates with the current code
# and compare to the stored values. Run from the repository root.
src <- readLines("validation/lta-simulation.R")
src <- src[seq_len(grep("^started <- Sys.time", src) - 1L)]
env <- new.env()
eval(parse(text = src), envir = env)
reps <- as.integer(strsplit(Sys.getenv("REPS", "1,2,3,250,500,750,1000"), ",")[[1]])
jobs <- expand.grid(cell = names(env$cells), replicate = reps, stringsAsFactors = FALSE)
new <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), \(k)
  env$run_one(jobs$cell[k], jobs$replicate[k]), mc.cores = 4L))
old <- read.csv("validation/lta-simulation-replicates.csv", stringsAsFactors = FALSE)
m <- merge(old, new, by = c("cell", "replicate", "key"), suffixes = c(".old", ".new"))
cols <- c("estimate", "se_observed", "se_robust", "low_observed", "high_observed")
diffs <- vapply(cols, \(v) max(abs(m[[paste0(v, ".old")]] - m[[paste0(v, ".new")]])), numeric(1))
cat("rows compared:", nrow(m), "of new", nrow(new), "\n")
print(signif(diffs, 3))
by_cell <- aggregate(abs(m$estimate.old - m$estimate.new), list(cell = m$cell), max)
print(by_cell)
