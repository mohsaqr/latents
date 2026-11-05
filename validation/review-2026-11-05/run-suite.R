args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
sys.source <- args[[1L]]
out <- args[[2L]]
options(mc.cores = 2L, warn = 1)
Sys.setenv(NOT_CRAN = "true", VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")
cat("Review runtime:", R.version.string, "\n")
cat("Source:", sys.source, "\n")
pkgload::load_all(sys.source, quiet = TRUE)
cat("Loaded version:", as.character(utils::packageVersion("latents")), "\n")
cat("Worker library:", find.package("latents"), "\n")
start <- proc.time()
results <- testthat::test_dir(file.path(sys.source, "tests", "testthat"),
  load_package = "none", reporter = "summary", stop_on_failure = FALSE)
saveRDS(results, file.path(out, "full-suite-results.rds"))
tab <- as.data.frame(results)
utils::write.csv(tab[, setdiff(names(tab), "result"), drop = FALSE],
  file.path(out, "full-suite-results.csv"), row.names = FALSE)
cat("\nTOTALS\n")
print(colSums(tab[, intersect(c("failed", "error", "warning", "skipped", "passed"), names(tab)), drop = FALSE]))
print(proc.time() - start)
