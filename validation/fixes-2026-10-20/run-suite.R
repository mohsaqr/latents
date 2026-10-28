arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L, arguments %in% c("unit", "equivalence"))
kind <- arguments[[1L]]
destination <- "validation/fixes-2026-10-20"
pkgload::load_all(".", quiet = TRUE)
source_files <- list.files("R", "[.]R$", full.names = TRUE)
before <- tools::md5sum(source_files)
started <- Sys.time()
results <- testthat::test_dir(
  if (kind == "unit") "tests/testthat" else "tests/equivalence",
  reporter = "summary", load_package = "none", stop_on_failure = FALSE)
saveRDS(results, file.path(destination, paste0(kind, "-tests.rds")))
table <- as.data.frame(results)
write.csv(table[, setdiff(names(table), "result")],
  file.path(destination, paste0(kind, "-tests.csv")), row.names = FALSE)
totals <- colSums(table[, c("nb", "failed", "error", "warning", "skipped", "passed")])
cat("\nTOTALS\n")
print(totals)
cat("Cases:", nrow(table), "\nElapsed:", as.numeric(difftime(Sys.time(), started,
  units = "secs")), "seconds\n")
unchanged <- identical(before, tools::md5sum(source_files))
cat("Production source unchanged during run:", unchanged, "\n")
write.csv(data.frame(file = names(before), md5 = unname(before)),
  file.path(destination, paste0(kind, "-source-manifest.csv")), row.names = FALSE)
print(sessionInfo())
if (!unchanged || totals[["failed"]] > 0 || totals[["error"]] > 0) quit(status = 1L)
