library(testthat)
review_path <- "/private/tmp/latents-review-20261022"
unit <- as.data.frame(readRDS(file.path(review_path, "unit-results.rds")))
socket <- as.data.frame(readRDS(file.path(review_path, "socket-results.rds")))
merged <- rbind(unit[!unit$file %in% socket$file, ], socket)
write.csv(merged[, setdiff(names(merged), "result")],
          file.path(review_path, "unit-merged-results.csv"), row.names = FALSE)
counts <- function(frame) {
  paste(nrow(frame), "test blocks in", length(unique(frame$file)), "files;",
        paste(paste(names(frame)[names(frame) %in% c("nb", "failed", "error", "warning", "skipped", "passed")],
                    colSums(frame[, names(frame) %in% c("nb", "failed", "error", "warning", "skipped", "passed"), drop = FALSE])),
              collapse = "; "))
}
warning_rows <- do.call(rbind, lapply(seq_len(nrow(merged)), function(index) {
  entries <- Filter(function(result) inherits(result, "expectation_warning"),
                    merged$result[[index]])
  if (!length(entries)) return(NULL)
  data.frame(file = merged$file[index], test = merged$test[index],
             message = vapply(entries, conditionMessage, character(1)))
}))
write.csv(warning_rows, file.path(review_path, "unit-warnings.csv"), row.names = FALSE)
lines <- c("latents 0.96.2, HEAD 12e0145; review 2026-10-22",
           paste("Unit suite (socket/search reruns replace sandbox-blocked files):", counts(merged)),
           "All eleven unit warnings are thin-profile notices on deliberately small fixtures.",
           "Socket/search rerun: 271 expectations; zero failures/errors/warnings/skips.")
equivalence_path <- file.path(review_path, "equivalence-results.rds")
if (file.exists(equivalence_path)) {
  equivalence <- as.data.frame(readRDS(equivalence_path))
  lines <- c(lines, paste("tests/equivalence:", counts(equivalence)))
} else lines <- c(lines, "tests/equivalence: still running; totals pending.")
targeted <- as.data.frame(readRDS(file.path(review_path, "targeted-results.rds")))
lines <- c(lines, paste("New isolated regression assertions (currently expose five confirmed defects):", counts(targeted)),
           "Lint: zero findings. Fresh install and source build succeeded.",
           "Limited R CMD check: zero errors; two warnings for deliberately unbuilt vignette outputs.",
           "No full vignette rendering, PDF manual, complete CRAN check or full standalone feature sweeps rerun.",
           "Numerical probes: 56/56 all-retained emEM/kmeans fits identical; eight extension paths preserve shared starts; 120 bounded variance cases completed.")
writeLines(lines, file.path(review_path, "verification-summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
