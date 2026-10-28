Sys.setenv(NOT_CRAN = "true")
out <- "validation/audit-2026-10-13"
staging <- file.path("tmp", "audit-2026-10-13-coverage")
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
files <- c("DESCRIPTION", "NAMESPACE", ".Rbuildignore", list.files(c("R", "man", "data", "tests"), recursive = TRUE, full.names = TRUE))
invisible(lapply(files, function(p) {
  target <- file.path(staging, p)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(p, target, overwrite = TRUE))
}))
coverage <- covr::package_coverage(path = staging, type = "tests", quiet = FALSE,
                                  pre_clean = FALSE)
saveRDS(coverage, file.path(out, "coverage.rds"))
frame <- as.data.frame(coverage)
utils::write.csv(frame, file.path(out, "line-coverage.csv"), row.names = FALSE)
cat("Line coverage percent:", covr::percent_coverage(coverage), "\n")
print(names(frame))
print(head(frame))
print(covr::tally_coverage(coverage, by = "line"))
