pkgload::load_all(".", quiet = TRUE)
directory <- file.path("equivalence", "latentgold")
destination <- file.path("validation", "fixes-2026-10-20", "latentgold")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
source(file.path(directory, "cases.R"))
source(file.path(directory, "targets.R"))
source(file.path(directory, "compare-functions.R"))
lg_check_fingerprint(file.path(directory, "kit"), file.path(directory, "returned"))
cases <- lg_cases()
targets <- Map(function(name, case) {
  message("Recomputing retained reference case: ", name)
  target <- .targets(name, case, case$fit())
  saveRDS(target, file.path(destination, paste0(name, ".rds")))
  target
}, names(cases), cases)
comparison <- lg_compare_all(targets, file.path(directory, "returned"))
write.csv(comparison, file.path(destination, "comparison.csv"), row.names = FALSE)
print(lg_status_counts(comparison), row.names = FALSE)
cat(lg_totals(comparison), "\n")
cat("Compared with retained external output; no new external-program run.\n")
if (any(comparison$status == "disagree")) {
  print(comparison[comparison$status == "disagree", ])
  quit(status = 1L)
}
