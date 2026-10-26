# Run the equivalence suites and write the report.
#
# From the project root:
#   Rscript equivalence/harness/run.R              # every suite
#   Rscript equivalence/harness/run.R glca mclust  # named suites only
#
# Writes report.csv (one row per compared quantity), suites.csv (one row per
# suite, including anything skipped or failed), SESSION.txt (the full
# sessionInfo the run was produced under) and REPORT.md (the summary). All four
# are generated, so no value in the documentation is transcribed by hand and no
# report can be attributed to a version it was not produced under.
#
# Also writes coverage.csv (every output column of every compared fit, with how
# many values were compared or why none were) and known-open-status.csv.
#
# The script exits non-zero when a suite failed, a quantity disagreed without a
# registered open defect to explain it, an output was neither compared nor
# exempted, or a known-open entry no longer matches anything, so that a stale
# harness cannot pass quietly in CI.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("equivalence", "harness", "registry.R"))

directory <- file.path("equivalence", "harness")
requested <- commandArgs(trailingOnly = TRUE)

# Check every validation script against the loaded package's signatures before
# anything is fitted. This is the check that was missing when multilpa 0.8.0
# renamed `indicators` to `vars` and `group` to `id`: the harness then died at
# suite load with a bare "unused arguments" naming neither file nor line.
stale <- check_validation_api(root = ".")
if (nrow(stale) > 0L) {
  cat(sprintf("%d validation call site(s) do not match multilpa %s:\n\n",
              nrow(stale), utils::packageVersion("latents")))
  print(stale)
  cat("\nFix these before the harness can be believed.\n")
  quit(status = 1L, save = "no")
}
cat(sprintf("API check: every validation call matches multilpa %s.\n",
            utils::packageVersion("latents")))

run <- run_equivalence(if (length(requested) > 0L) requested else NULL,
                       check_api = FALSE)

print(run)
cat("\n")
print(as.data.frame(run, what = "summary"))

comparisons <- as.data.frame(run, what = "comparisons")
comparisons$multilpa_version <- equivalence_provenance(run, "multilpa")
utils::write.csv(comparisons, file.path(directory, "report.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(run, what = "suites"),
                 file.path(directory, "suites.csv"), row.names = FALSE)
coverage <- as.data.frame(run, what = "coverage")
utils::write.csv(coverage, file.path(directory, "coverage.csv"), row.names = FALSE)
register <- as.data.frame(run, what = "known_open")
utils::write.csv(register, file.path(directory, "known-open-status.csv"), row.names = FALSE)

new_failures <- subset(comparisons, !agrees & is.na(known_open))
cat(sprintf("\n%d quantities (%d values) compared, %d agreed, %d disagreed (%d explained by a known open defect)\n",
            nrow(comparisons), sum(comparisons$n_values), sum(comparisons$agrees),
            sum(!comparisons$agrees), sum(!comparisons$agrees & !is.na(comparisons$known_open))))
if (nrow(new_failures) > 0L) {
  cat("\nDisagreements not explained by a known open defect:\n")
  print(new_failures[c("suite", "quantity", "reference", "obtained", "difference", "tolerance")])
}
uncovered <- subset(coverage, status == "not covered")
if (nrow(uncovered) > 0L) {
  cat("\nOutputs neither compared nor exempted:\n")
  print(uncovered[c("suite", "fit", "output", "n_produced")])
}
stale <- subset(register, status == "stale")
if (nrow(stale) > 0L) {
  cat("\nKnown-open entries that no longer match a disagreement (fixed? remove them):\n")
  print(stale[c("id", "suite", "pattern")])
}

writeLines(.equivalence_markdown(run), file.path(directory, "REPORT.md"))
writeLines(c(sprintf("multilpa %s", equivalence_provenance(run, "multilpa")),
             sprintf("generated %s", equivalence_provenance(run, "generated")),
             sprintf("platform %s", equivalence_provenance(run, "platform")),
             "", utils::capture.output(utils::sessionInfo())),
           file.path(directory, "SESSION.txt"))

cat(sprintf("\nwrote report.csv, suites.csv, coverage.csv, known-open-status.csv, REPORT.md and SESSION.txt in %s\n",
            directory))

suites <- as.data.frame(run, what = "suites")
# A known open defect does not fail the run; a new disagreement, a failed
# suite, an output nobody accounted for, or a register entry that no longer
# matches anything does.
if (any(suites$status == "failed") || nrow(new_failures) > 0L ||
    nrow(uncovered) > 0L || nrow(stale) > 0L) {
  quit(status = 1L, save = "no")
}
