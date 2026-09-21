# Run the equivalence suites and write the report.
#
# From the project root:
#   Rscript validation/equivalence/run.R              # every suite
#   Rscript validation/equivalence/run.R glca mclust  # named suites only
#
# Writes report.csv (one row per compared quantity), suites.csv (one row per
# suite, including anything skipped or failed), SESSION.txt (the full
# sessionInfo the run was produced under) and REPORT.md (the summary). All four
# are generated, so no value in the documentation is transcribed by hand and no
# report can be attributed to a version it was not produced under.
#
# The script exits non-zero when any quantity disagreed or any suite failed, so
# that a stale harness cannot pass quietly in CI.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("validation", "equivalence", "registry.R"))

directory <- file.path("validation", "equivalence")
requested <- commandArgs(trailingOnly = TRUE)

# Check every validation script against the loaded package's signatures before
# anything is fitted. This is the check that was missing when multilpa 0.8.0
# renamed `indicators` to `vars` and `group` to `id`: the harness then died at
# suite load with a bare "unused arguments" naming neither file nor line.
stale <- check_validation_api(root = ".")
if (nrow(stale) > 0L) {
  cat(sprintf("%d validation call site(s) do not match multilpa %s:\n\n",
              nrow(stale), utils::packageVersion("multilpa")))
  print(stale)
  cat("\nFix these before the harness can be believed.\n")
  quit(status = 1L, save = "no")
}
cat(sprintf("API check: every validation call matches multilpa %s.\n",
            utils::packageVersion("multilpa")))

run <- run_equivalence(if (length(requested) > 0L) requested else NULL,
                       check_api = FALSE)

print(run)
cat("\n")
print(get_data(run, "summary"))

comparisons <- get_data(run, "comparisons")
comparisons$multilpa_version <- equivalence_provenance(run, "multilpa")
utils::write.csv(comparisons, file.path(directory, "report.csv"), row.names = FALSE)
utils::write.csv(get_data(run, "suites"),
                 file.path(directory, "suites.csv"), row.names = FALSE)

cat(sprintf("\n%d quantities compared, %d agreed, %d disagreed, worst difference %.3e\n",
            nrow(comparisons), sum(comparisons$agrees), sum(!comparisons$agrees),
            max(c(comparisons$difference, 0))))
if (!all(comparisons$agrees)) {
  cat("\nDisagreements:\n")
  print(subset(comparisons, !agrees,
               select = c(suite, quantity, reference, obtained, difference, tolerance)))
}

writeLines(.equivalence_markdown(run), file.path(directory, "REPORT.md"))
writeLines(c(sprintf("multilpa %s", equivalence_provenance(run, "multilpa")),
             sprintf("generated %s", equivalence_provenance(run, "generated")),
             sprintf("platform %s", equivalence_provenance(run, "platform")),
             "", utils::capture.output(utils::sessionInfo())),
           file.path(directory, "SESSION.txt"))

cat(sprintf("\nwrote %s, %s, %s and %s\n",
            file.path(directory, "report.csv"), file.path(directory, "suites.csv"),
            file.path(directory, "REPORT.md"), file.path(directory, "SESSION.txt")))

suites <- get_data(run, "suites")
if (any(suites$status == "failed") || !all(comparisons$agrees)) quit(status = 1L, save = "no")
