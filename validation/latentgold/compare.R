# Compare Latent GOLD's output against this package's targets.
#
# Run from the project root, after make-kit.R and after the Latent GOLD output
# has been copied into validation/latentgold/returned/:
#   Rscript validation/latentgold/compare.R
#
# Writes validation/latentgold/comparison.csv (one row per compared quantity)
# and validation/latentgold/REPORT.md. With no output returned yet it still
# runs, and reports every quantity as `missing`, which is the honest state.
#
# Status of each row:
#   agree            within the tolerance printed beside it
#   disagree         outside it: a finding
#   missing          Latent GOLD's file for this quantity is not in returned/
#   unparsed         the file is there but could not be read; the note says why
#   awaiting parser  only in the listing's tables, whose layout is not known
#                    until a real listing exists (see README.md, "Pilot")

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("validation", "latentgold")
source(file.path(directory, "compare-functions.R"))

target_files <- sort(list.files(file.path(directory, "targets"), pattern = "[.]rds$",
                                full.names = TRUE))
if (length(target_files) == 0L) {
  stop(errorCondition("No targets found. Run validation/latentgold/make-kit.R first.",
                      class = "multilpa_missing_targets", call = NULL))
}
returned <- file.path(directory, "returned")
targets <- lapply(target_files, readRDS)
built_with <- lg_target_versions(targets)
current <- as.character(utils::packageVersion("multilpa"))
if (!identical(built_with, current)) {
  warning(warningCondition(sprintf(
    "Targets were built with multilpa %s but %s is loaded; rebuild the kit if the fitting code changed.",
    paste(built_with, collapse = ", "), current), class = "multilpa_stale_targets"))
}

# The retained output is evidence for the kit it was produced from, and for no
# other; `record-run.R` fingerprinted that kit.
lg_check_fingerprint(file.path(directory, "kit"), returned)

comparison <- lg_compare_all(targets, returned)
lg_write_report(comparison, directory, current)

print(lg_status_counts(comparison), row.names = FALSE)
cat("\nEvery compared quantity:\n")
print(lg_compared(comparison), row.names = FALSE, digits = 8)
cat(sprintf("\n%s.\n", lg_totals(comparison)))
if (any(comparison$status == "disagree")) quit(save = "no", status = 1L)
