# Run the equivalence suites and write the report.
#
# From the project root:
#   Rscript validation/equivalence/run.R              # every suite
#   Rscript validation/equivalence/run.R glca mclust  # named suites only
#
# Writes report.csv (one row per compared quantity) and REPORT.md (the summary
# table). Both are generated, so no value in the documentation is transcribed
# by hand.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("validation", "equivalence", "registry.R"))

requested <- commandArgs(trailingOnly = TRUE)
report <- run_equivalence(if (length(requested) > 0L) requested else NULL)
directory <- file.path("validation", "equivalence")
utils::write.csv(report, file.path(directory, "report.csv"), row.names = FALSE)

summary_table <- equivalence_summary(report)
print(summary_table)
cat(sprintf("\n%d quantities compared, %d agreed, worst difference %.3e\n",
            nrow(report), sum(report$agrees), max(report$difference)))
if (!all(report$agrees)) {
  cat("\nDisagreements:\n")
  print(subset(report, !agrees,
               select = c(suite, quantity, reference, obtained, difference)))
}
writeLines(.equivalence_markdown(report, summary_table),
           file.path(directory, "REPORT.md"))
cat(sprintf("\nwrote %s and %s\n", file.path(directory, "report.csv"),
            file.path(directory, "REPORT.md")))
