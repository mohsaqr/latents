# Runs every latent-transition extension comparison and writes REPORT.md.
# Run from the project root with JStats checked out beside this repository.
scripts <- c("compare-ex814.R", "compare-ex813.R", "compare-mplus-measurement.R",
             "compare-depmix-covariate.R", "compare-lmest-occasion.R")
base <- file.path("equivalence", "lta-extensions")
invisible(lapply(scripts, function(script) {
  system2("Rscript", file.path(base, script), stdout = FALSE, stderr = FALSE)
}))
results <- do.call(rbind, lapply(list.files(base, "^result-.*\\.rds$", full.names = TRUE),
                                 readRDS))
utils::write.csv(results, file.path(base, "comparison.csv"), row.names = FALSE)
lines <- c("# Latent transition extensions: external comparison", "",
  sprintf("Run %s; %s. **%d of %d quantities agree.**", format(Sys.Date()),
          R.version.string, sum(results$agree), nrow(results)), "",
  "| Case | Quantity | latents | Reference | Difference | Tolerance | Agree |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %.6g | %.6g | %.2e | %.0e | %s |", results$case,
          results$quantity, results$latents, results$reference, results$difference,
          results$tolerance, ifelse(results$agree, "yes", "**no**")))
writeLines(lines, file.path(base, "REPORT.md"))
cat(sum(results$agree), "of", nrow(results), "agree\n")
