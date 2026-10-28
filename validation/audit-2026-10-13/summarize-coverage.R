out <- "validation/audit-2026-10-13"
stopifnot(requireNamespace("covr", quietly = TRUE))
coverage <- readRDS(file.path(out, "coverage.rds"))
entries <- as.data.frame(coverage)
inventory <- read.csv(file.path(out, "function-inventory.csv"), stringsAsFactors = FALSE)
derived <- grepl("^(instrumented_expressions|expressions_hit|expression_coverage_percent|executed_in_installed_unit_suite)([.][0-9]+)?$", names(inventory))
inventory <- inventory[, !derived, drop = FALSE]
rows <- lapply(seq_len(nrow(inventory)), function(i) {
  body <- entries[entries$filename == inventory$file[i] & entries$functions == inventory$name[i], , drop = FALSE]
  data.frame(instrumented_expressions = nrow(body), expressions_hit = sum(body$value > 0),
             expression_coverage_percent = if (nrow(body)) 100 * mean(body$value > 0) else NA_real_,
             executed_in_installed_unit_suite = any(body$value > 0))
})
inventory <- cbind(inventory, do.call(rbind, rows))
inventory$review_status <- "whole-source structural screening; installed unit-suite coverage recorded; semantic findings and scope in REVIEW.md"
recent <- read.csv(file.path(out, "recent-api", "function-inventory.csv"), stringsAsFactors = FALSE)
core <- read.csv(file.path(out, "model-core", "function-inventory.csv"), stringsAsFactors = FALSE)
plot_rows <- readLines(file.path(out, "plots-docs", "function-inventory.txt"))
plot_names <- vapply(strsplit(plot_rows[grepl("^R/.*:[0-9]+", plot_rows)], "\t"), function(x) x[2L], character(1))
inventory$specialist_report <- ifelse(inventory$name %in% recent$function_name, "recent-api/findings.txt",
                                      ifelse(inventory$name %in% core$function_name, "model-core/README.md",
                                             ifelse(inventory$name %in% plot_names, "plots-docs/review.txt", "REVIEW.md shared services and coverage")))
if (file.exists(file.path(out, "usage-check.csv"))) {
  usage <- read.csv(file.path(out, "usage-check.csv"), stringsAsFactors = FALSE)
  inventory$usage_diagnostics <- usage$messages[match(paste(inventory$file, inventory$name), paste(usage$file, usage$name))]
}
utils::write.csv(inventory, file.path(out, "function-inventory.csv"), row.names = FALSE)
unhit <- inventory[inventory$definition == "function" & !inventory$executed_in_installed_unit_suite, ]
utils::write.csv(unhit, file.path(out, "functions-without-unit-hits.csv"), row.names = FALSE)
by_file <- lapply(unique(entries$filename), function(p) {
  x <- entries[entries$filename == p, , drop = FALSE]
  data.frame(file = p, expressions = nrow(x), hit_expressions = sum(x$value > 0),
             expression_percent = 100 * mean(x$value > 0))
})
utils::write.csv(do.call(rbind, by_file), file.path(out, "file-coverage.csv"), row.names = FALSE)
summary <- c(sprintf("Installed unit-suite line coverage: %.4f%%", covr::percent_coverage(coverage)),
             sprintf("Function definitions with observed expression hits: %d / %d", sum(inventory$executed_in_installed_unit_suite & inventory$definition == "function"), sum(inventory$definition == "function")),
             sprintf("Function definitions without hits in this coverage run: %d", nrow(unhit)),
             "Coverage measures installed unit tests only, not equivalence tests, help examples or audit probes.",
             "Passing instrumented statements do not establish branch coverage or numerical correctness.",
             "The raw coverage output succeeded, saved both results and reported coverage; the first wrapper ended on an unsupported by='file' summary argument. That argument is corrected in run-coverage.R. No tests were rerun to repair the summary.",
             "", "Functions without hits in the installed unit-suite coverage run", capture.output(print(unhit[c("file", "name")], row.names = FALSE)))
writeLines(summary, file.path(out, "coverage-summary.txt"))
cat(paste(summary[1:6], collapse = "\n"), "\n")
