pkgload::load_all(".", quiet = TRUE)
base <- "validation/release-fixes-2026-10-28/release"
run <- function(path, filter, name) {
  results <- testthat::test_dir(path, filter = filter, load_package = "none",
                                package = "latents", reporter = "summary",
                                stop_on_failure = FALSE)
  table <- as.data.frame(results)
  saveRDS(results, file.path(base, paste0(name, "-results.rds")))
  write.csv(table[setdiff(names(table), "result")],
            file.path(base, paste0(name, "-results.csv")), row.names = FALSE)
  cat("\n", name, "\n")
  print(colSums(table[c("nb", "failed", "skipped", "error", "warning", "passed")]))
  !any(table$failed > 0 | table$error)
}
reference_ok <- run("tests/equivalence", NULL, "external-tests")
retained_ok <- run("equivalence",
                   "^(latentgold.*|claims|mplus-(inference|robust|mixture-weights|lmr)|missing-patterns)$",
                   "retained-references")
quit(status = as.integer(!reference_ok || !retained_ok))
