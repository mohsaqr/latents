args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
cat(R.version.string, "\nlatents", as.character(utils::packageVersion("latents")),
    "\nLibrary", find.package("latents"), "\n")
results <- testthat::test_dir("tests/testthat", filter = args[[1L]],
  package = "latents", load_package = "installed", reporter = "summary",
  stop_on_failure = FALSE)
saveRDS(results, paste0(args[[2L]], "-results.rds"))
tab <- as.data.frame(results)
write.csv(tab[, setdiff(names(tab), "result"), drop = FALSE],
          paste0(args[[2L]], "-results.csv"), row.names = FALSE)
print(colSums(tab[, c("passed", "failed", "error", "warning", "skipped")]))
invisible(lapply(results, function(case) {
  unexpected <- Filter(function(x) inherits(x, c("expectation_failure",
    "expectation_error", "expectation_warning")), case$results)
  if (length(unexpected)) {
    cat("\n", case$file, ": ", case$test, "\n", sep = "")
    invisible(lapply(unexpected, function(x) {
      cat(paste(class(x), collapse = ", "), "\n", conditionMessage(x), "\n")
      if (!is.null(x$trace)) print(x$trace)
    }))
  }
}))
