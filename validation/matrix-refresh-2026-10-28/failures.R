# Lists every failed or erroring expectation of the rerun, with its message.
#   Rscript validation/matrix-refresh-2026-10-28/failures.R
base <- file.path("validation", "matrix-refresh-2026-10-28")

#' One row per failed expectation: file, test, class and the first lines of
#' its message
#' @noRd
failed_expectations <- function(results) {
  rows <- lapply(results, function(test) {
    bad <- Filter(function(e) inherits(e, c("expectation_failure", "expectation_error")),
                  test$results)
    if (!length(bad)) return(NULL)
    data.frame(file = test$file, test = test$test,
               class = vapply(bad, function(e) class(e)[[1L]], character(1)),
               message = vapply(bad, function(e) substr(conditionMessage(e), 1L, 600L),
                                character(1)))
  })
  do.call(rbind, rows)
}

failures <- failed_expectations(readRDS(file.path(base, "equivalence-results.rds")))
utils::write.csv(failures, file.path(base, "failures.csv"), row.names = FALSE)
invisible(lapply(seq_len(nrow(failures)), function(i) {
  cat(sprintf("\n== %s :: %s [%s]\n%s\n", failures$file[[i]], failures$test[[i]],
              failures$class[[i]], failures$message[[i]]))
}))
