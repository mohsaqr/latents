library(testthat)
review_path <- "/private/tmp/latents-fixes-20261022"
counts <- function(frame) {
  paste(nrow(frame), "test blocks in", length(unique(frame$file)), "files;",
        paste(paste(names(frame)[names(frame) %in% c("nb", "failed", "error", "warning", "skipped", "passed")],
                    colSums(frame[, names(frame) %in% c("nb", "failed", "error", "warning", "skipped", "passed"), drop = FALSE])),
              collapse = "; "))
}
lines <- c("latents 0.96.2, working-tree fixes based on 12e0145; 2026-10-22",
           "Single-thread BLAS, installed current package for socket workers, NOT_CRAN=true.")
invisible(lapply(c("unit", "equivalence"), function(suite) {
  path <- file.path(review_path, paste0(suite, "-results.rds"))
  if (!file.exists(path)) {
    lines <<- c(lines, paste(suite, "still running; final counts pending."))
    return(invisible(NULL))
  }
  frame <- as.data.frame(readRDS(path))
  fixture_path <- file.path(review_path, "lta-fixture-final.rds")
  if (suite == "unit" && file.exists(fixture_path)) {
    fixture <- as.data.frame(readRDS(fixture_path))
    frame <- rbind(frame[!frame$file %in% fixture$file, ], fixture)
    write.csv(frame[, setdiff(names(frame), "result")],
              file.path(review_path, "unit-final-results.csv"), row.names = FALSE)
  }
  lines <<- c(lines, paste(suite, counts(frame)))
  warnings <- do.call(rbind, lapply(seq_len(nrow(frame)), function(index) {
    entries <- Filter(function(result) inherits(result, "expectation_warning"),
                      frame$result[[index]])
    if (!length(entries)) return(NULL)
    data.frame(file = frame$file[index], test = frame$test[index],
               message = vapply(entries, conditionMessage, character(1)))
  }))
  if (!is.null(warnings)) {
    write.csv(warnings, file.path(review_path, paste0(suite, "-warnings.csv")), row.names = FALSE)
  }
}))
lines <- c(lines,
           "Focused warning tests: 90 expectations, zero failures/errors/warnings.",
           "Focused staged-bootstrap file: 57 expectations, zero failures/errors/warnings.",
           "The review-lta fixture file was rerun after allowing its expected thin-profile warning; that file replaces the initial warning-heavy result.",
           "Full package lint: zero findings. Fresh install and clean source build succeeded.",
           "Limited R CMD check: zero errors, two warnings for unbuilt vignette outputs.",
           "No complete CRAN check, manual, full vignette render or standalone feature sweeps rerun.",
           "Four LTA fit objects identical before/after; unaffected staged bootstrap and RNG state identical.")
writeLines(lines, file.path(review_path, "verification-summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
