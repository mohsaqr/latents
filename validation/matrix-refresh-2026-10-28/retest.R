# Rerun named equivalence/ test files after their failures were fixed and
# replace their rows in equivalence-results.csv (the earlier rows go to
# equivalence-results-first-run.csv once, for the record). Run from the
# project root with single-threaded BLAS:
#   VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript \
#     validation/matrix-refresh-2026-10-28/retest.R test-montecarlo.R ...
files <- commandArgs(trailingOnly = TRUE)
stopifnot("name at least one test file" = length(files) > 0L,
          "every file is in equivalence/" = all(file.exists(file.path("equivalence", files))))
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(".", quiet = TRUE)
base <- file.path("validation", "matrix-refresh-2026-10-28")

# Socket workers load latents from a library (see run-equivalence.R).
library_dir <- tempfile("latents-lib")
dir.create(library_dir)
status <- system2(file.path(R.home("bin"), "R"),
                  c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load",
                    paste0("--library=", shQuote(library_dir)), "."),
                  stdout = FALSE, stderr = FALSE)
stopifnot("installing the current latents for the socket workers failed" = status == 0L)
Sys.setenv(R_LIBS = paste(c(library_dir, Sys.getenv("R_LIBS")), collapse = .Platform$path.sep))

stems <- sub("^test-(.*)\\.R$", "\\1", files)
results <- testthat::test_dir(
  "equivalence", filter = paste0("^(", paste(gsub("([.+])", "\\\\\\1", stems), collapse = "|"), ")$"),
  load_package = "none", package = "latents", reporter = "summary", stop_on_failure = FALSE)

#' Replace the retested files' rows of the saved results; return the totals
#' of the retest
#' @noRd
merge_results <- function(results, files, base) {
  fresh <- as.data.frame(results)
  fresh <- fresh[setdiff(names(fresh), "result")]
  path <- file.path(base, "equivalence-results.csv")
  first <- file.path(base, "equivalence-results-first-run.csv")
  saved <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!file.exists(first)) file.copy(path, first)
  kept <- saved[!saved$file %in% files, ]
  merged <- rbind(kept, fresh[names(kept)])
  utils::write.csv(merged, path, row.names = FALSE)
  saveRDS(results, file.path(base, sprintf("retest-%s.rds", format(Sys.Date()))))
  # The retested files' rows of the first run, for the record.
  totals <- function(table) colSums(table[c("nb", "failed", "skipped", "error", "warning", "passed")])
  rbind(first_run = totals(saved[saved$file %in% files, ]), retest = totals(fresh),
        merged = totals(merged))
}
print(merge_results(results, files, base))
