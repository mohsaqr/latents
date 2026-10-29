# Rerun, on the current code, every equivalence/ test file whose latest saved
# result predates it (0.96.1) or that never saved one. Run from the project
# root with single-threaded BLAS:
#   VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript validation/matrix-refresh-2026-10-28/run-equivalence.R
# Left out: test-vlmr-runtime.R (3 workers, above the 2-core limit) and the
# runtime/infrastructure files that compare nothing with a reference.
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(".", quiet = TRUE)
base <- file.path("validation", "matrix-refresh-2026-10-28")

# Socket workers (macOS with Accelerate, Windows) load latents from a library,
# not from load_all(): without the current code installed where they look,
# they run whatever latents is installed (0.95.1 on 2026-10-29, which lacks
# .latents_collect_conditions(), so test-enumeration-mplus.R's two-core grid
# errored). The workers inherit R_LIBS from this process.
library_dir <- tempfile("latents-lib")
dir.create(library_dir)
status <- system2(file.path(R.home("bin"), "R"),
                  c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load",
                    paste0("--library=", shQuote(library_dir)), "."),
                  stdout = FALSE, stderr = FALSE)
stopifnot("installing the current latents for the socket workers failed" = status == 0L)
Sys.setenv(R_LIBS = paste(c(library_dir, Sys.getenv("R_LIBS")), collapse = .Platform$path.sep))

current <- utils::read.csv(file.path(
  "validation", "release-fixes-2026-10-28", "release", "retained-references-results.csv"))
excluded <- c("test-vlmr-runtime.R", "test-ml-three-step-runtime.R",
              "test-supporting-fork-preflight.R", "test-verification-artifact-export.R")
files <- setdiff(list.files("equivalence", "^test-.*\\.R$"), c(unique(current$file), excluded))
stems <- sub("^test-(.*)\\.R$", "\\1", files)
cat(sprintf("%d files:\n%s\n", length(files), paste(files, collapse = "\n")))

results <- testthat::test_dir(
  "equivalence", filter = paste0("^(", paste(gsub("([.+])", "\\\\\\1", stems), collapse = "|"), ")$"),
  load_package = "none", package = "latents", reporter = "summary", stop_on_failure = FALSE)
#' Save the results, their table and the session; return the totals
#' @noRd
save_results <- function(results, base) {
  table <- as.data.frame(results)
  saveRDS(results, file.path(base, "equivalence-results.rds"))
  utils::write.csv(table[setdiff(names(table), "result")],
                   file.path(base, "equivalence-results.csv"), row.names = FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(base, "SESSION.txt"))
  colSums(table[c("nb", "failed", "skipped", "error", "warning", "passed")])
}
print(save_results(results, base))
