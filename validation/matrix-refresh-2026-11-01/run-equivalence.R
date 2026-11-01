# Run, on the current code, the equivalence/ test files added or changed
# since the 2026-10-29 matrix refresh: the likelihood-ratio test (1.1.8),
# categorical distal outcomes against Latent GOLD, covariate models in every
# covariance structure against MoEClust, and the claims registry. Run from the
# project root with single-threaded BLAS:
#   VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript validation/matrix-refresh-2026-11-01/run-equivalence.R
# The runtime source is captured before the run and bound to the result
# (equivalence-results.csv.provenance.csv), as the matrix requires.
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(".", quiet = TRUE)
source(file.path("equivalence", "verification-matrix", "provenance.R"))
base <- file.path("validation", "matrix-refresh-2026-11-01")
snapshot <- matrix_capture_source(".")

files <- c("test-lrt.R", "test-latentgold-categorical-distal.R",
           "test-moeclust-covariate-structures.R", "test-claims.R")
stems <- sub("^test-(.*)\\.R$", "\\1", files)
results <- testthat::test_dir(
  "equivalence", filter = paste0("^(", paste(gsub("([.+])", "\\\\\\1", stems), collapse = "|"), ")$"),
  load_package = "none", package = "latents", reporter = "summary", stop_on_failure = FALSE)
table <- as.data.frame(results)
path <- file.path(base, "equivalence-results.csv")
utils::write.csv(subset(table, select = -result), path, row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(base, "SESSION.txt"))
matrix_record_evidence(path, snapshot, ".")
print(colSums(subset(table, select = c(nb, failed, skipped, error, warning, passed))))
