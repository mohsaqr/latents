# Lint the named files with the project's .lintr, on the ci VM:
#   mini-ci -- Rscript equivalence/argument-sweep/diagnostics/lint-files.R R/search.R R/predict.R
# Exits non-zero when any lint is found.
paths <- commandArgs(trailingOnly = TRUE)
stopifnot("name at least one file" = length(paths) > 0L, "every path must exist" = all(file.exists(paths)))
if (!requireNamespace("lintr", quietly = TRUE)) utils::install.packages("lintr")
lints <- do.call(c, lapply(paths, lintr::lint))
print(lints)
cat(sprintf("%d lints in %d files\n", length(lints), length(paths)))
if (length(lints)) quit(status = 1L)
