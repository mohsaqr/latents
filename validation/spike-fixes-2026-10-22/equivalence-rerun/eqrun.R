# Run equivalence/ test files against a given package tree, never stopping early.
# Usage (from the repo root): Rscript eqrun.R <package_tree> <out.rds> [filter]
args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(args[[1]], quiet = TRUE)
testthat::set_max_fails(Inf)
filter <- if (length(args) >= 3L) args[[3]] else NULL
res <- testthat::test_dir("equivalence", filter = filter, package = "latents",
                          load_package = "none", stop_on_failure = FALSE,
                          reporter = testthat::ProgressReporter$new(show_praise = FALSE))
df <- as.data.frame(res)
saveRDS(df, args[[2]])
cat(sprintf("\nFILES %d TESTS %d EXPECTATIONS %d FAILED %d ERRORS %d SKIPPED %d WARNINGS %d\n",
            length(unique(df$file)), nrow(df), sum(df$nb), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
