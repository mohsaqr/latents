args <- commandArgs(trailingOnly = TRUE)
if (identical(args[[1L]], "old")) {
  library(latents)
  environment <- new.env(parent = asNamespace("latents"))
} else {
  pkgload::load_all(".", quiet = TRUE)
  environment <- new.env(parent = asNamespace("latents"))
}
results <- testthat::test_file("tests/testthat/test-astra-diagnostics.R", env = environment,
                              reporter = "summary", stop_on_failure = FALSE)
print(as.data.frame(results)[c("test", "failed", "error", "passed")])
saveRDS(results, file.path("validation/astra-fixes-2026-10-27/diagnostics",
                           paste0(args[[1L]], "-results.rds")))
