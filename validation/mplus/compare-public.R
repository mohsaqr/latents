# Run from the package root: Rscript validation/mplus/compare-public.R
# Official examples 7.9 and 7.10 validate ONLY the single-level LPA limit.
# This script never fits the source column c: it is a generating class label.
source(file.path("R", "fit-ml-lpa.R"))
source(file.path("R", "gaussian-moments.R"))
source(file.path("tests", "testthat", "helper-mplus-public.R"))
dir.create(file.path("tests", "fixtures", "mplus"), recursive = TRUE, showWarnings = FALSE)

# Exercise the parser on simulated output before reading official data/results.
testthat::test_file(file.path("tests", "testthat", "test-mplus-public-parser.R"))

#' Reproduce one published single-level Mplus example from the original artifacts
#' @param example Official example identifier, 7.9 or 7.10.
#' @return Discrepancy table; also saves an offline regression fixture.
run_public_comparison <- function(example) {
  stopifnot(is.character(example), length(example) == 1L, example %in% c("7.9", "7.10"))
  paths <- file.path("validation", "mplus", "public", paste0("ex", example, c(".dat", ".inp", ".html")))
  stopifnot(all(file.exists(paths)))
  data <- tryCatch(utils::read.table(paths[1L], col.names = c("y1", "y2", "y3", "y4", "c")),
                   error = function(error) stop(sprintf("Cannot read example %s: %s", example, conditionMessage(error))))
  str(data)
  print(head(data))
  print(summary(data))
  print(vapply(data, class, character(1)))
  print(dim(data))
  stopifnot(identical(dim(data), c(500L, 5L)), !anyNA(data), !anyDuplicated(data),
            all(vapply(data, is.numeric, logical(1))), all(data$c %in% c(1L, 2L)))
  input <- readLines(paths[2L])
  print(input)
  expected <- parse_public_mplus(readLines(paths[3L]))
  variance_model <- if (example == "7.9") "equal" else "varying"
  data$cluster <- 1L
  stopifnot(!anyNA(data), !anyDuplicated(data), all(data$cluster == 1L))
  fit <- fit_multilpa(data, c("y1", "y2", "y3", "y4"), "cluster", 2L, 1L,
                    variance_model = variance_model, n_starts = 20L, max_iter = 1000L,
                    tol = 1e-12, seed = 912)
  comparison <- compare_public_mplus(fit, expected)
  print(comparison, digits = 9)
  stopifnot(fit$converged, !fit$boundary, fit$n_failed_starts == 0L,
            fit$n_best_replicated >= 2L, all(comparison$passed))
  provenance <- list(
    urls = paste0("https://www.statmodel.com/usersguide/chap7/ex", example, c(".dat", ".inp", ".html")),
    md5 = tools::md5sum(paths), mplus_version = "8.8", output_date = "2022-04-19",
    retrieved = "2026-09-17", scope = "single-level LPA only",
    tolerance_note = "Absolute half printed unit, except effective counts: 0.001 subject optimization allowance.")
  saveRDS(list(data = data, expected = expected, variance_model = variance_model,
               provenance = provenance), file.path("tests", "fixtures", "mplus", paste0("public-", example, ".rds")))
  comparison$example <- example
  comparison
}

comparisons <- do.call(rbind, lapply(c("7.9", "7.10"), run_public_comparison))
print(comparisons, row.names = FALSE, digits = 9)
