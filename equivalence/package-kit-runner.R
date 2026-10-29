# compare.R runs from the project root; testthat sources this from equivalence/.
provenance_path <- file.path("equivalence", "verification-matrix", "provenance.R")
if (!file.exists(provenance_path)) provenance_path <- file.path("verification-matrix", "provenance.R")
source(provenance_path, local = TRUE)

# Shared runner for the package comparison kits (lca-packages/, lpa-packages/,
# ...): mark disagreements that a kit's known-differences.csv explains, fail on
# any other disagreement or on a register entry that explains nothing, and
# write comparison.csv and RUN.txt. Sourced from each kit's compare.R.

#' Validate a comparison before agreement or exception handling
#'
#' Aggregate comparisons deliberately omit individual values, so their finite
#' difference and tolerance are the contract. A reference that could not fit
#' is an explicit failed row with zero values, never a successful comparison.
#' @param comparison Nonempty table of comparison results.
#' @return The input invisibly, or an error for incomplete/inconsistent results.
validate_package_comparison <- function(comparison) {
  required <- c("case", "quantity", "agrees", "difference", "tolerance", "n_values")
  stopifnot("comparison must contain results and required columns" =
              is.data.frame(comparison) && nrow(comparison) > 0L &&
              all(required %in% names(comparison)))
  stopifnot("agreement flags must be logical and complete" =
              is.logical(comparison$agrees) && !anyNA(comparison$agrees),
            "comparison value counts must be finite nonnegative integers" =
              is.numeric(comparison$n_values) && all(is.finite(comparison$n_values)) &&
              all(comparison$n_values >= 0 & comparison$n_values == floor(comparison$n_values)),
            "comparison differences and tolerances must be numeric" =
              is.numeric(comparison$difference) && is.numeric(comparison$tolerance))
  measured <- comparison$n_values > 0L
  stopifnot("measured comparisons require finite nonnegative differences and tolerances" =
              all(is.finite(comparison$difference[measured])) &&
              all(is.finite(comparison$tolerance[measured])) &&
              all(comparison$difference[measured] >= 0 & comparison$tolerance[measured] >= 0),
            "agreement flags must match the numerical comparisons" =
              all(comparison$agrees[measured] ==
                    (comparison$difference[measured] <= comparison$tolerance[measured])),
            "a missing reference solution cannot agree" = !any(comparison$agrees[!measured]))
  invisible(comparison)
}

#' Mark disagreements the register explains; list entries that explain nothing
#' @param comparison Rows with `case`, `quantity`, `agrees` and comparison values.
#' @param known The register: `case`, `quantity` (a regex), `reason`, and optionally
#'   `reference` (required on comparison rows when present in the register).
#' @return A list of the comparison with a `status` column, and the stale rows.
apply_register <- function(comparison, known) {
  validate_package_comparison(comparison)
  stopifnot("known differences require case, quantity and reason" =
              all(c("case", "quantity", "reason") %in% names(known)),
            "reference-specific exceptions require comparison references" =
              !"reference" %in% names(known) || "reference" %in% names(comparison))
  # A different optimum cannot explain an incorrect model dimension. Never
  # let an exception (including a broad legacy regex) waive parameter counts.
  count <- grepl("number of parameters|parameter count|n_parameters", comparison$quantity,
                  ignore.case = TRUE)
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    if (comparison$agrees[[i]] || count[[i]]) return(NA_integer_)
    reference_matches <- if ("reference" %in% names(known)) {
      known$reference == comparison$reference[[i]]
    } else rep(TRUE, nrow(known))
    hits <- which(known$case == comparison$case[[i]] & reference_matches &
                    vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
    if (length(hits)) hits[[1L]] else NA_integer_
  }, integer(1))
  comparison$status <- ifelse(comparison$agrees, "agrees",
                              ifelse(is.na(explained), "disagrees", "known difference"))
  list(comparison = comparison, stale = known[setdiff(seq_len(nrow(known)), explained), ],
       explained = explained)
}

#' Run a kit's comparison, write its files and report
#' @param directory The kit folder.
#' @param compare A function returning the comparison rows.
#' @param packages The reference packages, for RUN.txt.
#' @return `TRUE` when nothing disagrees unexplained and no entry is stale.
run_package_kit <- function(directory, compare, packages) {
  source_snapshot <- matrix_capture_source()
  known <- utils::read.csv(file.path(directory, "known-differences.csv"), stringsAsFactors = FALSE)
  raw <- compare()
  result <- apply_register(raw, known)
  comparison <- result$comparison
  utils::write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
  matrix_record_evidence(file.path(directory, "comparison.csv"), source_snapshot)
  # A kit may record its own search beside the comparison (not pass/fail).
  if (!is.null(attr(raw, "search"))) {
    utils::write.csv(attr(raw, "search"), file.path(directory, "search.csv"), row.names = FALSE)
    print(attr(raw, "search"), row.names = FALSE)
  }
  # ... and further named tables, each written as <name>.csv.
  tables <- attr(raw, "tables") %||% list()
  invisible(lapply(names(tables), function(name) {
    utils::write.csv(tables[[name]], file.path(directory, paste0(name, ".csv")), row.names = FALSE)
    print(tables[[name]], row.names = FALSE)
  }))
  writeLines(c(sprintf("run: %s", format(Sys.Date())),
               sprintf("latents: %s", utils::packageVersion("latents")),
               sprintf("packages: %s", paste(vapply(packages, function(p)
                 sprintf("%s %s", p, utils::packageVersion(p)), character(1)), collapse = ", ")),
               sprintf("R: %s", R.version.string)),
             file.path(directory, "RUN.txt"))
  by_reference <- stats::aggregate(cbind(compared = 1L, agreed = status == "agrees",
                                         known = status == "known difference",
                                         values = n_values) ~ reference, comparison, sum)
  print(by_reference, row.names = FALSE)
  cat(sprintf("\n%d quantities compared: %d agree, %d known differences, %d disagree\n",
              nrow(comparison), sum(comparison$status == "agrees"),
              sum(comparison$status == "known difference"), sum(comparison$status == "disagrees")))
  if (nrow(result$stale)) {
    cat("Registered differences that no longer match a disagreement:\n")
    print(result$stale, row.names = FALSE)
  }
  disagreeing <- subset(comparison, status == "disagrees")
  if (nrow(disagreeing)) print(disagreeing, row.names = FALSE)
  !nrow(disagreeing) && !nrow(result$stale)
}

#' Apply a kit's register to its saved comparison without refitting
#'
#' For a register edited after a run: the compared values and their
#' agreement are those of the saved comparison.csv; only `status` changes.
#' @param directory The kit folder.
#' @return `TRUE` when nothing disagrees unexplained and no entry is stale.
reapply_package_register <- function(directory) {
  # Reclassifying an old result is not a new comparison run. The original
  # sidecar deliberately becomes stale when the result bytes change.
  known <- utils::read.csv(file.path(directory, "known-differences.csv"), stringsAsFactors = FALSE)
  saved <- utils::read.csv(file.path(directory, "comparison.csv"), stringsAsFactors = FALSE)
  result <- apply_register(saved[setdiff(names(saved), "status")], known)
  comparison <- result$comparison
  utils::write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
  cat(sprintf("%d quantities compared: %d agree, %d known differences, %d disagree\n",
              nrow(comparison), sum(comparison$status == "agrees"),
              sum(comparison$status == "known difference"), sum(comparison$status == "disagrees")))
  if (nrow(result$stale)) print(result$stale, row.names = FALSE)
  !any(comparison$status == "disagrees") && !nrow(result$stale)
}
