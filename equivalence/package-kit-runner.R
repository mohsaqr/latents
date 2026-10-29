# Shared runner for the package comparison kits (lca-packages/, lpa-packages/,
# ...): mark disagreements that a kit's known-differences.csv explains, fail on
# any other disagreement or on a register entry that explains nothing, and
# write comparison.csv and RUN.txt. Sourced from each kit's compare.R.

#' Mark disagreements the register explains; list entries that explain nothing
#' @param comparison Rows with `case`, `reference`, `quantity`, `agrees`.
#' @param known The register: `case`, `reference`, `quantity` (a regex), `reason`.
#' @return A list of the comparison with a `status` column, and the stale rows.
apply_register <- function(comparison, known) {
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    hits <- which(known$case == comparison$case[[i]] & known$reference == comparison$reference[[i]] &
                    vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
    if (length(hits) && !comparison$agrees[[i]]) hits[[1L]] else NA_integer_
  }, integer(1))
  comparison$status <- ifelse(comparison$agrees, "agrees",
                              ifelse(is.na(explained), "disagrees", "known difference"))
  list(comparison = comparison, stale = known[setdiff(seq_len(nrow(known)), explained), ])
}

#' Run a kit's comparison, write its files and report
#' @param directory The kit folder.
#' @param compare A function returning the comparison rows.
#' @param packages The reference packages, for RUN.txt.
#' @return `TRUE` when nothing disagrees unexplained and no entry is stale.
run_package_kit <- function(directory, compare, packages) {
  known <- utils::read.csv(file.path(directory, "known-differences.csv"), stringsAsFactors = FALSE)
  raw <- compare()
  result <- apply_register(raw, known)
  comparison <- result$comparison
  utils::write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
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
