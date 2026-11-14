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
#'   `reference` (required on comparison rows when present in the register) and
#'   `occurs`: `"always"` (the default) or `"sometimes"`, for a reference whose
#'   result differs by platform or version. A `"sometimes"` entry explains a
#'   difference when one occurs and is not stale when none does.
#' @return A list of the comparison with a `status` column, and the stale rows.
apply_register <- function(comparison, known) {
  validate_package_comparison(comparison)
  stopifnot("known differences require case, quantity and reason" =
              all(c("case", "quantity", "reason") %in% names(known)),
            "reference-specific exceptions require comparison references" =
              !"reference" %in% names(known) || "reference" %in% names(comparison),
            "`occurs` must be \"always\" or \"sometimes\"" =
              !"occurs" %in% names(known) || all(known$occurs %in% c("always", "sometimes")))
  optional <- if ("occurs" %in% names(known)) known$occurs == "sometimes" else
    rep(FALSE, nrow(known))
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
  list(comparison = comparison, stale = known[setdiff(which(!optional), explained), ],
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

#' The reference group of each comparison row
#'
#' "R packages", "published analyses", or the commercial program a row's
#' reference names ("Mplus", "Latent GOLD", "Stata", "SAS"); NA for a row of another
#' kind. One group is one reference axis of CVS and the verification matrix.
#' @param comparison Comparison rows; their `kind` column, or `kind` for all.
#' @param kind Kind of every row when the comparison has no `kind` column.
#' @return Character vector, one group per row.
kit_reference_group <- function(comparison, kind = NULL) {
  kinds <- comparison$kind %||% rep(kind %||% NA_character_, nrow(comparison))
  programs <- c("Latent GOLD", "Mplus", "Stata", "SAS")
  vapply(seq_len(nrow(comparison)), function(i) {
    switch(kinds[[i]],
      `R package` = "R packages",
      published = "published analyses",
      `commercial program` = {
        named <- programs[vapply(programs, grepl, logical(1), x = comparison$reference[[i]],
                                 fixed = TRUE)]
        if (length(named) == 1L) named else NA_character_
      },
      NA_character_)
  }, character(1))
}

#' The testthat blocks of a package kit, one per reference group
#'
#' The first block checks that the kit compared every case and that its
#' register explains only what occurs; then one block per reference group
#' (`kit_reference_group()`), so each block is evidence for one reference axis:
#' Save a kit's comparison, computed by a test, as verification-matrix evidence
#'
#' Writes tmp/kits/<kit>/comparison.csv under the project root with its
#' provenance sidecar (verification-matrix/provenance.R), so the comparison a
#' ci run computes from the current package becomes current evidence once it is
#' copied over the kit's comparison.csv (equivalence/README.md, "Evidence for
#' the matrix"). Call from a testthat file in equivalence/.
#' @param kit Kit folder, e.g. "latentgold-wider".
#' @param comparison The comparison rows, in the layout of the kit's own
#'   comparison.csv.
#' @param snapshot `matrix_capture_source("..")`, taken before the comparison
#'   was computed.
#' @return The path written, invisibly.
kit_save_comparison <- function(kit, comparison, snapshot) {
  stopifnot("`kit` is one kit folder" = is.character(kit) && length(kit) == 1L && dir.exists(kit),
            "`comparison` has rows" = is.data.frame(comparison) && nrow(comparison) > 0L)
  output <- file.path("..", "tmp", "kits", kit)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output, "comparison.csv")
  utils::write.csv(comparison, path, row.names = FALSE)
  matrix_record_evidence(path, snapshot, "..")
  invisible(path)
}

#' CVS and the verification matrix take one reference per test, and a block
#' checking a package and a published table together counted for neither.
#' The comparison is computed once, by the first block that needs it, and
#' written with its provenance sidecar (verification-matrix/provenance.R) to
#' tmp/kits/<kit>/comparison.csv, where the verification matrix takes it from.
#' Call from a testthat file in equivalence/ (the working directory testthat
#' sets).
#' @param kit Kit folder, e.g. "stata-published"; also the start of each title.
#' @param compute Function with no arguments returning the comparison rows.
#' @param n_rows Rows the comparison must have.
#' @param groups The reference groups the kit must have.
#' @param kind Kind of every row when the comparison has no `kind` column.
#' @param coverage Optional function of the registered result (`comparison`,
#'   `stale`, `explained`, `known`) adding expectations to the first block.
#' @return Nothing useful; the blocks are registered with testthat.
kit_reference_tests <- function(kit, compute, n_rows, groups, kind = NULL, coverage = NULL) {
  stopifnot("`kit` is one kit folder" = is.character(kit) && length(kit) == 1L && dir.exists(kit),
            "`compute` is a function" = is.function(compute),
            "`groups` names reference groups" = is.character(groups) && length(groups) > 0L)
  state <- new.env()
  registered <- function() {
    if (is.null(state$value)) {
      snapshot <- matrix_capture_source("..")
      result <- kit_registered_comparison(kit, compute, kind)
      kit_save_comparison(kit, result$comparison[setdiff(names(result$comparison), "group")],
                          snapshot)
      state$value <- result
    }
    state$value
  }
  testthat::test_that(sprintf("%s compares every case, and its register explains only what occurs",
                              kit), {
    kit_expect_complete(registered(), n_rows, groups, coverage)
  })
  invisible(lapply(groups, function(group) {
    label <- if (group %in% c("R packages", "published analyses")) paste("the", group) else group
    testthat::test_that(sprintf("%s agrees with %s", kit, label), {
      kit_expect_group(registered(), group)
    })
  }))
}

#' A kit's comparison with its register applied and each row's reference group
#' @param kit Kit folder, holding known-differences.csv.
#' @param compute,kind As for [kit_reference_tests()].
#' @return The [apply_register()] result (`comparison` with a `group` column,
#'   `stale`, `explained`) plus the register as `known`.
kit_registered_comparison <- function(kit, compute, kind = NULL) {
  known <- utils::read.csv(file.path(kit, "known-differences.csv"), stringsAsFactors = FALSE)
  result <- apply_register(compute(), known)
  result$comparison$group <- kit_reference_group(result$comparison, kind)
  c(result, list(known = known))
}

#' Expectations of a kit's first block: every case compared, every group
#' present, no stale register entry, and the kit's own coverage floors
#' @param result From [kit_registered_comparison()].
#' @param n_rows,groups,coverage As for [kit_reference_tests()].
#' @return Nothing useful; called for its expectations.
kit_expect_complete <- function(result, n_rows, groups, coverage = NULL) {
  testthat::expect_identical(nrow(result$comparison), n_rows)
  testthat::expect_setequal(unique(result$comparison$group), groups)
  testthat::expect_identical(nrow(result$stale), 0L)
  if (!is.null(coverage)) coverage(result)
  invisible(NULL)
}

#' Expectations of one reference group: it has rows and none disagrees unexplained
#' @param result From [kit_registered_comparison()].
#' @param group One reference group ([kit_reference_group()]).
#' @return Nothing useful; called for its expectations.
kit_expect_group <- function(result, group) {
  rows <- result$comparison[result$comparison$group %in% group, , drop = FALSE]
  testthat::expect_gt(nrow(rows), 0L)
  unexplained <- rows[rows$status == "disagrees", , drop = FALSE]
  testthat::expect_identical(nrow(unexplained), 0L,
                             info = paste(unexplained$case, unexplained$reference,
                                          unexplained$quantity, collapse = "; "))
  invisible(NULL)
}
