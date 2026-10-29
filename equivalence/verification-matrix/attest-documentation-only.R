# Re-attest saved comparison results after a documentation-only source change.
#
# A result stays evidence for the current code only when the runtime source it
# was bound to differs from today's in roxygen comments alone. For every
# result sidecar this script requires that
#   1. the sidecar's source hashes equal the committed source at `--base`
#      (the result was current for that commit);
#   2. the result file still has the hash its sidecar recorded;
#   3. today's runtime files are the same set as the base's, NAMESPACE is
#      byte-identical, and every R file parses to identical code
#      (`parse(keep.source = FALSE)`, so comments and layout are ignored).
# Only then is the old sidecar archived and a new one written, labelled a
# retrospective source-equivalence attestation (not a rerun). Nothing is
# refitted. Without --write it only reports.
#
# Usage, from the repository root:
#   Rscript equivalence/verification-matrix/attest-documentation-only.R --base=HEAD
#   Rscript equivalence/verification-matrix/attest-documentation-only.R --base=HEAD --write

source(file.path("equivalence", "verification-matrix", "provenance.R"))

.git_show <- function(base, path, file = TRUE) {
  system2("git", c("show", sprintf("%s:%s", base, path)), stdout = file)
}

.git_md5 <- function(base, path) {
  scratch <- tempfile()
  on.exit(unlink(scratch), add = TRUE)
  if (!identical(.git_show(base, path, scratch), 0L)) return(NA_character_)
  unname(tools::md5sum(scratch))
}

.code_of <- function(text) deparse(parse(text = text, keep.source = FALSE))

# One row per runtime file of the base commit: its hashes there and now, and
# whether its bytes and its parsed code are identical.
.source_check <- function(base, current) {
  listed <- system2("git", c("ls-tree", "-r", "--name-only", base, "R/"), stdout = TRUE)
  paths <- c("NAMESPACE", sort(grep("^R/", listed, value = TRUE), method = "radix"))
  checks <- data.frame(path = paths,
                       base_md5 = vapply(paths, \(p) .git_md5(base, p), character(1)),
                       current_md5 = unname(current[paths]),
                       row.names = NULL, stringsAsFactors = FALSE)
  checks$bytes_identical <- checks$base_md5 == checks$current_md5
  checks$code_identical <- vapply(paths, \(path) {
    base_text <- .git_show(base, path, TRUE)
    if (identical(path, "NAMESPACE")) return(identical(readLines(path), base_text))
    identical(.code_of(readLines(path, warn = FALSE)), .code_of(base_text))
  }, logical(1), USE.NAMES = FALSE)
  checks
}

# Whether the change from the base commit is confined to comments.
.documentation_only <- function(checks, current) {
  namespace <- checks$path == "NAMESPACE"
  c(same_file_set = setequal(checks$path, names(current)),
    namespace_identical = isTRUE(checks$bytes_identical[namespace]),
    code_identical = all(checks$code_identical))
}

# One row per result sidecar: bound to the base source, result bytes unchanged.
.sidecar_check <- function(sidecar, checks) {
  result_path <- sub("[.]provenance[.]csv$", "", sidecar)
  rows <- utils::read.csv(sidecar, stringsAsFactors = FALSE)
  source_rows <- rows[rows$kind == "source", ]
  result_row <- rows[rows$kind == "result", ]
  base_md5 <- stats::setNames(checks$base_md5, checks$path)
  data.frame(
    result = result_path,
    bound_to_base = setequal(source_rows$path, checks$path) &&
      all(source_rows$md5 == base_md5[source_rows$path]),
    result_unchanged = nrow(result_row) == 1L && file.exists(result_path) &&
      identical(result_row$md5, unname(tools::md5sum(result_path))),
    previous_basis = paste(unique(rows$basis), collapse = "+"),
    stringsAsFactors = FALSE
  )
}

# Archive each old sidecar, then bind the unchanged result to the current source.
.reattest <- function(evidence, checks, current, base_commit) {
  stopifnot("every result must be eligible before any is re-attested" = all(evidence$eligible))
  archive <- file.path("validation", sprintf("documentation-only-attestation-%s", Sys.Date()))
  dir.create(file.path(archive, "previous-sidecars"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(checks, file.path(archive, "source-check.csv"), row.names = FALSE)
  utils::write.csv(evidence, file.path(archive, "results.csv"), row.names = FALSE)
  note <- sprintf(paste(
    "Retrospective attestation, not a rerun. The result was bound to the runtime source",
    "of commit %s. The current source differs from it only in roxygen comments: the same",
    "%d files, NAMESPACE byte-identical, every R file identical under parse(keep.source = FALSE)",
    "(%d files with changed comments). Result bytes unchanged. Previous sidecar and checks: %s."),
    base_commit, nrow(checks), sum(!checks$bytes_identical), archive)
  invisible(lapply(evidence$result, \(result_path) {
    sidecar <- paste0(result_path, ".provenance.csv")
    # A distinct suffix, so no scan for live sidecars finds the archived copy.
    kept <- file.path(archive, "previous-sidecars",
                      paste0(result_path, ".provenance.previous.csv"))
    dir.create(dirname(kept), recursive = TRUE, showWarnings = FALSE)
    stopifnot("previous sidecar archived" = file.copy(sidecar, kept, overwrite = FALSE))
    matrix_record_evidence(result_path, current,
                           basis = "retrospective source-equivalence attestation", note = note)
  }))
  stopifnot("every re-attested result is current" =
              all(vapply(evidence$result, matrix_evidence_current, logical(1))))
  archive
}

main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  base <- sub("^--base=", "", grep("^--base=", arguments, value = TRUE))
  base <- if (length(base)) base else "HEAD"
  stopifnot("`--base` must name one commit" = length(base) == 1L && nzchar(base))
  base_commit <- system2("git", c("rev-parse", "--short", base), stdout = TRUE)
  current <- matrix_capture_source(".")
  checks <- .source_check(base, current)
  verdict <- .documentation_only(checks, current)
  cat(sprintf("base %s: %d runtime files; %d with changed bytes, %d with changed code\n",
              base_commit, nrow(checks), sum(!checks$bytes_identical), sum(!checks$code_identical)))
  print(verdict)
  if (!all(verdict)) {
    stop("The runtime source changed beyond comments: results cannot be re-attested; rerun them.",
         call. = FALSE)
  }
  sidecars <- list.files(c("equivalence", "validation"), "[.]provenance[.]csv$",
                         recursive = TRUE, full.names = TRUE)
  evidence <- do.call(rbind, lapply(sidecars, .sidecar_check, checks = checks))
  evidence$eligible <- evidence$bound_to_base & evidence$result_unchanged
  print(evidence, row.names = FALSE)
  cat(sprintf("%d of %d results eligible\n", sum(evidence$eligible), nrow(evidence)))
  if ("--write" %in% arguments) {
    archive <- .reattest(evidence, checks, current, base_commit)
    cat(sprintf("re-attested %d results; archive in %s\n", nrow(evidence), archive))
  }
  invisible(evidence)
}

main()
