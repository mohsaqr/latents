# Evidence is current only while both its result and the complete runtime R
# source manifest match. MD5 is used for accidental-change detection, not as a
# security signature. Capture before a run; record only after writing results.

#' Capture every runtime R file and the namespace before a comparison
#' @param root Package source directory.
#' @return Named MD5 hashes, with the capture time and package version attached.
matrix_capture_source <- function(root = ".") {
  paths <- c("NAMESPACE", paste0("R/", sort(list.files(file.path(root, "R"), recursive = TRUE,
                                                     full.names = FALSE), method = "radix")))
  stopifnot(length(paths) > 1L, all(file.exists(file.path(root, paths))))
  hashes <- stats::setNames(unname(tools::md5sum(file.path(root, paths))), paths)
  stopifnot(!anyNA(hashes))
  attr(hashes, "captured_at") <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  attr(hashes, "version") <- read.dcf(file.path(root, "DESCRIPTION"), fields = "Version")[[1L]]
  hashes
}

#' Record a completed run, or an explicitly documented historical attestation
#' @param path Written result CSV.
#' @param snapshot Source hashes captured before the run.
#' @param root Package source directory.
#' @param basis Either a run or a retrospective source-equivalence attestation.
#' @param note Evidence supporting an attestation; never imply it is a rerun.
#' @return The sidecar path, invisibly.
matrix_record_evidence <- function(path, snapshot, root = ".", basis = "recorded run",
                                   note = "Source captured before this comparison run.") {
  stopifnot(is.character(snapshot), length(snapshot) > 1L, file.exists(path),
            basis %in% c("recorded run", "retrospective source-equivalence attestation"),
            is.character(note), length(note) == 1L, nzchar(note))
  current <- matrix_capture_source(root)
  stopifnot("runtime source changed during the comparison" =
              !anyDuplicated(names(snapshot)) && setequal(names(snapshot), names(current)) &&
              all(snapshot[names(current)] == current))
  rows <- data.frame(kind = c(rep("source", length(snapshot)), "result"),
                     path = c(names(snapshot), basename(path)),
                     md5 = c(unname(snapshot), unname(tools::md5sum(path))),
                     basis = basis, recorded_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                     captured_at = attr(snapshot, "captured_at"),
                     version = attr(snapshot, "version"), note = note,
                     stringsAsFactors = FALSE)
  sidecar <- paste0(path, ".provenance.csv")
  utils::write.csv(rows, sidecar, row.names = FALSE)
  invisible(sidecar)
}

#' Validate a result against its source and result digests
#' @param path Result CSV.
#' @param root Package source directory.
#' @return TRUE only for a valid result with an unchanged complete source tree.
matrix_evidence_current <- function(path, root = ".") {
  sidecar <- paste0(path, ".provenance.csv")
  if (!file.exists(path) || !file.exists(sidecar)) return(FALSE)
  tryCatch({
    rows <- utils::read.csv(sidecar, stringsAsFactors = FALSE)
    if (!all(c("kind", "path", "md5", "basis", "note") %in% names(rows))) return(FALSE)
    if (anyNA(rows) || any(!rows$basis %in%
        c("recorded run", "retrospective source-equivalence attestation"))) return(FALSE)
    source <- rows[rows$kind == "source", ]
    result <- rows[rows$kind == "result", ]
    current <- matrix_capture_source(root)
    isTRUE(nrow(result) == 1L && nrow(rows) == nrow(source) + 1L &&
             !anyDuplicated(source$path) && setequal(source$path, names(current)) &&
             all(source$md5[match(names(current), source$path)] == current) &&
             identical(result$path, basename(path)) &&
             identical(result$md5, unname(tools::md5sum(path))))
  }, error = function(e) FALSE)
}
