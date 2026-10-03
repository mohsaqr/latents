# Shared set-up for record.R and compare.R: load the package from source, pin
# the session settings a fingerprint depends on, and load the catalogue.

#' Directory of the running script (`Rscript path/to/script.R`)
#' @return An absolute path.
golden_script_dir <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                        value = TRUE)
  if (length(file_argument) == 0L) {
    return(normalizePath(file.path("equivalence", "engine-golden")))
  }
  dirname(normalizePath(sub("^--file=", "", file_argument[1L])))
}

#' Load the package and the catalogue into the global environment
#'
#' Printed tables depend on `width`, `digits` and `scipen`, and any sort of
#' labels on the collation locale; all are pinned so a fingerprint does not
#' change with the caller's profile. Plots are sent to a null device.
#' @param dir The engine-golden directory.
#' @return The package root, invisibly.
golden_setup <- function(dir) {
  root <- normalizePath(file.path(dir, "..", ".."))
  suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))
  options(width = 80L, digits = 7L, scipen = 0L, OutDec = ".",
          warn = 1L, useFancyQuotes = FALSE)
  Sys.setlocale("LC_COLLATE", "C")
  grDevices::pdf(NULL)
  sources <- c(file.path(dir, c("fingerprint.R", "catalogue.R")),
               sort(list.files(dir, pattern = "^catalogue-.*[.]R$",
                               full.names = TRUE)))
  invisible(lapply(sources, sys.source, envir = globalenv()))
  invisible(root)
}

#' Fit one case, time the fit, and fingerprint the result
#'
#' The fit is timed alone. A fit under one second is run twice and the
#' shorter time kept, because timing that short is mostly noise; the
#' fingerprint is taken from the first run.
#' @param case A catalogue case.
#' @return `list(fingerprint =, seconds =)`.
golden_run_case <- function(case) {
  seed <- golden_case_seed(case$name)
  timed <- function() {
    set.seed(seed)
    started <- proc.time()[["elapsed"]]
    object <- golden_probe(case$fit(), clean = FALSE)
    list(object = object, seconds = proc.time()[["elapsed"]] - started)
  }
  first <- timed()
  seconds <- first$seconds
  if (seconds < 1) seconds <- min(seconds, timed()$seconds)
  list(fingerprint = golden_fingerprint(first$object, case, seed),
       seconds = seconds)
}

#' Replace the rows of `table` for `cases` with `rows`, keeping the others
#' @param path CSV path (may not exist).
#' @param rows New rows, with a `case` column.
#' @return The merged table, written to `path`, invisibly.
golden_merge_csv <- function(path, rows) {
  old <- if (file.exists(path)) {
    utils::read.csv(path, stringsAsFactors = FALSE)
  } else {
    rows[0L, , drop = FALSE]
  }
  if (nrow(old) > 0L) old <- old[!old$case %in% rows$case, , drop = FALSE]
  merged <- rbind(old[names(rows)], rows)
  utils::write.csv(merged, path, row.names = FALSE)
  invisible(merged)
}

#' Output directory: GOLDEN_OUT if set (for scratch runs), else `dir`
#' @param dir The engine-golden directory.
#' @return A path.
golden_output_dir <- function(dir) {
  out <- Sys.getenv("GOLDEN_OUT", "")
  if (nzchar(out)) out else dir
}
