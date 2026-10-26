# Shared plumbing for the equivalence suite. Nothing in this folder ships with
# the package: it is listed in .Rbuildignore and run by equivalence/run.R.

# A recorded reference lives in exactly one of two stores. Fixtures that the
# shipped unit tests also read as input data stay in tests/fixtures/; fixtures
# that only an equivalence test reads live in equivalence/fixtures/. A name
# found in neither or in both is an error, never a silent pick.
equivalence_fixture <- function(...) {
  stores <- c(file.path("fixtures", ...),
              file.path("..", "tests", "fixtures", ...))
  found <- stores[file.exists(stores)]
  if (length(found) != 1L) {
    stop(sprintf("fixture %s found in %d stores; expected exactly one",
                 file.path(...), length(found)), call. = FALSE)
  }
  found
}

# Save fresh diagnostics only when the verification runner declares an output
# directory. Normal tests never rewrite retained reference reports.
equivalence_write_tables <- function(family, tables) {
  output <- getOption("latents.equivalence.output_dir", NULL)
  if (is.null(output)) return(invisible(NULL))
  stopifnot(is.character(output), length(output) == 1L, !is.na(output), nzchar(output),
            is.character(family), length(family) == 1L, nzchar(family),
            is.list(tables), !is.null(names(tables)),
            all(vapply(tables, is.data.frame, logical(1))))
  directory <- file.path(output, family)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  invisible(lapply(names(tables), function(name) {
    utils::write.csv(tables[[name]], file.path(directory, paste0(name, ".csv")),
                     row.names = FALSE)
  }))
}

# testthat runs every file with this folder as the working directory, so
# paths are relative to equivalence/. The in-house exhaustive-enumeration
# reference is shared with the in-package equivalence tests, so it is read
# from there rather than copied (it moved out of tests/testthat in 0.9.19).
source(file.path("..", "tests", "equivalence", "helper-independent-likelihood.R"),
       local = TRUE)
