# Where a validation script reads and writes an Mplus reference fixture.
#
# Fixtures whose data the shipped unit tests read stay in `tests/fixtures/`;
# fixtures that only the equivalence tests read live in `equivalence/fixtures/`
# (see `equivalence/README.md`). Commit 504ae1f made that split, and the scripts
# that write the fixtures kept writing to the old place, so a rerun recreated
# stale copies under `tests/fixtures/` that nothing read. The path is decided in
# one place here instead of at each call site.

#' Path to one Mplus reference fixture
#'
#' @param name File name of the fixture, e.g. `"public-7.9.rds"`.
#' @param root Project root directory.
#' @return The path under `tests/fixtures/mplus/` when the shipped tests hold
#'   that fixture, otherwise under `equivalence/fixtures/mplus/`.
mplus_fixture <- function(name, root = ".") {
  stopifnot("`name` must be a single file name" =
              is.character(name) && length(name) == 1L && !is.na(name))
  shipped <- file.path(root, "tests", "fixtures", "mplus", name)
  if (file.exists(shipped)) return(shipped)
  directory <- file.path(root, "equivalence", "fixtures", "mplus")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, name)
}
