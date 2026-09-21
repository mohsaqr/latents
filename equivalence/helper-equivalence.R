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

# testthat runs every file with this folder as the working directory, so
# paths are relative to equivalence/. The in-house exhaustive-enumeration
# reference is shared with the shipped tests, so it is read from there rather
# than copied.
source(file.path("..", "tests", "testthat", "helper-independent-likelihood.R"),
       local = TRUE)
