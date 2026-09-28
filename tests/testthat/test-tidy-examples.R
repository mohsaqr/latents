# Every help-page example is written as verbs with named arguments: no `$`,
# no bracket indexing and no `order()` on the public surface. An example
# teaches the idiom a user copies, so one that subsets a result, or builds its
# input by hand-editing columns, teaches the ritual the package exists to
# remove. This reads the examples as users see them, from the help pages.

.example_code <- function() {
  source_man <- test_path("..", "..", "man")
  database <- if (dir.exists(source_man)) {
    tools::Rd_db(dir = normalizePath(test_path("..", "..")))
  } else tools::Rd_db("latents")
  unlist(lapply(names(database), function(topic) {
    code <- utils::capture.output(tools::Rd2ex(database[[topic]], out = stdout()))
    code <- code[!grepl("^\\s*#", code) & nzchar(trimws(code))]
    if (length(code) == 0L) return(character())
    stats::setNames(code, rep(topic, length(code)))
  }))
}

test_that("help-page examples use verbs, not subsetting", {
  code <- .example_code()
  expect_gt(length(code), 100L)
  subsetting <- grepl("\\$[A-Za-z_.]", code) |
    grepl("[A-Za-z0-9_.)]\\[", code) |
    grepl("\\border\\(", code)
  offenders <- paste(names(code)[subsetting], code[subsetting], sep = ": ")
  expect_identical(offenders, character())
})
