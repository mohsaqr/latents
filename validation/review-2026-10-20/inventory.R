destination <- "validation/review-2026-10-20"
files <- list.files("R", "[.]R$", full.names = TRUE)
inventory <- do.call(rbind, lapply(files, function(path) {
  expressions <- parse(path, keep.source = TRUE)
  references <- attr(expressions, "srcref")
  definitions <- which(vapply(expressions, function(e) {
    is.call(e) && identical(e[[1L]], as.name("<-")) && length(e) == 3L &&
      is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))
  }, logical(1)))
  if (!length(definitions)) return(NULL)
  do.call(rbind, lapply(definitions, function(i) {
    expression <- expressions[[i]]
    data.frame(file = path, name = as.character(expression[[2L]]),
      start = references[[i]][1L], end = references[[i]][3L],
      arguments = paste(names(expression[[3L]][[2L]]), collapse = ","))
  }))
}))
write.csv(inventory, file.path(destination, "function-inventory.csv"), row.names = FALSE)
manifest <- data.frame(file = files,
  lines = vapply(files, function(f) length(readLines(f, warn = FALSE)), integer(1)),
  md5 = unname(tools::md5sum(files)))
write.csv(manifest, file.path(destination, "source-manifest.csv"), row.names = FALSE)
namespace <- readLines("NAMESPACE")
unit <- list.files("tests/testthat", "^test-.*[.]R$", full.names = TRUE)
equivalence <- list.files("tests/equivalence", "^test-.*[.]R$", full.names = TRUE)
cat(sprintf("Production files: %d\nProduction lines: %d\nTop-level function definitions: %d\nExports: %d\nRegistered S3 methods: %d\nUnit test files: %d\nEquivalence test files: %d\n",
  nrow(manifest), sum(manifest$lines), nrow(inventory),
  sum(startsWith(namespace, "export(")), sum(startsWith(namespace, "S3method(")),
  length(unit), length(equivalence)))
print(sessionInfo())
