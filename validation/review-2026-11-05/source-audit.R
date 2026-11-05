args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
scratch <- args[[1L]]
out <- args[[2L]]
#' Parse top-level named function definitions without comments
#' @param root A package source directory.
#' @return A named list of file names, formals and bodies.
functions_in <- function(root) {
  stopifnot(is.character(root), length(root) == 1L, dir.exists(file.path(root, "R")))
  files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  result <- list()
  lapply(files, function(file) {
    expressions <- parse(file)
    lapply(expressions, function(expr) {
      if (is.call(expr) && is.symbol(expr[[1L]]) &&
          as.character(expr[[1L]]) %in% c("<-", "=") &&
          length(expr) == 3L && is.symbol(expr[[2L]]) &&
          is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))) {
        result[[as.character(expr[[2L]])]] <<- list(
          file = basename(file), formals = expr[[3L]][[2L]], body = expr[[3L]][[3L]])
      }
    })
  })
  result
}
#' Read explicit exports from a package namespace
#' @param root A package source directory.
#' @return A character vector of exported names.
exports_in <- function(root) {
  stopifnot(is.character(root), length(root) == 1L,
            file.exists(file.path(root, "NAMESPACE")))
  text <- readLines(file.path(root, "NAMESPACE"))
  sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", text, value = TRUE))
}
current <- functions_in(file.path(scratch, "source"))
current_exports <- exports_in(file.path(scratch, "source"))
baseline_results <- lapply(c("a10a919b", "8518dffe", "HEAD"), function(base) {
  old <- functions_in(file.path(scratch, base))
  old_exports <- exports_in(file.path(scratch, base))
  common <- intersect(names(old), names(current))
  changed <- common[!vapply(common, function(name)
    identical(old[[name]]$body, current[[name]]$body) &&
      identical(old[[name]]$formals, current[[name]]$formals), logical(1L))]
  formal_changes <- common[!vapply(common, function(name)
    identical(old[[name]]$formals, current[[name]]$formals), logical(1L))]
  list(baseline = base, functions_before = length(old), functions_current = length(current),
    added = setdiff(names(current), names(old)), removed = setdiff(names(old), names(current)),
    changed = changed, exported_removed = setdiff(old_exports, current_exports),
    exported_added = setdiff(current_exports, old_exports),
    changed_export_formals = lapply(intersect(formal_changes, old_exports), function(name)
      list(name = name, before = old[[name]]$formals, after = current[[name]]$formals)))
})
saveRDS(baseline_results, file.path(out, "source-audit.rds"))
lapply(baseline_results, function(x) {
  cat("\nBASELINE", x$baseline, "\n")
  cat("functions before/current:", x$functions_before, x$functions_current, "\n")
  cat("added/removed/changed:", length(x$added), length(x$removed), length(x$changed), "\n")
  cat("exports removed:", paste(x$exported_removed, collapse = ", "), "\n")
  cat("exports added:", paste(x$exported_added, collapse = ", "), "\n")
  print(x$changed_export_formals)
  writeLines(x$changed, file.path(out, paste0("changed-functions-", x$baseline, ".txt")))
})
