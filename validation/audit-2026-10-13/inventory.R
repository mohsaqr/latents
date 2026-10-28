out <- "validation/audit-2026-10-13"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source_files <- sort(list.files("R", pattern = "[.]R$", full.names = TRUE))
test_files <- list.files(c("tests/testthat", "tests/equivalence"), pattern = "[.]R$", full.names = TRUE)
test_text <- vapply(test_files, function(p) paste(readLines(p, warn = FALSE), collapse = "\n"), character(1))
test_symbols <- lapply(test_files, function(p) unique(all.names(parse(p), functions = TRUE)))
names(test_symbols) <- test_files
calls_in <- function(e) {
  if (!is.call(e) && !is.pairlist(e) && !is.expression(e)) return(character())
  current <- if (is.call(e)) {
    if (is.symbol(e[[1L]])) as.character(e[[1L]]) else
      if (is.call(e[[1L]]) && identical(e[[1L]][[1L]], as.name("::"))) paste0(as.character(e[[1L]][[2L]]), "::", as.character(e[[1L]][[3L]])) else "<computed call>"
  } else character()
  nested <- lapply(as.list(e), function(part) {
    if (is.symbol(part) && identical(part, quote(expr = ))) character() else calls_in(part)
  })
  unique(c(current, unlist(nested, use.names = FALSE)))
}
namespace <- readLines("NAMESPACE")
exports <- sub("^export[(](.*)[)]$", "\\1", namespace[grepl("^export[(]", namespace)])
methods <- sub("^S3method[(](.*),(.*)[)]$", "\\1.\\2", namespace[grepl("^S3method[(]", namespace)])
rd_files <- sort(list.files("man", pattern = "[.]Rd$", full.names = TRUE))
rd_text <- vapply(rd_files, function(p) paste(readLines(p, warn = FALSE), collapse = "\n"), character(1))
rows <- lapply(source_files, function(file) {
  parsed <- parse(file, keep.source = TRUE)
  text <- readLines(file, warn = FALSE)
  baseline <- suppressWarnings(system2("git", c("show", paste0("HEAD:", file)), stdout = TRUE, stderr = FALSE))
  old <- if (is.null(attr(baseline, "status"))) parse(text = baseline) else expression()
  old_definitions <- Filter(function(e) is.call(e) && is.symbol(e[[1L]]) && identical(e[[1L]], as.name("<-")) && is.symbol(e[[2L]]), as.list(old))
  old_by_name <- stats::setNames(lapply(old_definitions, function(e) e[[3L]]), vapply(old_definitions, function(e) as.character(e[[2L]]), character(1)))
  refs <- attr(parsed, "srcref")
  assignments <- which(vapply(parsed, function(e) is.call(e) && is.symbol(e[[1L]]) && as.character(e[[1L]]) %in% c("<-", "=") && is.symbol(e[[2L]]), logical(1)))
  lapply(assignments, function(i) {
    expr <- parsed[[i]]
    name <- as.character(expr[[2L]])
    rhs <- expr[[3L]]
    is_function <- is.call(rhs) && identical(rhs[[1L]], as.name("function"))
    is_alias <- is.symbol(rhs)
    if (!is_function && !is_alias) return(NULL)
    start <- as.integer(refs[[i]])[1L]
    end <- as.integer(refs[[i]])[3L]
    literal <- paste0("\\alias{", name, "}")
    docs <- names(rd_text)[vapply(rd_text, grepl, logical(1), pattern = literal, fixed = TRUE)]
    needle <- paste0(name, "(")
    test_hits <- names(test_text)[vapply(test_text, grepl, logical(1), pattern = needle, fixed = TRUE)]
    symbols <- names(test_symbols)[vapply(test_symbols, function(s) name %in% s, logical(1))]
    calls <- sort(calls_in(rhs))
    roxy_start <- start - 1L
    preceding <- if (roxy_start > 0L) text[seq_len(roxy_start)] else character()
    last_blank <- max(c(0L, which(!grepl("^#", preceding))))
    roxy <- preceding[seq.int(last_blank + 1L, length(preceding))]
    purpose <- roxy[grepl("^#' [^@]", roxy)]
    purpose <- if (length(purpose)) sub("^#' *", "", purpose[1L]) else ""
    validation <- intersect(calls, c("stopifnot", "match.arg", ".multilpa_check_seed", ".multilpa_check_alignment", ".multilpa_check_row_count", ".multilpa_check_regularity"))
    conditions <- intersect(calls, c("stop", "warning", "message", "tryCatch", "withCallingHandlers", "errorCondition", "warningCondition"))
    random <- calls[grepl("(^|::)(rnorm|runif|rbinom|rpois|rnbinom|sample|sample.int|set.seed)$", calls)]
    io <- calls[grepl("(^|::)(writeLines|write.csv|saveRDS|readRDS|read.csv|system|system2|unlink|file.copy|setwd|options|par|assign)$", calls)]
    body_changed <- !name %in% names(old_by_name) || !identical(deparse(rhs), deparse(old_by_name[[name]]))
    data.frame(file = file, name = name, line = start, end_line = end,
               definition = if (is_function) "function" else "alias_or_constant",
               public = name %in% c(exports, methods),
               documented_in = paste(docs, collapse = ";"),
               direct_test_text_files = paste(test_hits, collapse = ";"),
               direct_test_text_count = length(test_hits),
               purpose = purpose,
               arguments = if (is_function) paste(names(rhs[[2L]]), collapse = ";") else "",
               called_functions = paste(calls, collapse = ";"),
               validation_calls = paste(validation, collapse = ";"),
               condition_calls = paste(conditions, collapse = ";"),
               random_calls = paste(random, collapse = ";"),
               state_or_io_calls = paste(io, collapse = ";"),
               direct_test_symbol_files = paste(symbols, collapse = ";"),
               direct_test_symbol_count = length(symbols),
               body_changed_since_head = body_changed,
               review_status = "parsed and indexed; specialist review in separate reports",
               stringsAsFactors = FALSE)
  })
})
inventory <- do.call(rbind, unlist(rows, recursive = FALSE))
utils::write.csv(inventory, file.path(out, "function-inventory.csv"), row.names = FALSE)
file_rows <- lapply(source_files, function(p) {
  data.frame(file = p, lines = length(readLines(p, warn = FALSE)),
             functions = sum(inventory$file == p & inventory$definition == "function"),
             public_definitions = sum(inventory$file == p & inventory$public),
             md5 = unname(tools::md5sum(p)), stringsAsFactors = FALSE)
})
utils::write.csv(do.call(rbind, file_rows), file.path(out, "source-files.csv"), row.names = FALSE)
doc_files <- sort(unique(c(rd_files, "README.md", "README.Rmd", "NEWS.md", "DESCRIPTION", "NAMESPACE", "_pkgdown.yml", list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE, recursive = TRUE))))
utils::write.csv(data.frame(file = doc_files, lines = vapply(doc_files, function(p) length(readLines(p, warn = FALSE)), integer(1)), md5 = unname(tools::md5sum(doc_files))), file.path(out, "documentation-inventory.csv"), row.names = FALSE)
print(table(inventory$definition, inventory$public))
cat("Source files:", length(source_files), "lines:", sum(vapply(source_files, function(p) length(readLines(p, warn = FALSE)), integer(1))), "\n")
cat("Documentation files:", length(doc_files), "\n")
cat("Public definitions without Rd alias:\n")
print(inventory[inventory$public & inventory$documented_in == "", c("file", "name", "line")], row.names = FALSE)
other_docs <- sort(list.files(c("validation", "equivalence", "future", "data-raw", "tools", "Tutorial_docs"), pattern = "[.](md|MD|Rmd)$", full.names = TRUE, recursive = TRUE))
other_docs <- other_docs[!startsWith(other_docs, paste0(out, "/"))]
other_docs <- unique(c(other_docs, intersect(c("ROADMAP.md", "todo.MD"), list.files())))
utils::write.csv(data.frame(file = other_docs, lines = vapply(other_docs, function(p) length(readLines(p, warn = FALSE)), integer(1)), md5 = unname(tools::md5sum(other_docs)), scope = "historical validation, reference studies, plans or tutorial sources; inventoried, not current public API documentation"), file.path(out, "supporting-documentation-inventory.csv"), row.names = FALSE)
