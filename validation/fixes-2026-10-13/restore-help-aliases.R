old_files <- system2("git", c("ls-tree", "--name-only", "a10a919:man"), stdout = TRUE)
old_files <- old_files[grepl("[.]Rd$", old_files)]
aliases <- unique(unlist(lapply(old_files, function(p) {
  text <- system2("git", c("show", paste0("a10a919:man/", p)), stdout = TRUE)
  doc <- tools::parse_Rd(textConnection(text))
  blocks <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"), doc)
  vapply(blocks, function(x) paste(x, collapse = ""), character(1))
})))
current_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
current_aliases <- unique(unlist(lapply(current_files, function(p) {
  doc <- tools::parse_Rd(p)
  blocks <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"), doc)
  vapply(blocks, function(x) paste(x, collapse = ""), character(1))
})))
restored <- sort(setdiff(aliases[grepl("^(get_results|print|summary|as[.]data[.]frame|coef|vcov|confint|logLik|nobs)[.]", aliases)], current_aliases))
source <- readLines("R/get-data.R")
marker <- which(grepl("^get_results <- function", source))
stopifnot(length(marker) == 1L)
insert <- marker - 1L
stopifnot(source[insert] == "#' @export")
addition <- sprintf("#' @aliases %s", restored)
writeLines(append(source, addition, after = insert - 1L), "R/get-data.R")
writeLines(restored, "validation/fixes-2026-10-13/restored-help-aliases.txt")
cat(length(restored), "released method help aliases restored to the shared results topic\n")
