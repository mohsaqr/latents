S <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad"
stems <- setdiff(list.dirs(file.path(S, "out", "cur"), recursive = FALSE, full.names = FALSE),
                 "class-enumeration-orig")
tokens <- \(w, s) {
  out <- grep("^#>", readLines(file.path(S, "out", w, s, paste0(s, ".md"))), value = TRUE)
  table(unlist(regmatches(out, gregexpr("-?[0-9]+\\.?[0-9]*(e[-+]?[0-9]+)?", out))))
}
rows <- lapply(stems, \(s) {
  b <- tokens("base", s); c <- tokens("cur", s)
  have <- c[names(b)]; have[is.na(have)] <- 0
  lost <- names(b)[have < b]
  data.frame(vignette = s, base_numbers = sum(b), cur_numbers = sum(c),
             base_numbers_missing_in_cur = length(lost),
             examples = paste(head(lost, 8), collapse = " "))
})
print(do.call(rbind, rows), row.names = FALSE)
