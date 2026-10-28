S <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/eq"
rows <- subset(readRDS(file.path(S, "golden-pair.rds")), !is.na(path))
get_path <- function(x, path) Reduce(\(o, k) if (grepl("^[0-9]+$", k)) o[[as.integer(k)]] else o[[k]],
                                     strsplit(sub("^/", "", path), "/")[[1]], x)
# Text outputs: drop the new eigenvalue line and the gaussian_diagnostics block.
clean <- function(txt) {
  txt <- grep("Smallest standardized Gaussian eigenvalue|Gaussian covariance ratios unavailable", txt, value = TRUE, invert = TRUE)
  start <- grep("gaussian_diagnostics", txt)
  if (length(start)) {
    nxt <- grep("^-- ", txt); nxt <- nxt[nxt > start[1]]
    stop_at <- if (length(nxt)) nxt[1] - 1L else length(txt)
    txt <- txt[-(start[1]:stop_at)]
  }
  # Table counts in summary footers ("N tables above") move by the new table.
  txt <- gsub("[0-9]+ tables above", "N tables above", txt)
  txt
}
text_rows <- subset(rows, grepl("(summary|print)/value$", path))
leftover <- Filter(Negate(is.null), lapply(seq_len(nrow(text_rows)), \(i) {
  r <- text_rows[i, ]
  a <- get_path(readRDS(file.path(S, "rec-base", "golden", paste0(r$case, ".rds"))), r$path)
  b <- get_path(readRDS(file.path(S, "rec-cur", "golden", paste0(r$case, ".rds"))), r$path)
  if (identical(clean(a), clean(b))) NULL else list(case = r$case, path = r$path,
    base_only = setdiff(clean(a), clean(b)), cur_only = setdiff(clean(b), clean(a)))
}))
cat("print/summary diffs:", nrow(text_rows), "; still different after removing the additive lines:", length(leftover), "\n")
for (l in head(leftover, 6)) { cat("\n#", l$case, l$path, "\n- base:\n"); writeLines(head(l$base_only, 6)); cat("+ cur:\n"); writeLines(head(l$cur_only, 6)) }
msg_rows <- subset(rows, grepl("message$", path))
msgs <- unique(do.call(rbind, lapply(seq_len(nrow(msg_rows)), \(i) {
  r <- msg_rows[i, ]
  data.frame(base = paste(get_path(readRDS(file.path(S, "rec-base", "golden", paste0(r$case, ".rds"))), r$path), collapse = " "),
             cur = paste(get_path(readRDS(file.path(S, "rec-cur", "golden", paste0(r$case, ".rds"))), r$path), collapse = " "))
})))
cat("\nDistinct message changes:", nrow(msgs), "\n")
for (i in seq_len(nrow(msgs))) cat(sprintf("\n- base: %s\n+ cur:  %s\n", strtrim(msgs$base[i], 230), strtrim(msgs$cur[i], 230)))
