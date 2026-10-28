S <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/eq"
# Drop printed vector indices and the new gaussian_diagnostics block.
strip <- function(x) {
  x <- sub("^ *\\[[0-9]+\\] *", "", x)
  s <- grep("gaussian_diagnostics -", x)
  if (length(s)) {
    e <- grep("^\"-- ", x)
    e <- e[e > s[1]][1] - 1L
    x <- x[-(s[1]:e)]
  }
  x
}
for (case in c("svc_diagnostics_report", "profiles_lpa_single_varying")) {
  a <- readRDS(file.path(S, "rec-base", "golden", paste0(case, ".rds")))$services$report$value
  b <- readRDS(file.path(S, "rec-cur", "golden", paste0(case, ".rds")))$services$report$value
  cat(case, ": print identical after removing indices and the new table:",
      identical(strip(a$print$value), strip(b$print$value)), "\n")
  va <- a$value$value; vb <- b$value$value
  cat("   report value: names only in current:", paste(setdiff(names(vb), names(va)), collapse = ","),
      "; shared parts identical:", identical(va, vb[names(va)]), "\n")
}
a <- readRDS(file.path(S, "rec-base", "golden", "svc_diagnostics_report.rds"))$services$report$value
b <- readRDS(file.path(S, "rec-cur", "golden", "svc_diagnostics_report.rds"))$services$report$value
cat("\nstr(base report):\n"); str(a, max.level = 3, list.len = 8, vec.len = 1)
pa <- strip(a$print$value); pb <- strip(b$print$value)
cat("\nlengths after strip:", length(pa), length(pb), "\n")
cat("only base:\n"); print(head(setdiff(pa, pb), 10)); cat("only current:\n"); print(head(setdiff(pb, pa), 10))
cat("\nvalue diff:\n"); print(all.equal(a$value, b$value))
