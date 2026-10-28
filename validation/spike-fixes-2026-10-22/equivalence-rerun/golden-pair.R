# Compare golden fingerprints recorded under 12e0145 and under the working tree.
S <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/eq"
# Leaf paths at which two nested objects differ (doubles within 1e-10 relative).
diff_paths <- function(a, b, path = "") {
  if (is.list(a) && is.list(b) && !is.data.frame(a)) {
    keys <- union(names(a) %||% seq_along(a), names(b) %||% seq_along(b))
    return(unlist(lapply(keys, \(k) {
      ka <- if (is.character(k)) a[[k]] else a[[as.integer(k)]]
      kb <- if (is.character(k)) b[[k]] else b[[as.integer(k)]]
      if (is.null(ka) || is.null(kb)) return(sprintf("%s/%s [%s]", path, k,
        if (is.null(ka)) "only current" else "only base"))
      diff_paths(ka, kb, paste0(path, "/", k))
    })))
  }
  if (identical(a, b)) return(character())
  if (is.numeric(a) && is.numeric(b) && length(a) == length(b) &&
      isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))) return(character())
  if (is.data.frame(a) && is.data.frame(b)) {
    extra <- setdiff(names(b), names(a)); gone <- setdiff(names(a), names(b))
    common <- intersect(names(a), names(b))
    inner <- if (nrow(a) == nrow(b)) unlist(lapply(common, \(k) diff_paths(a[[k]], b[[k]], paste0(path, "$", k)))) else
      sprintf("%s [rows %d -> %d]", path, nrow(a), nrow(b))
    return(c(if (length(extra)) sprintf("%s [new columns: %s]", path, paste(extra, collapse = ",")),
             if (length(gone)) sprintf("%s [removed columns: %s]", path, paste(gone, collapse = ",")), inner))
  }
  path
}
cases <- sub("[.]rds$", "", list.files(file.path(S, "rec-cur", "golden")))
res <- lapply(cases, \(case) {
  a <- readRDS(file.path(S, "rec-base", "golden", paste0(case, ".rds")))
  b <- readRDS(file.path(S, "rec-cur", "golden", paste0(case, ".rds")))
  if (identical(a, b)) return(data.frame(case = case, path = NA_character_))
  p <- diff_paths(a, b)
  data.frame(case = case, path = if (length(p)) p else NA_character_)
})
all_rows <- do.call(rbind, res)
saveRDS(all_rows, file.path(S, "golden-pair.rds"))
changed <- subset(all_rows, !is.na(path))
cat("cases:", length(cases), " identical:", length(cases) - length(unique(changed$case)),
    " differing:", length(unique(changed$case)), "\n\n")
# Generalise the paths so the kinds of change can be counted.
kind <- gsub("^[^/]+/", "", sub("^/", "", changed$path))
kind <- gsub("/[0-9]+(/|$)", "/#\\1", kind)
print(sort(table(kind), decreasing = TRUE))
