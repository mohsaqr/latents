# Compare baseline vs current knits of every vignette: status, saved
# objects (after stripping the purely additive diagnostics), markdown lines
# and figure files.
S <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad"
new_columns <- c("min_standardized_eigenvalue", "min_gaussian_mass")

# Remove the additive fields so that what remains is the pre-existing result.
strip_new <- function(o) {
  if (is.data.frame(o)) return(o[, setdiff(names(o), new_columns), drop = FALSE])
  if (is.list(o) && !is.null(names(o))) {
    o <- o[setdiff(names(o), "gaussian_diagnostics")]
    o[] <- lapply(o, strip_new)
    return(o)
  }
  if (is.list(o)) o[] <- lapply(o, strip_new)
  o
}

read_status <- function(path) {
  if (!file.exists(path)) return(c(status = "MISSING", elapsed = NA))
  lines <- strsplit(readLines(path), "\t")
  setNames(vapply(lines, `[[`, "", 2L), vapply(lines, `[[`, "", 1L))
}

compare_one <- function(stem) {
  dirs <- file.path(S, "out", c("base", "cur"), stem)
  st <- lapply(file.path(dirs, "status.tsv"), read_status)
  env <- lapply(file.path(dirs, "env.rds"), \(p) if (file.exists(p)) readRDS(p) else list())
  common <- intersect(names(env[[1]]), names(env[[2]]))
  obj_check <- vapply(common, \(nm) {
    a <- strip_new(env[[1]][[nm]]); b <- strip_new(env[[2]][[nm]])
    if (identical(a, b)) return("identical")
    ae <- all.equal(a, b, tolerance = 1e-8, check.environment = FALSE)
    if (isTRUE(ae)) "equal_1e-8" else paste(head(ae, 3), collapse = " | ")
  }, character(1))
  md <- lapply(file.path(dirs, paste0(stem, ".md")),
               \(p) if (file.exists(p)) readLines(p, warn = FALSE) else character())
  only_cur <- setdiff(md[[2]], md[[1]]); only_base <- setdiff(md[[1]], md[[2]])
  figs <- lapply(file.path(dirs, "figure"), \(d)
    tools::md5sum(sort(list.files(d, full.names = TRUE))) |> setNames(nm = NULL))
  fig_names <- lapply(file.path(dirs, "figure"), \(d) sort(list.files(d)))
  list(stem = stem, status_base = st[[1]][["status"]], status_cur = st[[2]][["status"]],
       elapsed_base = st[[1]][["elapsed"]], elapsed_cur = st[[2]][["elapsed"]],
       objects = obj_check, only_cur = only_cur, only_base = only_base,
       n_fig = length(fig_names[[2]]),
       figs_differ = if (identical(fig_names[[1]], fig_names[[2]]))
         fig_names[[2]][figs[[1]] != figs[[2]]] else "FIGURE SET DIFFERS")
}

stems <- list.dirs(file.path(S, "out", "cur"), recursive = FALSE, full.names = FALSE)
results <- lapply(stems, compare_one)
saveRDS(results, file.path(S, "comparison.rds"))
summary_table <- do.call(rbind, lapply(results, \(r) data.frame(
  vignette = r$stem, base = r$status_base, cur = r$status_cur,
  sec_base = r$elapsed_base, sec_cur = r$elapsed_cur,
  objects = length(r$objects), objects_changed = sum(!r$objects %in% c("identical", "equal_1e-8")),
  md_lines_new = length(r$only_cur), md_lines_gone = length(r$only_base),
  figures = r$n_fig, figures_changed = length(r$figs_differ))))
print(summary_table, row.names = FALSE)
invisible(lapply(results, \(r) {
  changed <- r$objects[!r$objects %in% c("identical", "equal_1e-8")]
  if (length(changed) == 0L && length(r$only_cur) == 0L && length(r$only_base) == 0L &&
      length(r$figs_differ) == 0L) return(NULL)
  cat("\n=====", r$stem, "=====\n")
  if (length(changed)) { cat("-- objects changed:\n"); print(changed) }
  if (length(r$only_base)) { cat("-- md lines only in BASE:\n"); writeLines(head(r$only_base, 40)) }
  if (length(r$only_cur)) { cat("-- md lines only in CURRENT:\n"); writeLines(head(r$only_cur, 40)) }
  if (length(r$figs_differ)) { cat("-- figures differing:\n"); print(r$figs_differ) }
}))
