# Rscript compare.R <dir_a> <dir_b>
# Compares two recordings of the identity panel case by case: the result
# (every element, recursively, with identical()), the RNG state on exit and
# the warnings. Prints the paths that differ; exits non-zero when any does.
args <- commandArgs(trailingOnly = TRUE)
dir_a <- args[1]
dir_b <- args[2]
suppressPackageStartupMessages(library(latents))

# Paths whose values differ between two runs of the same code (wall-clock
# timings); found by recording the same library twice.
volatile <- c("seconds", "elapsed", "time", "timing", "started", "finished")

differences <- function(a, b, path = "") {
  if (identical(a, b, ignore.environment = TRUE, ignore.bytecode = TRUE,
                ignore.srcref = TRUE)) return(character())
  leaf <- sub(".*\\$", "", path)
  if (leaf %in% volatile) return(character())
  if (is.function(a) && is.function(b)) return(character())
  # A field the recorded release did not have is allowed (new controls);
  # every field it had must be identical, and in the same order.
  if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.null(names(a)) &&
      all(names(a) %in% names(b)) &&
      identical(names(a), intersect(names(b), names(a)))) {
    same_attributes <- identical(attributes(a)[setdiff(names(attributes(a)), "names")],
                                 attributes(b)[setdiff(names(attributes(b)), "names")])
    return(c(if (!same_attributes) paste0(path, " (attributes)"),
             unlist(lapply(names(a), function(key) {
               differences(a[[key]], b[[key]], paste0(path, "$", key))
             }))))
  }
  if (is.list(a) && is.list(b) && !is.data.frame(a) && is.null(names(a)) &&
      length(a) == length(b)) {
    return(unlist(lapply(seq_along(a), function(i) {
      differences(a[[i]], b[[i]], paste0(path, "[[", i, "]]"))
    })))
  }
  # The panel calls verbs through do.call() with the function itself, so a
  # recorded call carries the function, whose formals gain the new arguments;
  # the call's arguments are what must match.
  if (is.call(a) && is.call(b)) {
    head_a <- a[[1L]]
    head_b <- b[[1L]]
    heads_match <- (is.function(head_a) && is.function(head_b)) ||
      identical(head_a, head_b)
    # A build whose default is emEM is recorded with `search = "kmeans"`
    # passed explicitly, which is the call the release made by default.
    arguments <- function(call) {
      listed <- as.list(call)[-1L]
      if (identical(listed$search, "kmeans")) listed$search <- NULL
      listed
    }
    if (heads_match && identical(arguments(a), arguments(b))) return(character())
    return(paste0(path, ": ", deparse1(as.list(a)[-1L]), " vs ",
                  deparse1(as.list(b)[-1L])))
  }
  if (is.data.frame(a) && is.data.frame(b) && identical(names(a), names(b)) &&
      nrow(a) == nrow(b)) {
    return(unlist(lapply(names(a), function(column) {
      differences(a[[column]], b[[column]], paste0(path, "$", column))
    })))
  }
  detail <- if (is.numeric(a) && is.numeric(b) && length(a) == length(b)) {
    sprintf(" (max abs diff %.3g)", max(abs(a - b), na.rm = TRUE))
  } else ""
  paste0(path, detail)
}

files <- list.files(dir_a, "\\.rds$")
report <- lapply(files, function(file) {
  a <- readRDS(file.path(dir_a, file))
  b_path <- file.path(dir_b, file)
  if (!file.exists(b_path)) return(c(case = sub("\\.rds$", "", file), what = "missing in b"))
  b <- readRDS(b_path)
  what <- c(differences(a$result, b$result, "result"),
            if (!identical(a$random_seed, b$random_seed)) "random_seed",
            if (!identical(a$warnings, b$warnings)) "warnings")
  if (length(what)) cbind(case = sub("\\.rds$", "", file), what = what)
})
report <- do.call(rbind, report)
if (is.null(report)) {
  cat(sprintf("identical: %d cases\n", length(files)))
} else {
  print(as.data.frame(report), row.names = FALSE, right = FALSE)
  quit(status = 1L)
}
