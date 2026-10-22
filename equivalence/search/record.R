# Rscript record.R <out_dir> [case ...]
# Runs the identity panel under the latents found on .libPaths() and saves,
# per case, the result, the RNG state on exit and the warnings. Run under
# R_LIBS pointing at the library to record.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(latents))
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "cases.R"))
# PANEL_SEARCH=kmeans runs every case with `search = "kmeans"`, for a build
# whose default is emEM.
if (nzchar(Sys.getenv("PANEL_SEARCH"))) options(latents.panel_search = Sys.getenv("PANEL_SEARCH"))
chosen <- if (length(args) > 1L) args[-1] else names(cases)

record <- function(name) {
  # Every case starts from the same RNG state, so `seed = NULL` cases are
  # compared on the draws the search itself makes.
  set.seed(20261021L)
  warnings <- character()
  result <- tryCatch(withCallingHandlers(cases[[name]](),
    warning = function(w) {
      warnings <<- c(warnings, paste(class(w)[1], conditionMessage(w)))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) structure(list(message = conditionMessage(e),
                                       class = class(e)), class = "panel_error"))
  entry <- list(result = result, random_seed = get(".Random.seed", globalenv()),
                warnings = warnings, version = as.character(packageVersion("latents")),
                r = R.version.string, rng = RNGkind())
  saveRDS(entry, file.path(out, paste0(name, ".rds")))
  cat(sprintf("%-26s %s\n", name,
              if (inherits(result, "panel_error")) paste("ERROR:", result$message) else "ok"))
}
invisible(lapply(chosen, record))
cat(capture.output(sessionInfo()), file = file.path(out, "session.txt"), sep = "\n")
