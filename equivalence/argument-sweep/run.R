# Argument sweep: every fitting argument, across its options and several
# datasets, against a reference that does not use the package's likelihood
# code (a closed form, another package, or the mixture likelihood rebuilt from
# the fit's public result tables), plus invariances that any correct
# implementation must satisfy.
#
# From the project root:
#   Rscript equivalence/argument-sweep/run.R            # every module
#   Rscript equivalence/argument-sweep/run.R gaussian   # modules whose name matches
#
# Writes tmp/argument-sweep/results.csv (one row per check) and
# tmp/argument-sweep/SUMMARY.md, and exits non-zero when any check disagrees
# or any case stops with an error.
requested <- commandArgs(trailingOnly = TRUE)
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "argument-sweep")
source(file.path(directory, "helpers.R"))
source(file.path(directory, "datasets.R"))
module_files <- list.files(file.path(directory, "modules"), pattern = "[.]R$", full.names = TRUE)
if (length(requested)) {
  module_files <- module_files[grepl(paste(requested, collapse = "|"), basename(module_files))]
}
stopifnot("no module matches the request" = length(module_files) > 0L)

# Each module file defines one function, named after the file, returning cases.
modules <- lapply(module_files, function(path) {
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  name <- sub("[.]R$", "", basename(path))
  list(name = name, cases = get(sprintf("sweep_module_%s", gsub("-", "_", name)), environment)())
})
jobs <- unlist(lapply(modules, function(module) {
  lapply(module$cases, function(definition) list(module = module$name, definition = definition))
}), recursive = FALSE)
cores <- as.integer(Sys.getenv("SWEEP_CORES", "6"))
message(sprintf("%d cases in %d modules on %d cores", length(jobs), length(modules), cores))
RNGkind("L'Ecuyer-CMRG")
set.seed(20261111)
outcomes <- parallel::mclapply(jobs, function(job) {
  sweep_run_case(job$module, job$definition)
}, mc.cores = cores, mc.set.seed = TRUE, mc.preschedule = FALSE)
crashed <- !vapply(outcomes, is.data.frame, logical(1))
outcomes[crashed] <- lapply(which(crashed), function(i) {
  data.frame(module = jobs[[i]]$module, case = jobs[[i]]$definition$case,
             quantity = "(case ran)", difference = Inf, tolerance = 0, n_values = 0L,
             status = "error", note = paste("worker failed:", paste(format(outcomes[[i]]), collapse = " ")),
             conditions = NA_character_, seconds = NA_real_, stringsAsFactors = FALSE)
})
results <- do.call(rbind, outcomes)

# A disagreement the register explains is a known difference; an entry that
# explains nothing any more is stale and fails the run, so the register cannot
# outlive what it describes. Only a complete run can call an entry stale.
known <- utils::read.csv(file.path(directory, "known-differences.csv"), stringsAsFactors = FALSE)
key <- function(table) paste(table$module, table$case, table$quantity, sep = "\r")
explained <- results$status == "disagrees" & key(results) %in% key(known)
results$status[explained] <- "known difference"
results$note[explained] <- known$reason[match(key(results)[explained], key(known))]
stale <- if (length(requested)) known[0L, ] else known[!key(known) %in% key(results)[explained], ]

output <- file.path("tmp", "argument-sweep")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
suffix <- if (length(requested)) paste0("-", paste(gsub("[^a-z0-9]+", "", requested), collapse = "-")) else ""
utils::write.csv(results, file.path(output, sprintf("results%s.csv", suffix)), row.names = FALSE)
tally <- as.data.frame(table(module = results$module, status = results$status),
                       stringsAsFactors = FALSE)
tally <- tally[tally$Freq > 0L, ]
problems <- results[results$status %in% c("disagrees", "error"), ]
if (nrow(stale)) {
  problems <- rbind(problems, data.frame(
    module = stale$module, case = stale$case, quantity = stale$quantity, difference = NA_real_,
    tolerance = NA_real_, n_values = 0L, status = "stale known difference",
    note = "the register explains a disagreement that no longer occurs; remove the entry",
    conditions = NA_character_, seconds = NA_real_, stringsAsFactors = FALSE))
}
summary_lines <- c(
  sprintf("# Argument sweep, latents %s", utils::packageVersion("latents")),
  "", sprintf("%d cases, %d checks.", length(jobs), nrow(results)), "",
  "| module | status | checks |", "|---|---|---|",
  sprintf("| %s | %s | %d |", tally$module, tally$status, tally$Freq), "",
  sprintf("## Disagreements and errors (%d)", nrow(problems)), "",
  sprintf("- `%s` / `%s` / %s: difference %.3g (tolerance %.3g). %s", problems$module,
          problems$case, problems$quantity, problems$difference, problems$tolerance,
          ifelse(is.na(problems$note), "", problems$note)))
writeLines(summary_lines, file.path(output, sprintf("SUMMARY%s.md", suffix)))
writeLines(utils::capture.output(utils::sessionInfo()), file.path(output, "SESSION.txt"))
cat(summary_lines, sep = "\n")
if (nrow(problems)) quit(status = 1L)
