# Record the golden baseline.
#
#   Rscript equivalence/engine-golden/record.R            # every case
#   Rscript equivalence/engine-golden/record.R '^lta_'    # cases matching a regex
#
# Writes golden/<case>.rds (the fingerprint), and merges the recorded cases
# into timing.csv (case, seconds), unexpected-errors.csv (any error that is
# not one of the package's classed refusals) and catalogue.csv (what each
# case covers). Set GOLDEN_OUT to write somewhere else (a scratch run).

local({
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                        value = TRUE)
  dir <- if (length(file_argument) == 0L) {
    normalizePath(file.path("equivalence", "engine-golden"))
  } else {
    dirname(normalizePath(sub("^--file=", "", file_argument[1L])))
  }
  source(file.path(dir, "setup.R"))
  golden_setup(dir)
  assign(".golden_dir", dir, envir = globalenv())
})

pattern <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(pattern)) pattern <- NULL
out_dir <- golden_output_dir(.golden_dir)
rds_dir <- file.path(out_dir, "golden")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

cases <- golden_catalogue(pattern)
if (length(cases) == 0L) stop("No case matches the pattern.", call. = FALSE)
cat(sprintf("Recording %d cases into %s\n", length(cases), rds_dir))

run_started <- proc.time()[["elapsed"]]
results <- lapply(cases, function(case) {
  outcome <- tryCatch(golden_run_case(case), error = function(e) e)
  if (inherits(outcome, "error")) {
    cat(sprintf("%-45s ERROR (harness): %s\n", case$name,
                conditionMessage(outcome)))
    return(list(case = case, seconds = NA_real_, status = "ERROR",
                unexpected = data.frame(
                  path = case$name, error_class = "harness",
                  message = conditionMessage(outcome),
                  stringsAsFactors = FALSE)))
  }
  saveRDS(outcome$fingerprint, file.path(rds_dir, paste0(case$name, ".rds")),
          version = 3L)
  unexpected <- golden_unexpected_errors(outcome$fingerprint)
  status <- if (!is.null(outcome$fingerprint$refusal)) "refused" else "fitted"
  cat(sprintf("%-45s %-8s %7.2f s  unexpected errors: %d\n", case$name, status,
              outcome$seconds, nrow(unexpected)))
  list(case = case, seconds = outcome$seconds, status = status,
       unexpected = unexpected)
})
total <- proc.time()[["elapsed"]] - run_started

timing <- data.frame(
  case = vapply(results, function(r) r$case$name, character(1L)),
  seconds = vapply(results, function(r) r$seconds, numeric(1L)),
  stringsAsFactors = FALSE)
golden_merge_csv(file.path(out_dir, "timing.csv"), timing)

catalogue_rows <- data.frame(
  case = timing$case,
  engine = vapply(results, function(r) r$case$engine, character(1L)),
  status = vapply(results, function(r) r$status, character(1L)),
  services = vapply(results, function(r) paste(r$case$services, collapse = " "),
                    character(1L)),
  covers = vapply(results, function(r) r$case$covers, character(1L)),
  stringsAsFactors = FALSE)
golden_merge_csv(file.path(out_dir, "catalogue.csv"), catalogue_rows)

unexpected <- do.call(rbind, lapply(results, function(r) {
  if (nrow(r$unexpected) == 0L) return(NULL)
  cbind(case = r$case$name, r$unexpected, stringsAsFactors = FALSE)
}))
if (is.null(unexpected)) {
  unexpected <- data.frame(case = character(), path = character(),
                           error_class = character(), message = character(),
                           stringsAsFactors = FALSE)
}
# Cases recorded now replace their old rows, including a row that is now gone.
unexpected_path <- file.path(out_dir, "unexpected-errors.csv")
old_unexpected <- if (file.exists(unexpected_path)) {
  utils::read.csv(unexpected_path, stringsAsFactors = FALSE)
} else {
  unexpected[0L, , drop = FALSE]
}
old_unexpected <- old_unexpected[!old_unexpected$case %in% timing$case, ,
                                 drop = FALSE]
utils::write.csv(rbind(old_unexpected[names(unexpected)], unexpected),
                 unexpected_path, row.names = FALSE)

writeLines(c(sprintf("latents %s", as.character(utils::packageVersion("latents"))),
             utils::capture.output(print(utils::sessionInfo()))),
           file.path(out_dir, "session-info.txt"))

cat(sprintf("\nRecorded %d cases in %.1f s (fits alone: %.1f s).\n",
            length(results), total, sum(timing$seconds, na.rm = TRUE)))
cat(sprintf("Refused fits (recorded refusals): %d\n",
            sum(catalogue_rows$status == "refused")))
if (nrow(unexpected) > 0L) {
  cat(sprintf("\nUnexpected (non-refusal) errors in probes: %d\n",
              nrow(unexpected)))
  print(unexpected, row.names = FALSE, right = FALSE)
}
if (any(catalogue_rows$status == "ERROR")) quit(status = 1L)
