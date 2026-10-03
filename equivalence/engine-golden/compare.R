# Compare the current code against the golden baseline.
#
#   Rscript equivalence/engine-golden/compare.R            # every case
#   Rscript equivalence/engine-golden/compare.R '^lta_'    # cases matching a regex
#
# Re-fits each case, fingerprints it, and compares with golden/<case>.rds:
# text, logicals, integers and names identical; doubles within 1e-10 relative
# for log likelihoods and 1e-8 for everything else (an absolute difference of
# at most 1e-10 always passes). Prints one row per case and writes
# last-compare.csv. Exits non-zero when any case is DIFFERENT or ERROR. A fit
# more than 1.2 times slower than its baseline, and slower by more than 0.05 s
# (below that, a ratio of milliseconds is clock noise), is flagged in `slow`,
# which is a warning, not a failure: timing is noisy.

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
timing_path <- file.path(out_dir, "timing.csv")
baseline_timing <- if (file.exists(timing_path)) {
  utils::read.csv(timing_path, stringsAsFactors = FALSE)
} else {
  data.frame(case = character(), seconds = numeric())
}

cases <- golden_catalogue(pattern)
if (length(cases) == 0L) stop("No case matches the pattern.", call. = FALSE)
status_labels <- c("identical", "within_tolerance", "DIFFERENT")
slow_ratio <- 1.2
slow_floor <- 0.05

rows <- lapply(cases, function(case) {
  baseline_seconds <- baseline_timing$seconds[match(case$name,
                                                    baseline_timing$case)]
  row <- function(status, max_rel = NA_real_, path = NA_character_,
                  detail = NA_character_, seconds = NA_real_) {
    ratio <- seconds / baseline_seconds
    data.frame(case = case$name, status = status, max_rel_diff = max_rel,
               first_diff_path = path, seconds_baseline = baseline_seconds,
               seconds_now = seconds, slowdown = ratio,
               slow = isTRUE(ratio > slow_ratio &&
                               seconds - baseline_seconds > slow_floor),
               detail = detail,
               stringsAsFactors = FALSE)
  }
  rds <- file.path(rds_dir, paste0(case$name, ".rds"))
  if (!file.exists(rds)) {
    result <- row("ERROR", detail = "no baseline recorded")
  } else {
    baseline <- readRDS(rds)
    run <- tryCatch(golden_run_case(case), error = function(e) e)
    result <- if (inherits(run, "error")) {
      row("ERROR", detail = paste("harness:", conditionMessage(run)))
    } else {
      outcome <- golden_compare(baseline, run$fingerprint, case$name)
      row(status_labels[outcome$status + 1L], outcome$max_rel, outcome$path,
          outcome$detail, run$seconds)
    }
  }
  cat(sprintf("%-45s %-16s max_rel %-9s %6.2fs (baseline %6.2fs)%s\n",
              result$case, result$status,
              formatC(result$max_rel_diff, format = "g", digits = 3),
              result$seconds_now, result$seconds_baseline,
              if (result$slow) "  SLOW" else ""))
  result
})
report <- do.call(rbind, rows)
utils::write.csv(report, file.path(out_dir, "last-compare.csv"),
                 row.names = FALSE)

cat("\n")
print(report[c("case", "status", "max_rel_diff", "first_diff_path",
               "seconds_now", "seconds_baseline", "slowdown", "slow")],
      row.names = FALSE, digits = 3)
counts <- table(factor(report$status, c(status_labels, "ERROR")))
cat(sprintf("\n%s\n", paste(sprintf("%s: %d", names(counts), counts),
                            collapse = "; ")))
cat(sprintf("Slower than %.1fx baseline (warning only): %d\n", slow_ratio,
            sum(report$slow)))
failed <- report$status %in% c("DIFFERENT", "ERROR")
if (any(failed)) {
  cat("\nFirst differences:\n")
  print(report[failed, c("case", "first_diff_path", "detail")],
        row.names = FALSE, right = FALSE)
  quit(status = 1L)
}
