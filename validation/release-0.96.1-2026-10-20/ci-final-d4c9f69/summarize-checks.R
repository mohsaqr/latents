output <- file.path("validation", "release-0.96.1-2026-10-20", "ci-final-d4c9f69")
expected_sha <- "d4c9f699525d9e07959463010e004b6362dcdaaf"
run <- jsonlite::fromJSON(file.path(output, "39456155958-latest.json"),
                           simplifyVector = FALSE)
stopifnot(identical(run$headSha, expected_sha), length(run$jobs) == 6L)
completed <- Filter(function(job) identical(job$status, "completed"), run$jobs)
rows <- lapply(completed, function(job) {
  path <- file.path(output, sprintf("%.0f-%s.log", job$databaseId, job$conclusion))
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  status_lines <- grep("Status: ", lines, value = TRUE)
  status <- if (length(status_lines)) sub(".*Status: ", "", tail(status_lines, 1L)) else NA_character_
  count_status <- function(kind) {
    if (is.na(status)) return(NA_integer_)
    if (!grepl(kind, status, fixed = TRUE)) return(0L)
    as.integer(sub(sprintf(".*?([0-9]+) %s.*", kind), "\\1", status))
  }
  totals <- grep("\\[ FAIL [0-9]+ \\| WARN [0-9]+ \\| SKIP [0-9]+ \\| PASS [0-9]+ \\]",
                 lines, value = TRUE)
  totals <- if (length(totals)) tail(totals, 1L) else NA_character_
  count_test <- function(kind) {
    if (is.na(totals)) return(NA_integer_)
    as.integer(sub(sprintf(".*%s ([0-9]+).*", kind), "\\1", totals))
  }
  notes <- grep("checking .*NOTE", lines, value = TRUE)
  note_indices <- grep("Z \\* checking .*NOTE", lines)
  boundaries <- grep("Z \\* (checking|DONE)", lines)
  note_details <- paste(vapply(note_indices, function(index) {
    end <- c(boundaries[boundaries > index], length(lines) + 1L)[[1L]] - 1L
    paste(trimws(sub(".*Z ", "", lines[seq.int(index, end)])), collapse = " | ")
  }, character(1)), collapse = "; ")
  versions <- grep("\\* using R (version|Under development)", lines, value = TRUE)
  actual_r <- if (length(versions)) sub(".*\\* using ", "", versions[[1L]]) else NA_character_
  skip_headers <- grep("Skipped tests ", lines)
  skip_reasons <- if (length(skip_headers)) {
    start <- tail(skip_headers, 1L) + 1L
    endings <- grep("\\[ FAIL ", lines)
    sections <- grep("Z ══", lines)
    end <- min(c(endings[endings >= start], sections[sections >= start])) - 1L
    block <- trimws(sub(".*Z ", "", lines[seq.int(start, end)]))
    paste(block[nzchar(block)], collapse = "; ")
  } else ""
  data.frame(platform = job$name, head_sha = expected_sha,
    actual_r = actual_r,
    job_id = as.character(job$databaseId), conclusion = job$conclusion,
    check_status = status, errors = count_status("ERROR"),
    warnings = count_status("WARNING"), notes = count_status("NOTE"),
    unit_pass = count_test("PASS"), unit_fail = count_test("FAIL"),
    unit_warn = count_test("WARN"), unit_skip = count_test("SKIP"),
    skip_reasons = skip_reasons,
    note_topics = paste(sub(".*Z ", "", notes), collapse = "; "),
    note_details = note_details,
    log = basename(path), url = job$url)
})
result <- do.call(rbind, rows)
if (is.null(result)) stop("No completed matrix logs available yet.")
utils::write.csv(result, file.path(output, "check-summary.csv"), row.names = FALSE)
platforms <- result[c("platform", "actual_r", "job_id", "errors", "warnings", "notes",
                       "unit_pass", "unit_skip", "skip_reasons", "url")]
names(platforms) <- c("config", "actual_r_version", "job_id", "errors", "warnings", "notes",
                       "pass", "skipped", "skip_reasons", "url")
utils::write.csv(platforms, file.path(output, "platform-summary.csv"), row.names = FALSE)
equivalence_path <- file.path(output, "39456160271-complete.log")
equivalence_totals <- if (file.exists(equivalence_path)) {
  lines <- readLines(equivalence_path, warn = FALSE)
  unique(sub(".*Z ", "", grep("[[] FAIL [0-9]+", lines, value = TRUE)))
} else "Complete equivalence log not available yet."
lint_path <- file.path(output, "39456165092-complete.log")
lint_summary <- if (file.exists(lint_path) &&
  any(grepl("No lints found", readLines(lint_path, warn = FALSE), fixed = TRUE))) {
  "Successful; No lints found."
} else "Successful lint log not available yet."
summary <- c(
  "Release 0.96.1 final CI platform evidence",
  paste("Exact commit:", expected_sha),
  sprintf("Completed matrix job logs: %d/6", nrow(result)),
  if (any(result$conclusion != "success")) "Candidate has a failed platform check and requires a superseding fix."
    else if (nrow(result) == 6L) "All six platform checks passed."
    else "Platform checks are still running; this is a partial summary.",
  "All values below are read from the actual per-platform job logs.",
  "",
  unlist(lapply(seq_len(nrow(result)), function(index) {
    row <- result[index, ]
    c(sprintf("%s — %s", row$platform, row$actual_r),
      sprintf("  Job %s: %s; check %s; %d errors, %d warnings, %d NOTEs.",
        row$job_id, row$conclusion, row$check_status, row$errors, row$warnings, row$notes),
      sprintf("  Unit expectations: %d passed, %d failed, %d warnings; %d tests skipped.",
        row$unit_pass, row$unit_fail, row$unit_warn, row$unit_skip),
      paste("  Skips:", row$skip_reasons),
      if (nzchar(row$note_topics)) paste("  NOTE topics:", row$note_topics),
      if (nzchar(row$note_details)) paste("  NOTE details:", row$note_details),
      paste("  URL:", row$url), "")
  })),
  "Equivalence run 39456160271 (inner suite then retained references/claims):",
  equivalence_totals,
  paste("Lint run 39456165092:", lint_summary),
  "Prior failed/cancelled candidates are preserved in separate directories.")
writeLines(summary, file.path(output, "SUMMARY.txt"))
print(result[c("platform", "conclusion", "check_status", "errors", "warnings", "notes",
                "unit_pass", "unit_fail", "unit_warn", "unit_skip")], row.names = FALSE)
