output <- file.path("validation", "release-0.96.1-2026-10-20", "ci-final-cf1127f")
expected_sha <- "cf1127fa66283826ed421d028999b5ac23e1901d"
run_ids <- c("39449966461", "39449971798", "39449977056")
observed_at <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
runs <- lapply(run_ids, function(id) {
  path <- file.path(output, paste0(id, "-latest.json"))
  run <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  stopifnot(identical(run$headSha, expected_sha))
  file.copy(path, file.path(output, paste0(id, "-", observed_at, ".json")))
  cat(sprintf("%s: %s %s; exact commit verified\n", run$name,
               run$status, run$conclusion))
  invisible(lapply(run$jobs, function(job) {
    active <- Filter(function(step) identical(step$status, "in_progress"), job$steps)
    step <- if (length(active)) active[[1L]]$name else ""
    cat(sprintf("  %s: %s %s %s\n", job$name, job$status, job$conclusion, step))
  }))
  run
})
jobs <- do.call(rbind, lapply(runs, function(run) {
  do.call(rbind, lapply(run$jobs, function(job) {
    data.frame(run_id = as.character(run$databaseId), run = run$name,
                head_sha = run$headSha, job_id = as.character(job$databaseId),
                job = job$name, status = job$status, conclusion = job$conclusion,
                url = job$url)
  }))
}))
utils::write.csv(jobs, file.path(output, "jobs-latest.csv"), row.names = FALSE)
cat("Observed UTC:", observed_at, "\n")
