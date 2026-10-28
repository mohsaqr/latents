output <- file.path("validation", "release-0.96.1-2026-10-20", "ci-final-cf1127f")
expected_sha <- "cf1127fa66283826ed421d028999b5ac23e1901d"
run_ids <- c("39449966461", "39449971798", "39449977056")
previous <- NULL

read_run <- function(id) {
  path <- file.path(output, paste0(id, "-latest.json"))
  status <- system2("gh", c("run", "view", id, "--json",
    "databaseId,name,headSha,status,conclusion,createdAt,updatedAt,jobs,url"),
    stdout = path, stderr = file.path(output, paste0(id, "-read.stderr")))
  stopifnot(status == 0L)
  result <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  stopifnot(identical(result$headSha, expected_sha))
  observed_at <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  file.copy(path, file.path(output, paste0(id, "-", observed_at, ".json")))
  result
}

save_completed_logs <- function(run) {
  invisible(lapply(run$jobs, function(job) {
    if (!identical(job$status, "completed")) return(NULL)
    path <- file.path(output, sprintf("%.0f-%s.log", job$databaseId, job$conclusion))
    if (file.exists(path) && file.info(path)$size > 0) return(NULL)
    status <- system2("gh", c("api", sprintf(
      "repos/mohsaqr/latents/actions/jobs/%.0f/logs", job$databaseId)),
      stdout = path, stderr = paste0(path, ".stderr"))
    if (status != 0L) unlink(path)
    cat(sprintf("LOG: %s %s %s (download status %d)\n", run$name,
                 job$name, job$conclusion, status))
  }))
  if (identical(run$status, "completed")) {
    path <- file.path(output, sprintf("%.0f-complete.log", run$databaseId))
    if (!file.exists(path) || file.info(path)$size == 0) {
      status <- system2("gh", c("run", "view", as.character(run$databaseId), "--log"),
                         stdout = path, stderr = paste0(path, ".stderr"))
      if (status != 0L) unlink(path)
    }
  }
}

monitor_once <- function() {
  runs <- lapply(run_ids, read_run)
  jobs <- do.call(rbind, lapply(runs, function(run) {
    do.call(rbind, lapply(run$jobs, function(job) {
      active <- Filter(function(step) identical(step$status, "in_progress"), job$steps)
      step <- if (length(active)) active[[1L]]$name else ""
      data.frame(run_id = as.character(run$databaseId), run = run$name,
        head_sha = run$headSha, job_id = as.character(job$databaseId), job = job$name,
        status = job$status, conclusion = job$conclusion, step = step, url = job$url)
    }))
  }))
  utils::write.csv(jobs, file.path(output, "jobs-latest.csv"), row.names = FALSE)
  current <- paste(jobs$job, jobs$status, jobs$conclusion, jobs$step)
  if (!identical(previous, current)) {
    cat("UTC:", format(Sys.time(), tz = "UTC"), "head SHA verified\n")
    cat(paste(current, collapse = "\n"), "\n")
    previous <<- current
  }
  invisible(lapply(runs, save_completed_logs))
  flush.console()
  all(vapply(runs, function(run) identical(run$status, "completed"), logical(1)))
}

repeat {
  outcome <- tryCatch(monitor_once(), error = function(error) {
    cat("MONITOR READ ERROR:", conditionMessage(error), "\n")
    flush.console()
    FALSE
  })
  if (isTRUE(outcome)) break
  Sys.sleep(45)
}
cat("ALL FINAL RUNS TERMINAL\n")
