library(latents)
stopifnot(as.character(packageVersion("latents")) == "0.96.1")
cat("Installed version:", as.character(packageVersion("latents")), "\n")
cat("Installed library:", find.package("latents"), "\n")
source("validation/fixes-2026-10-20/original-public-proofs.R")

# Confirm the built archive preserves the typed worker failure and closes sockets.
connections_before <- nrow(showConnections(all = TRUE))
warning_messages <- character()
withCallingHandlers({
  worker_failure <- tryCatch(latents:::.latents_map_candidates(1:3,
    function(i) {
      if (i == 2L) tools::pskill(Sys.getpid(), tools::SIGTERM)
      i
    }, cores = 2L, backend = "socket"), error = identity)
  stopifnot(inherits(worker_failure, "latents_parallel_failure"))
  invisible(gc())
  stopifnot(nrow(showConnections(all = TRUE)) == connections_before)
}, warning = function(condition) {
  warning_messages <<- c(warning_messages, conditionMessage(condition))
  invokeRestart("muffleWarning")
})
stopifnot(length(warning_messages) == 0L)
cat("\nInstalled archive: typed worker failure, zero cleanup warnings, no leaked sockets.\n")
