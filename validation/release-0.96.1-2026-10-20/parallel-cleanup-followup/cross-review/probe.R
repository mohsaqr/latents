pkgload::load_all(".", quiet = TRUE)
baseline <- nrow(showConnections(all = TRUE))
capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }), error = identity)
  withCallingHandlers(invisible(gc()), warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  })
  connections <- nrow(showConnections(all = TRUE))
  print(list(value = value, warnings = warnings, connections = connections))
  stopifnot(!length(warnings), connections == baseline)
  value
}
normal <- capture(.latents_map_candidates(1:3, identity, 2L, "socket"))
stopifnot(identical(normal, list(1L, 2L, 3L)))
remote <- capture(.latents_map_candidates(1:3, function(i) {
  if (i == 2L) stop("independent original worker error")
  i
}, 2L, "socket"))
stopifnot(inherits(remote, "error"),
          grepl("independent original worker error", conditionMessage(remote)))
dead <- capture(.latents_map_candidates(1:3, function(i) {
  if (i == 2L) tools::pskill(Sys.getpid(), tools::SIGTERM)
  i
}, 2L, "socket"))
stopifnot(inherits(dead, "latents_parallel_failure"))
exercise_mocked_cleanup <- function() {
  owned <- lapply(1:3, function(i) rawConnection(raw(), "w"))
  unrelated <- rawConnection(raw(), "w")
  on.exit(tryCatch(close(unrelated), error = function(condition) NULL), add = TRUE)
  cluster <- structure(lapply(owned, function(connection) list(con = connection)),
                        class = c("SOCKcluster", "cluster"))
  attempts <- 0L
  testthat::local_mocked_bindings(stopCluster = function(cluster) {
    attempts <<- attempts + 1L
    stop("forced independent cleanup failure")
  }, .package = "parallel")
  .latents_stop_cluster(cluster)
  stopifnot(attempts == 3L, isOpen(unrelated),
    all(vapply(owned, function(connection) {
      tryCatch({ isOpen(connection); FALSE }, error = function(condition) TRUE)
    }, logical(1))))
  .latents_stop_cluster(cluster)
  stopifnot(attempts == 6L, isOpen(unrelated))
}
capture(exercise_mocked_cleanup())
cat("Independent cleanup probes passed: normal/remote-error/dead-worker/mocked-shutdown,")
cat(" no warnings, no leaked connections, unrelated connection preserved.\n")
