#' Combine worker-loss reports from the phases of one fitting call
#'
#' Each mapper retries independently. Screening and finalist pools share this
#' reporter so a public call emits one qualification, including on error.
#' @return A collector for zero-argument callbacks and an exit-time reporter.
#' @noRd
.latents_worker_reporter <- function() {
  reports <- list()
  list(collect = function(fun) {
    stopifnot("fun is a function" = is.function(fun))
    force(fun)
    withCallingHandlers(fun(), latents_worker_lost = function(condition) {
      reports[[length(reports) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    })
  }, signal = function() {
    if (length(reports)) {
      .latents_warn_worker_lost(
        sum(vapply(reports, `[[`, integer(1), "lost")),
        sum(vapply(reports, `[[`, integer(1), "recovered")),
        max(vapply(reports, `[[`, integer(1), "retry_cores")))
    }
    invisible(NULL)
  })
}
