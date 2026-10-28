Sys.setenv(VECLIB_MAXIMUM_THREADS = 1, OMP_NUM_THREADS = 1,
           OPENBLAS_NUM_THREADS = 1)
pkgload::load_all(quiet = TRUE)
quiet_fit <- function(...) {
  withCallingHandlers(lpa(...), warning = function(w) {
    cat("warning:", class(w)[1L], conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
}
items <- names(srl)
fits <- lapply(1:3, function(k) {
  fit <- quiet_fit(srl, items, 4, model = "EEE", n_starts = 1,
                   search = "emEM", n_screen = 10, seed = k)
  cat("\nSEED", k, "LL", fit$log_likelihood, "\n")
  print(get_results(fit, "starts"))
  fit
})
saveRDS(fits, "validation/review-fit-search-2026-10-22/root-search-probes.rds")
z <- readRDS("validation/review-fit-search-2026-10-22/zero-weight-starts.rds")
stopifnot(z$zero_at_zero$converged, z$zero_extreme$converged,
          z$zero_extreme100$converged,
          abs(z$zero_at_zero$log_likelihood - z$good$log_likelihood) < 1e-7,
          abs(z$emem$log_likelihood - z$good$log_likelihood) < 1e-7,
          z$good$log_likelihood - z$zero_extreme$log_likelihood > 100)
cat("\nZero-weight saved-object assertions passed\n")
d <- readRDS("validation/review-fit-search-2026-10-22/spike-diagnostic-probes.rds")
boundary <- getFromNamespace(".multilpa_covariance_boundary", "latents")
occasion <- d$lta[[2L]]$fit
stopifnot(occasion$converged, !occasion$boundary,
          all(vapply(occasion$measurement, boundary, logical(1),
                     min_variance = 1e-6)))
cat("LTA covariance-boundary mismatch independently verified\n")
