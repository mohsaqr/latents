Sys.setenv(VECLIB_MAXIMUM_THREADS = 1, OMP_NUM_THREADS = 1,
           OPENBLAS_NUM_THREADS = 1)
pkgload::load_all(quiet = TRUE)
fit_search <- function() lpa(srl, names(srl), 4, model = "EEE",
                             n_starts = 1, search = "emEM", n_screen = 10,
                             seed = 1)
baseline <- fit_search()
trace(".latents_run_screened", tracer = quote(screen$protect <- integer()),
      print = FALSE, where = asNamespace("latents"))
ranked <- fit_search()
untrace(".latents_run_screened", where = asNamespace("latents"))
print(data.frame(protected = c(TRUE, FALSE),
                 log_likelihood = c(baseline$log_likelihood, ranked$log_likelihood),
                 converged = c(baseline$converged, ranked$converged)))
print(get_results(ranked, "starts"))
stopifnot(baseline$log_likelihood + 2.8 < ranked$log_likelihood,
          baseline$converged, ranked$converged,
          identical(baseline$starts$screened_log_likelihood,
                    ranked$starts$screened_log_likelihood))
saveRDS(list(baseline = baseline, unprotected = ranked),
        "validation/review-fit-search-2026-10-22/root-emem-protection.rds")
