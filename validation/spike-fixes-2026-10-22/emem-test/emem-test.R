# Does a warm start from a less-constrained model's posteriors (JStats'
# strategy) reach the maxima the default search misses on the 9 VLMR cases?
# The package is unchanged: the warm start is passed through `start =`.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
here <- file.path("equivalence", "features", "vlmr")
cases <- readRDS(file.path(here, "kit", "cases.rds"))
comparisons <- utils::read.csv(file.path(here, "comparisons.csv"))
comparisons$k <- comparisons$null_k + 1L
missed <- data.frame(
  case = c("eei_k2_02", "eei_k3poor_02", "eei_k3poor_04", "eei_large_02", "eei_k1_01",
           "eei_k1_02", "eei_k1_05", "eee_k2_04", "eee_k2_04"),
  k = c(3L, 4L, 3L, 4L, 4L, 3L, 4L, 3L, 4L))
missed <- merge(missed, subset(comparisons, select = c(case, k, lg_ll1, mp_h1_ll)),
                by = c("case", "k"))
missed$reference <- pmax(missed$lg_ll1, missed$mp_h1_ll, na.rm = TRUE)
stopifnot("9 missed candidates with a finite reference" =
            nrow(missed) == 9L && all(is.finite(missed$reference)))

fit_ll <- function(case, k, ...) {
  spec <- case$spec
  warned <- character()
  seconds <- system.time(fit <- withCallingHandlers(
    lpa(case$data, spec$vars, k, model = spec$model, tol = 1e-12, max_iter = 5000, ...),
    warning = \(w) {
      warned <<- c(warned, class(w)[[1L]])
      invokeRestart("muffleWarning")
    }))[["elapsed"]]
  data.frame(log_likelihood = fit$log_likelihood, boundary = fit$boundary,
             warnings = paste(unique(warned), collapse = "; "), seconds = round(seconds, 1))
}

runs <- do.call(rbind, lapply(seq_len(nrow(missed)), \(i) {
  row <- missed[i, ]
  case <- cases[[row$case]]
  study <- if (case$spec$n >= 3000L) 30L else 60L
  settings <- rbind(
    data.frame(setting = "emEM, study starts", n_starts = study, n_screen = 5L * study, seed = 1:2),
    data.frame(setting = "emEM, Mplus default (20 4)", n_starts = 4L, n_screen = 20L, seed = 1:3))
  do.call(rbind, lapply(seq_len(nrow(settings)), \(j) {
    s <- settings[j, ]
    data.frame(case = row$case, k = row$k, model = case$spec$model, reference = row$reference, s,
               fit_ll(case, row$k, search = "emEM", n_starts = s$n_starts,
                      n_screen = s$n_screen, seed = s$seed))
  }))
}))
runs$gap <- runs$log_likelihood - runs$reference
runs$reached <- runs$gap >= -1e-3
options(width = 220)
print(subset(runs, select = -c(n_screen)), row.names = FALSE)
cat("\nReached, by setting:\n")
print(aggregate(reached ~ setting, data = runs, FUN = \(x) sprintf("%d of %d", sum(x), length(x))))
cat("\nCases reached by at least one seed, by setting:\n")
print(aggregate(reached ~ setting, data = aggregate(reached ~ setting + case + k, data = runs, FUN = any),
                FUN = \(x) sprintf("%d of %d", sum(x), length(x))))
cat("\nTotal seconds by setting:\n"); print(aggregate(seconds ~ setting, data = runs, FUN = sum))
saveRDS(runs, "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/emem-test.rds")
