# The 9 VLMR equal-covariance misses at the DEFAULT budget (n_starts = 10,
# n_screen = 50): what a user gets from the new default search (search = NULL,
# emEM for EEI/EEE) against the old default (kmeans), 3 seeds each.
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
runs <- do.call(rbind, lapply(seq_len(nrow(missed)), \(i) {
  row <- missed[i, ]
  case <- cases[[row$case]]
  do.call(rbind, lapply(c("default", "kmeans"), \(setting) {
    do.call(rbind, lapply(1:3, \(seed) {
      seconds <- system.time(fit <- withCallingHandlers(
        if (setting == "default") lpa(case$data, case$spec$vars, row$k, model = case$spec$model,
                                      tol = 1e-12, max_iter = 5000, seed = seed)
        else lpa(case$data, case$spec$vars, row$k, model = case$spec$model,
                 tol = 1e-12, max_iter = 5000, seed = seed, search = "kmeans"),
        warning = \(w) invokeRestart("muffleWarning")))[["elapsed"]]
      data.frame(case = row$case, k = row$k, model = case$spec$model, n = case$spec$n,
                 setting, seed, search = fit$search, gap = fit$log_likelihood - row$reference,
                 seconds = round(seconds, 1))
    }))
  }))
}))
runs$reached <- runs$gap >= -1e-3
options(width = 200)
print(runs, row.names = FALSE)
cat("\nRuns reaching the reference, by setting:\n")
print(aggregate(reached ~ setting, data = runs, FUN = \(x) sprintf("%d of %d", sum(x), length(x))))
cat("\nCases reached by at least one seed:\n")
print(aggregate(reached ~ setting, data = aggregate(reached ~ setting + case + k, data = runs, FUN = any),
                FUN = \(x) sprintf("%d of %d", sum(x), length(x))))
cat("\nSeconds, total by setting:\n"); print(aggregate(seconds ~ setting, data = runs, FUN = sum))
saveRDS(runs, file.path("validation", "spike-fixes-2026-10-22", "emem-test", "emem-default-budget.rds"))
