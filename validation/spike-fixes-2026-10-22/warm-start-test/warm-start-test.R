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

# One M-step of the target equal-covariance model from given posteriors.
start_from_posteriors <- function(x, posteriors, model) {
  weights <- colSums(posteriors)
  means <- crossprod(posteriors, x) / weights
  pooled <- Reduce(`+`, lapply(seq_len(ncol(posteriors)), \(j) {
    centred <- sweep(x, 2L, means[j, ], "-")
    crossprod(centred * posteriors[, j], centred)
  })) / nrow(x)
  k <- ncol(posteriors)
  start <- list(means = unname(means),
                variances = matrix(diag(pooled), k, ncol(x), byrow = TRUE),
                profile_probabilities = matrix(weights / sum(weights), 1L),
                group_probabilities = 1)
  if (identical(model, "EEE")) start$covariances <- array(pooled, c(ncol(x), ncol(x), k))
  start
}

quiet_lpa <- function(...) {
  warned <- character()
  fit <- withCallingHandlers(lpa(...), warning = \(w) {
    warned <<- c(warned, class(w)[[1L]])
    invokeRestart("muffleWarning")
  })
  fit$warned <- unique(warned)
  fit
}

rows <- lapply(seq_len(nrow(missed)), \(i) {
  row <- missed[i, ]
  case <- cases[[row$case]]
  spec <- case$spec
  x <- as.matrix(case$data[spec$vars])
  n_starts <- if (spec$n >= 3000L) 30L else 60L
  controls <- list(n_starts = n_starts, seed = 1L, tol = 1e-12, max_iter = 5000)
  default <- do.call(quiet_lpa, c(list(case$data, spec$vars, row$k, model = spec$model), controls))
  warm_models <- if (identical(spec$model, "EEE")) "VVV" else c("VVV", "EEE")
  warm <- lapply(warm_models, \(warm_model) {
    base <- do.call(quiet_lpa, c(list(case$data, spec$vars, row$k, model = warm_model), controls))
    start <- start_from_posteriors(x, base$subject_posteriors, spec$model)
    refit <- quiet_lpa(case$data, spec$vars, row$k, model = spec$model, start = start,
                       n_starts = 1L, seed = 1L, tol = 1e-12, max_iter = 5000)
    data.frame(warm_model = warm_model, base_converged = base$converged,
               warm_ll = refit$log_likelihood, warm_converged = refit$converged)
  })
  warm <- do.call(rbind, warm)
  best_warm <- max(warm$warm_ll)
  data.frame(case = row$case, k = row$k, model = spec$model, n = spec$n,
             reference = row$reference, default_ll = default$log_likelihood,
             warm_from = paste(sprintf("%s:%.3f", warm$warm_model, warm$warm_ll), collapse = " "),
             best_ll = max(default$log_likelihood, best_warm),
             default_gap = default$log_likelihood - row$reference,
             with_warm_gap = max(default$log_likelihood, best_warm) - row$reference)
})
result <- do.call(rbind, rows)
result$default_reaches <- result$default_gap >= -1e-3
result$with_warm_reaches <- result$with_warm_gap >= -1e-3
options(width = 220)
print(result, row.names = FALSE)
cat(sprintf("\nReach the Mplus/Latent GOLD maximum: default %d of 9, default + warm starts %d of 9\n",
            sum(result$default_reaches), sum(result$with_warm_reaches)))
saveRDS(result, "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/warm-start-test.rds")
