# heavy_02 under search = "emEM": does the opt-in search reach the maximum
# that Latent GOLD and Mplus report (-3098.4868) where the default search
# stops at -3106.9860 (finding 1)? Same model and controls as
# heavy-02-lower-maximum.R, plus a typical budget of 10 starts from 5 seeds.
# From the project root (after make-kit.R):
#   Rscript equivalence/features/direct-effects/repro/heavy-02-emem.R
pkgload::load_all(".", quiet = TRUE)
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
sys.source(file.path(root, "reference.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
reference <- -3098.4868

settings <- rbind(
  data.frame(search = "kmeans", n_starts = 200L, seed = 1:3),
  data.frame(search = "emEM", n_starts = 200L, seed = 1:3),
  data.frame(search = "kmeans", n_starts = 10L, seed = 1:5),
  data.frame(search = "emEM", n_starts = 10L, seed = 1:5))

rows <- lapply(seq_len(nrow(settings)), \(i) {
  setting <- settings[i, ]
  warned <- character()
  seconds <- system.time(fit <- withCallingHandlers(
    lpa(case$data, case$spec$vars, 2, model = "VVI", n_starts = setting$n_starts,
        seed = setting$seed, tol = 1e-10, max_iter = 5000,
        direct_effects = "y3 ~ z", search = setting$search),
    warning = \(w) {
      warned <<- c(warned, class(w)[[1L]])
      invokeRestart("muffleWarning")
    }))[["elapsed"]]
  starts <- get_results(fit, "starts")
  finalists <- if (is.null(starts$stage)) starts else subset(starts, stage == "final")
  reaching <- subset(finalists, abs(log_likelihood - fit$log_likelihood) < 1e-4)
  diagnostics <- get_results(fit, "gaussian_diagnostics")
  data.frame(setting, log_likelihood = fit$log_likelihood,
             gap_to_reference = fit$log_likelihood - reference,
             finalists = nrow(finalists), finalists_reaching = nrow(reaching),
             winning_kinds = paste(sort(unique(reaching$kind)), collapse = "+"),
             min_standardized_eigenvalue = min(diagnostics$min_standardized_eigenvalue),
             boundary = fit$boundary,
             warnings = paste(unique(warned), collapse = "; "),
             seconds = round(seconds, 1))
})
result <- do.call(rbind, rows)
print(result, row.names = FALSE)
saveRDS(result, file.path(root, "repro", "heavy-02-emem.rds"))
