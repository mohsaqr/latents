# Type-I error of bootstrap_lrt() under the emEM search, observed fits and
# bootstrap refits searched alike (the refits inherit the observed fits'
# controls). Same null data as the 2026-10-21 BLRT study (VVI, 2 against 3
# profiles, n = 300; and EEE, 1 against 2, srl-like), smaller: the question is
# whether emEM with its admissibility rule keeps the test calibrated where
# Mplus with heavy starts reached 0.33 on spikes.
#
# Rscript equivalence/search/blrt-emem.R <out_dir> <condition> <first> <last> <iter> <cores>
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
condition <- args[2]
replicates <- seq(as.integer(args[3]), as.integer(args[4]))
iter <- as.integer(args[5])
cores <- as.integer(args[6])
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({
  library(latents)
  library(parallel)
})
srl_mean <- colMeans(srl)
srl_cov <- stats::cov(srl)
conditions <- list(
  vvi_2v3 = list(base = 200000L, k0 = 2L, model = "VVI", make = function(seed) {
    d <- as.data.frame(simulab::simulate_lpa(300, means = rbind(rep(0, 4), rep(1.5, 4)),
      sds = rbind(rep(1, 4), c(0.7, 1.3, 0.8, 1.2)), proportions = c(0.6, 0.4), seed = seed))
    stats::setNames(d[paste0("indicator_", 1:4)], paste0("y", 1:4))
  }),
  eee_1v2 = list(base = 100000L, k0 = 1L, model = "EEE", make = function(seed) {
    set.seed(seed)
    z <- matrix(stats::rnorm(300 * 5), 300, 5) %*% chol(srl_cov)
    stats::setNames(as.data.frame(sweep(z, 2L, srl_mean, "+")), paste0("y", 1:5))
  }))
spec <- conditions[[condition]]

run_one <- function(replicate) {
  seed <- spec$base + replicate
  data <- spec$make(seed)
  vars <- names(data)
  quiet <- function(expr) withCallingHandlers(expr,
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    null <- quiet(lpa(data, vars, spec$k0, model = spec$model, seed = seed, search = "emEM"))
    alternative <- quiet(lpa(data, vars, spec$k0 + 1L, model = spec$model, seed = seed,
                             search = "emEM"))
    test <- quiet(bootstrap_lrt(null, alternative, iter = iter, seed = seed))
    data.frame(condition = condition, replicate = replicate, p_value = test$p_value,
               statistic = test$statistic, null_ll = null$log_likelihood,
               alternative_ll = alternative$log_likelihood,
               refit_search = paste(unique(test$refit_controls$search), collapse = "/"),
               error = NA_character_)
  }, error = function(e) data.frame(condition = condition, replicate = replicate,
                                    p_value = NA_real_, statistic = NA_real_,
                                    null_ll = NA_real_, alternative_ll = NA_real_,
                                    refit_search = NA_character_,
                                    error = paste(class(e)[1], conditionMessage(e))))
  result$seconds <- proc.time()[["elapsed"]] - started
  saveRDS(result, file.path(out, sprintf("%s_%03d.rds", condition, replicate)))
  result
}
invisible(mclapply(replicates, run_one, mc.cores = cores, mc.preschedule = FALSE))
cat(condition, "done\n")
