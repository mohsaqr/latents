# Pilot of the start-search benchmark (SEARCH-PLAN.md, section 6): how often
# each search misses the best known maximum, what it costs, and how often the
# solution it reports is a spike (a near-empty class or a variance near 0).
#
# Rscript equivalence/search/pilot.R <out_dir> <condition> <first> <last> <cores>
# Run under R_LIBS pointing at a build with the emEM search. PILOT_FITS, a
# comma-separated list of labels (e.g. "emem_10"), runs only those fits.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                value = TRUE))),
                 "..", "..", "validation", "fork-runtime.R"))
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
condition <- args[2]
replicates <- seq(as.integer(args[3]), as.integer(args[4]))
cores <- as.integer(args[5])
invisible(.validation_fork_cores(cores))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({
  library(latents)
  library(simulab)
  library(parallel)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "spike-measure.R"))
only <- strsplit(Sys.getenv("PILOT_FITS"), ",", fixed = TRUE)[[1]]

conditions <- list(
  lpa_vvi_k3 = list(family = "lpa", model = "VVI", n = 300, fit_k = 3:4, d = 4L,
                    means = rbind(c(0, 0, 0, 0), c(1.5, 1.5, 0, 0), c(0, 1.5, 1.5, 1.5)),
                    sds = rbind(rep(1, 4), rep(0.7, 4), rep(1.3, 4))),
  lpa_vvv_k2 = list(family = "lpa", model = "VVV", n = 300, fit_k = 2:3, d = 5L,
                    means = rbind(rep(0, 5), c(1.2, 1.2, 1.2, 0, 0)),
                    sds = rbind(rep(1, 5), rep(0.8, 5))),
  lpa_eee_k3 = list(family = "lpa", model = "EEE", n = 500, fit_k = 3:4, d = 5L,
                    means = rbind(rep(0, 5), c(1.3, 1.3, 0, 0, 0), c(0, 0, 1.3, 1.3, 1.3)),
                    sds = rbind(rep(1, 5), rep(1, 5), rep(1, 5))),
  lca_k3 = list(family = "lca", n = 500, fit_k = 3:4, d = 6L,
                yes = rbind(c(.85, .85, .85, .2, .2, .2), c(.2, .2, .85, .85, .85, .85),
                            c(.15, .15, .15, .15, .15, .15))),
  ml_lpa = list(family = "ml", model = "VVI", clusters = 60L, cluster_size = 8L,
                fit_k = 3L, fit_g = 2L, d = 4L,
                means = rbind(c(0, 0, 0, 0), c(1.6, 1.6, 0, 0), c(0, 1.6, 1.6, 1.6)),
                profile_probabilities = rbind(c(.6, .3, .1), c(.1, .3, .6)))
)
spec <- conditions[[condition]]

simulate <- function(spec, seed) {
  if (identical(spec$family, "lpa")) {
    d <- simulate_lpa(n = spec$n, means = spec$means, sds = spec$sds, seed = seed)
    return(d[grep("^indicator_", names(d))])
  }
  if (identical(spec$family, "lca")) {
    k <- nrow(spec$yes)
    probabilities <- array(c(1 - spec$yes, spec$yes), dim = c(k, ncol(spec$yes), 2L),
                           dimnames = list(class = paste0("c", seq_len(k)),
                                           indicator = paste0("item_", seq_len(ncol(spec$yes))),
                                           category = c("no", "yes")))
    d <- simulate_lca(n = spec$n, probabilities = probabilities,
                      category_labels = c("no", "yes"), seed = seed)
    return(d[grep("^item_", names(d))])
  }
  d <- simulate_ml_lpa(clusters = spec$clusters, cluster_size = spec$cluster_size,
                       means = spec$means, profile_probabilities = spec$profile_probabilities,
                       seed = seed)
  d
}

# The diagnostics a spike shows: the smallest expected class size, and the
# smallest class variance relative to the indicator's sample variance (on the
# standardized scale for full covariances). These benchmark measurements do
# not filter production solutions; current finalist tables may lack the metric.
describe <- function(fit, data, vars, seconds, label) {
  has_gaussian <- !identical(spec$family, "lca")
  if (inherits(fit, "error")) {
    return(list(row = data.frame(fit = label, log_likelihood = NA_real_, converged = NA,
                      boundary = NA, min_class = NA_real_, min_relative_variance = NA_real_,
                      replicated = NA_integer_, seconds = seconds,
                      error = conditionMessage(fit), has_gaussian = has_gaussian), finals = NULL))
  }
  sizes <- colSums(fit$subject_posteriors)
  numeric_vars <- vars[vapply(data[vars], is.numeric, TRUE)]
  overall <- vapply(data[numeric_vars], function(v) mean((v - mean(v))^2), numeric(1))
  ratio <- if (length(numeric_vars) == length(vars)) {
    spike_measure(fit, overall)
  } else NA_real_
  starts <- get_results(fit, "starts")
  finals <- if ("stage" %in% names(starts)) {
    # Missing Gaussian finalist metrics are unmeasured, while an LCA metric
    # is inapplicable. Preserve that distinction explicitly in saved records.
    if (!"min_relative_variance" %in% names(starts)) {
      starts$min_relative_variance <- NA_real_
      starts$admissible <- NA
    }
    starts$has_gaussian <- has_gaussian
    subset(starts, stage == "final", select = c(start, kind, log_likelihood, converged,
                                                min_relative_variance, admissible, has_gaussian))
  }
  list(row = data.frame(fit = label, log_likelihood = fit$log_likelihood,
             converged = isTRUE(fit$converged), boundary = isTRUE(fit$boundary),
             min_class = min(sizes), min_relative_variance = ratio,
             replicated = fit$n_best_replicated %||% NA_integer_,
             seconds = seconds, error = NA_character_, has_gaussian = has_gaussian),
       finals = if (!is.null(finals)) cbind(fit = label, finals))
}

run_one <- function(replicate) {
  seed <- 20261021L + 1000L * match(condition, names(conditions)) + replicate
  data <- simulate(spec, seed)
  vars <- if (identical(spec$family, "ml")) grep("^indicator_", names(data), value = TRUE) else
    names(data)
  ks <- spec$fit_k
  rows <- lapply(ks, function(k) {
    fit_with <- function(label, ...) {
      if (length(only) && !label %in% only) return(NULL)
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(withCallingHandlers(
        if (identical(spec$family, "lpa")) lpa(data, vars, k, model = spec$model, seed = seed, ...)
        else if (identical(spec$family, "lca")) lca(data, vars, k, seed = seed, ...)
        else multilpa(data, vars, "cluster", k, spec$fit_g, model = spec$model, seed = seed, ...),
        warning = function(w) invokeRestart("muffleWarning"),
        message = function(m) invokeRestart("muffleMessage")),
        error = function(e) e)
      describe(fit, data, vars, proc.time()[["elapsed"]] - started, label)
    }
    found <- list(
      fit_with("kmeans_10", search = "kmeans"),
      fit_with("kmeans_20", search = "kmeans", n_starts = 20),
      fit_with("kmeans_100", search = "kmeans", n_starts = 100),
      fit_with("emem_10", search = "emEM"),
      fit_with("emem_10_pool100", search = "emEM", n_screen = 100),
      fit_with("emem_20", search = "emEM", n_starts = 20),
      fit_with("reference_emem", search = "emEM", n_starts = 50, n_screen = 1000))
    found <- Filter(Negate(is.null), found)
    list(rows = cbind(condition = condition, replicate = replicate, seed = seed, k = k,
                      do.call(rbind, lapply(found, `[[`, "row"))),
         finals = cbind(condition = condition, replicate = replicate, k = k,
                        do.call(rbind, lapply(found, `[[`, "finals"))))
  })
  result <- list(rows = do.call(rbind, lapply(rows, `[[`, "rows")),
                 finals = do.call(rbind, lapply(rows, `[[`, "finals")))
  saveRDS(result, file.path(out, sprintf("%s_%03d.rds", condition, replicate)))
  result
}

invisible(mclapply(replicates, run_one, mc.cores = .validation_fork_cores(cores),
                   mc.preschedule = FALSE))
cat(condition, "done\n")
