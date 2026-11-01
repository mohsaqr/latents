# Size of the likelihood-ratio test of compare_models(test = "lrt") under the
# null. Run from the package root:
#   Rscript validation/lrt-size.R          # 300 data sets per arm, 2 cores
#   Rscript validation/lrt-size.R 5        # a quick pilot
#
# In every data set the restricted model is true, so the p-values should be
# uniform and a test at level a should reject in about a proportion a of the
# data sets. Declared before running: with 300 data sets the Monte Carlo
# standard error of a 5% rejection rate is 0.0126, so a rate between 0.025 and
# 0.075 (two standard errors) is consistent with the chi-square reference.
#   structure:   EEI against VVI, three Gaussian indicators, two profiles
#                (means 0 and 2 on each), 400 rows; df = 3.
#   correlation: VVI against VVV, the same design with unequal variances;
#                df = 6.
#   covariate:   lca() without and with a membership covariate that has no
#                effect, five binary items, two classes, 400 rows; df = 1.
#   invariance:  lta() with measurement "invariant" against "occasion", three
#                binary items, two states, two occasions, 400 persons; df = 6.
# Fits with a variance at its bound, refused pairs and failed fits are kept as
# rows with their reason, and the rates are computed over the usable p-values.
# Two socket workers (forked workers are unsafe with Apple's Accelerate BLAS),
# each with its own L'Ecuyer-CMRG stream and a fixed share of the data sets,
# so a rerun with two workers reproduces the results.
pkgload::load_all(".", quiet = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(arguments) > 0L) as.integer(arguments[[1L]]) else 300L

quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(condition) invokeRestart("muffleWarning"))
}

gaussian_rows <- function(standard_deviations) {
  state <- stats::rbinom(400L, 1L, 0.5) + 1L
  centre <- c(0, 2)[state]
  data.frame(y1 = stats::rnorm(400L, centre, standard_deviations[state, 1L]),
             y2 = stats::rnorm(400L, centre, standard_deviations[state, 2L]),
             y3 = stats::rnorm(400L, centre, standard_deviations[state, 3L]))
}

binary_answers <- function(state, n_items) {
  probability <- c(0.2, 0.8)[state]
  as.data.frame(lapply(stats::setNames(seq_len(n_items), paste0("u", seq_len(n_items))),
                       function(item) factor(stats::rbinom(length(state), 1L, probability),
                                             levels = 0:1)))
}

# The pair of fits for one arm: the restricted model first.
fit_pair <- function(arm) {
  if (arm %in% c("structure", "correlation")) {
    spread <- if (identical(arm, "structure")) matrix(1, 2L, 3L) else
      rbind(c(1, 0.7, 1.2), c(0.8, 1.3, 0.9))
    data <- gaussian_rows(spread)
    models <- if (identical(arm, "structure")) c("EEI", "VVI") else c("VVI", "VVV")
    return(lapply(models, function(model) {
      lpa(data, c("y1", "y2", "y3"), n_profiles = 2, model = model, n_starts = 5, seed = 1)
    }))
  }
  if (identical(arm, "covariate")) {
    data <- cbind(binary_answers(stats::rbinom(400L, 1L, 0.5) + 1L, 5L), z = stats::rnorm(400L))
    items <- paste0("u", 1:5)
    return(list(lca(data, items, n_classes = 2, n_starts = 5, seed = 1),
                lca(data, items, n_classes = 2, n_starts = 5, seed = 1,
                    profile_covariates = "z")))
  }
  first <- stats::rbinom(400L, 1L, 0.5) + 1L
  second <- ifelse(stats::runif(400L) < 0.8, first, 3L - first)
  panel <- rbind(cbind(person = seq_len(400L), occasion = 1L, binary_answers(first, 3L)),
                 cbind(person = seq_len(400L), occasion = 2L, binary_answers(second, 3L)))
  items <- paste0("u", 1:3)
  lapply(c("invariant", "occasion"), function(measurement) {
    lta(panel, items, "person", n_profiles = 2, time = "occasion", categorical = items,
        measurement = measurement, n_starts = 5, seed = 1)
  })
}

replicate_once <- function(replication, arm) {
  record <- function(statistic = NA_real_, df = NA_integer_, p_value = NA_real_, reason = NA) {
    data.frame(arm = arm, replication = replication, statistic = statistic, df = df,
               p_value = p_value, reason = reason)
  }
  fits <- tryCatch(quiet(fit_pair(arm)), error = function(condition) conditionMessage(condition))
  if (is.character(fits)) return(record(reason = paste("fit failed:", fits)))
  if (any(vapply(fits, function(fit) isTRUE(fit$boundary), logical(1)))) {
    return(record(reason = "variance at bound"))
  }
  tested <- tryCatch(quiet(compare_models(restricted = fits[[1L]], general = fits[[2L]],
                                          test = "lrt")),
                     error = function(condition) class(condition)[1L])
  if (is.character(tested)) return(record(reason = tested))
  row <- subset(as.data.frame(tested), model == "general")
  with(row, record(lr_statistic, lr_df, lr_p_value))
}

# One row per arm: how many data sets gave a p-value, and how often it fell
# at or below 0.05 and 0.10.
summarise_size <- function(results) {
  do.call(rbind, lapply(split(results, results$arm), function(part) {
    usable <- subset(part, !is.na(p_value))
    with(usable, data.frame(
      arm = part$arm[1L], df = stats::median(df), datasets = nrow(part),
      usable = nrow(usable), reject_05 = mean(p_value <= 0.05),
      reject_10 = mean(p_value <= 0.10), mean_statistic = mean(statistic),
      within_two_se = abs(mean(p_value <= 0.05) - 0.05) <=
        2 * sqrt(0.05 * 0.95 / nrow(usable))))
  }))
}

arms <- c("structure", "correlation", "covariate", "invariance")
started <- proc.time()[["elapsed"]]
cluster <- parallel::makePSOCKcluster(2L)
parallel::clusterSetRNGStream(cluster, 20261031)
invisible(parallel::clusterEvalQ(cluster, {
  Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  pkgload::load_all(".", quiet = TRUE)
}))
parallel::clusterExport(cluster, c("quiet", "gaussian_rows", "binary_answers", "fit_pair",
                                   "replicate_once"))
results <- do.call(rbind, lapply(arms, function(arm) {
  do.call(rbind, parallel::parLapply(cluster, seq_len(datasets), replicate_once, arm = arm))
}))
parallel::stopCluster(cluster)
summary_table <- summarise_size(results)
print(summary_table, digits = 3, row.names = FALSE)
print(table(arm = results$arm, reason = results$reason, useNA = "no"))
cat(sprintf("%.1f min\n", (proc.time()[["elapsed"]] - started) / 60))
if (datasets >= 100L) {
  utils::write.csv(results, "validation/lrt-size-replicates.csv", row.names = FALSE)
  utils::write.csv(summary_table, "validation/lrt-size.csv", row.names = FALSE)
}
