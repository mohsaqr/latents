# Recovery and interval coverage of a membership slope when its covariate is
# missing at random, comparing complete-case analysis with multiple imputation
# pooled by `pool_imputations()`. Run from the package root (needs mice):
#   Rscript validation/pool-imputations-coverage.R
#
# Data: 60 groups of eight rows, two profiles whose membership has logit slope
# 1.2 on a row-level covariate `z`, two Gaussian indicators. `z` is deleted
# with probability plogis(-2 + 0.8 y1), so rows in the high-mean profile lose
# it more often: missing at random given the indicators, which the imputation
# model includes and complete-case analysis does not. Four analyses of each
# data set: complete cases; and ten imputations each from mice's predictive
# mean matching (`pmm`, its default) and Bayesian linear regression (`norm`),
# the group identifier excluded as a predictor, and from an approximate conditional
# distribution of z given the indicators under the generating model (`exact`,
# sampling-importance-resampling with the true parameters). The legacy `exact`
# label denotes a finite 400-proposal approximation, not an exact sampler.
# This oracle arm holds imputation parameters at their known values; it omits
# their estimation uncertainty and does not establish proper MI coverage.
# A linear imputation model can miss the latent structure, but this comparison
# alone does not isolate the cause of a coverage difference. Profiles
# are aligned to the generating labels by the first indicator's mean.
pkgload::load_all(".", quiet = TRUE)
stopifnot("this study needs the mice package" = requireNamespace("mice", quietly = TRUE))
RNGkind("L'Ecuyer-CMRG")
set.seed(20260927)
replications <- 200L
truth <- 1.2

replicate_once <- function(replication) {
  n <- 480L
  data <- data.frame(g = rep(seq_len(60L), each = 8L), z = stats::rnorm(n))
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(truth * data$z))
  data$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  data$y2 <- stats::rnorm(n, c(0, 2)[state])
  data$z[stats::runif(n) < stats::plogis(-2 + 0.8 * data$y1)] <- NA
  quiet <- function(expr) {
    withCallingHandlers(expr,
      latents_unconverged = function(condition) invokeRestart("muffleWarning"),
      latents_extreme_coefficients = function(condition) invokeRestart("muffleWarning"))
  }
  # The slope of the profile with the lower first-indicator mean, whatever the
  # fit called it: relabelling the two profiles negates the logit.
  slope_row <- function(table, means) {
    row <- table[table$level == "profile" & table$term == "z", ]
    if (means[1L, 1L] > means[2L, 1L]) {
      row$estimate <- -row$estimate
      bounds <- -c(row$conf_high, row$conf_low)
      row$conf_low <- bounds[1L]
      row$conf_high <- bounds[2L]
    }
    row
  }
  complete <- quiet(multilpa(data[!is.na(data$z), ], c("y1", "y2"), "g", 2L, 1L,
                             profile_covariates = "z", n_starts = 2, seed = 1))
  missing_rows <- which(is.na(data$z))
  predictors <- mice::make.predictorMatrix(data)
  predictors[, "g"] <- 0
  # z given the indicators under the generating model: prior N(0, 1) times the
  # mixture likelihood of the row's indicators, resampled from 400 draws.
  exact <- lapply(seq_len(10L), function(index) {
    completed <- data
    completed$z[missing_rows] <- vapply(missing_rows, function(row) {
      candidates <- stats::rnorm(400L)
      first <- stats::plogis(truth * candidates)
      weight <- first * stats::dnorm(data$y1[row], 0) * stats::dnorm(data$y2[row], 0) +
        (1 - first) * stats::dnorm(data$y1[row], 2.5) * stats::dnorm(data$y2[row], 2)
      candidates[sample.int(400L, 1L, prob = weight)]
    }, numeric(1))
    completed
  })
  imputations <- list(
    pmm = mice::mice(data, m = 10, predictorMatrix = predictors, printFlag = FALSE),
    norm = mice::mice(data, m = 10, predictorMatrix = predictors, method = "norm",
                      printFlag = FALSE),
    exact = exact)
  pooled_rows <- lapply(imputations, function(imputed) {
    pooled <- as.data.frame(quiet(pool_imputations(
      imputed, c("y1", "y2"), "g", 2L, n_group_classes = 1L,
      profile_covariates = "z", n_starts = 2, seed = 1)))
    means <- matrix(pooled$estimate[pooled$parameter == "mean" & pooled$term == "y1"],
                    ncol = 1L)
    slope_row(pooled, means)
  })
  rows <- c(list(complete_case = slope_row(parameter_inference(complete),
                                           complete$means)),
            pooled_rows)
  do.call(rbind, lapply(names(rows), function(method) {
    row <- rows[[method]]
    data.frame(replication = replication, method = method,
               estimate = row$estimate, standard_error = row$standard_error,
               covered = row$conf_low <= truth && truth <= row$conf_high)
  }))
}

# One row per method.
summarise_recovery <- function(results) {
  do.call(rbind, lapply(split(results, results$method), function(part) {
    data.frame(method = part$method[1L], usable = nrow(part),
               bias = mean(part$estimate) - truth,
               empirical_sd = stats::sd(part$estimate),
               mean_standard_error = mean(part$standard_error),
               coverage = mean(part$covered))
  }))
}

results <- do.call(rbind, parallel::mclapply(
  seq_len(replications), replicate_once, mc.set.seed = TRUE,
  mc.cores = max(1L, parallel::detectCores() - 2L)))
summary_table <- summarise_recovery(results)
print(summary_table, digits = 3, row.names = FALSE)
utils::write.csv(summary_table, "validation/pool-imputations-coverage.csv",
                 row.names = FALSE)
