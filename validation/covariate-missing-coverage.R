# Monte Carlo interval coverage for a membership-covariate fit with missing
# indicators (`multilpa(profile_covariates = , missing = "fiml")`), against
# the same design with complete data. Run from the package root:
#   Rscript validation/covariate-missing-coverage.R
#
# Data: 60 groups of eight rows, two profiles whose membership has logit
# slope 1 on a row-level covariate, two Gaussian indicators and one binary
# categorical indicator. Each indicator value is deleted completely at random
# with probability `share`. Profiles are aligned to the generating labels by
# the first indicator's mean. Coverage is of the nominal 95% Wald interval
# under the observed information.
pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260927)
replications <- 500L
truth <- c(slope = 1, mean_y1 = 0, probability_a = 0.7)

replicate_once <- function(replication, share) {
  n_groups <- 60L
  size <- 8L
  n <- n_groups * size
  data <- data.frame(g = rep(seq_len(n_groups), each = size), z = stats::rnorm(n))
  profile <- 2L - stats::rbinom(n, 1L, stats::plogis(truth[["slope"]] * data$z))
  data$y1 <- stats::rnorm(n, c(0, 2.5)[profile])
  data$y2 <- stats::rnorm(n, c(0, 2)[profile])
  data$q <- ifelse(stats::runif(n) < c(0.7, 0.2)[profile], "a", "b")
  blank <- function(value) replace(value, stats::runif(n) < share, NA)
  data[c("y1", "y2", "q")] <- lapply(data[c("y1", "y2", "q")], blank)
  fit <- withCallingHandlers(
    multilpa(data, c("y1", "y2", "q"), "g", 2L, 1L, categorical = "q",
             missing = if (share > 0) "fiml" else "error",
             profile_covariates = "z", n_starts = 2, seed = 1),
    latents_unconverged = function(condition) invokeRestart("muffleWarning"))
  if (!fit$converged) return(NULL)
  inference <- parameter_inference(fit)
  flipped <- fit$means[1L, 1L] > fit$means[2L, 1L]
  first <- if (flipped) "profile_2" else "profile_1"
  pick <- function(outcome, term, parameter) {
    inference[inference$outcome == outcome & inference$term == term &
                inference$parameter == parameter, ]
  }
  covers <- function(row, value, sign = 1) {
    bounds <- sort(sign * c(row$conf_low, row$conf_high))
    bounds[1L] <= value && value <= bounds[2L]
  }
  slope <- pick("profile_1", "z", "coefficient")
  data.frame(share = share, replication = replication,
             slope = covers(slope, truth[["slope"]], if (flipped) -1 else 1),
             mean_y1 = covers(pick(first, "y1", "mean"), truth[["mean_y1"]]),
             probability_a = covers(pick(first, "q:a", "response"),
                                    truth[["probability_a"]]),
             slope_estimate = if (flipped) -slope$estimate else slope$estimate,
             slope_standard_error = slope$standard_error)
}

results <- do.call(rbind, lapply(c(0, 0.2), function(share) {
  do.call(rbind, parallel::mclapply(seq_len(replications), replicate_once,
                                    share = share, mc.set.seed = TRUE,
                                    mc.cores = max(1L, parallel::detectCores() - 2L)))
}))
summary_table <- do.call(rbind, lapply(split(results, results$share), function(part) {
  data.frame(share = part$share[1L], usable = nrow(part),
             slope = mean(part$slope), mean_y1 = mean(part$mean_y1),
             probability_a = mean(part$probability_a),
             slope_bias = mean(part$slope_estimate) - truth[["slope"]],
             slope_sd = stats::sd(part$slope_estimate),
             slope_mean_standard_error = mean(part$slope_standard_error))
}))
print(summary_table, digits = 3, row.names = FALSE)
utils::write.csv(summary_table, "validation/covariate-missing-coverage.csv",
                 row.names = FALSE)
