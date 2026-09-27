# Monte Carlo recovery and interval coverage for profile-covariate slopes that
# differ across group classes (`multilpa(profile_slopes = "group_class")`).
# Run from the package root:
#   Rscript validation/covariate-slopes-coverage.R
#
# Data: groups of eight rows in two latent group classes of equal share. The
# logit of profile 1 is 0.5 + 1.5 z in group class 1 and -0.5 - 1.5 z in group
# class 2, so the covariate's effect reverses between the two kinds of group.
# Two Gaussian indicators. Labels are aligned to the generating ones by the
# first indicator's mean (profiles) and by the sign of the slope (group
# classes). Each slope is estimated from the groups in its class only, so the
# study is run at two group counts. Coverage is of the nominal 95% Wald
# interval under the observed information and under the cluster-robust
# sandwich.
pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260927)
replications <- 300L
truth <- c(intercept_1 = 0.5, intercept_2 = -0.5, slope_1 = 1.5, slope_2 = -1.5)

replicate_once <- function(replication, n_groups) {
  size <- 8L
  n <- n_groups * size
  data <- data.frame(g = rep(seq_len(n_groups), each = size), z = stats::rnorm(n))
  group_class <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size) + 1L
  logit <- truth[c("intercept_1", "intercept_2")][group_class] +
    truth[c("slope_1", "slope_2")][group_class] * data$z
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(logit))
  data$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  data$y2 <- stats::rnorm(n, c(0, 2)[state])
  fit <- withCallingHandlers(
    multilpa(data, c("y1", "y2"), "g", 2L, 2L, profile_covariates = "z",
             profile_slopes = "group_class", n_starts = 3, seed = 1),
    latents_unconverged = function(condition) invokeRestart("muffleWarning"),
    latents_extreme_coefficients = function(condition) invokeRestart("muffleWarning"))
  if (!fit$converged) return(NULL)
  observed <- parameter_inference(fit)
  robust <- parameter_inference(fit, vcov_type = "robust")
  # Profile 1 of the truth has the lower first-indicator mean; flipping the
  # profiles negates every profile-1 logit.
  profile_sign <- if (fit$means[1L, 1L] > fit$means[2L, 1L]) -1 else 1
  terms <- c("group_class_1", "group_class_2", "z:group_class_1", "z:group_class_2")
  pick <- function(table) {
    rows <- table[table$level == "profile" & table$term %in% terms, ]
    rows[match(terms, rows$term), ]
  }
  estimate_rows <- pick(observed)
  estimates <- profile_sign * estimate_rows$estimate
  # Group class 1 of the truth is the one with the positive slope.
  order <- if (estimates[3L] > estimates[4L]) c(1L, 2L, 3L, 4L) else c(2L, 1L, 4L, 3L)
  covers <- function(table) {
    rows <- pick(table)[order, ]
    bounds <- cbind(profile_sign * rows$conf_low, profile_sign * rows$conf_high)
    low <- pmin(bounds[, 1L], bounds[, 2L])
    high <- pmax(bounds[, 1L], bounds[, 2L])
    low <= truth & truth <= high
  }
  data.frame(n_groups = n_groups, replication = replication,
             parameter = names(truth), estimate = estimates[order],
             standard_error = estimate_rows$standard_error[order],
             covered_observed = covers(observed), covered_robust = covers(robust))
}

results <- do.call(rbind, lapply(c(60L, 120L), function(n_groups) {
  do.call(rbind, parallel::mclapply(seq_len(replications), replicate_once,
                                    n_groups = n_groups, mc.set.seed = TRUE,
                                    mc.cores = max(1L, parallel::detectCores() - 2L)))
}))
# One row per group count and parameter, ordered by both.
summarise_coverage <- function(results) {
  table <- do.call(rbind, lapply(
    split(results, list(results$n_groups, results$parameter), drop = TRUE),
    function(part) {
      data.frame(n_groups = part$n_groups[1L], parameter = part$parameter[1L],
                 usable = nrow(part),
                 bias = mean(part$estimate) - truth[[part$parameter[1L]]],
                 empirical_sd = stats::sd(part$estimate),
                 mean_standard_error = mean(part$standard_error),
                 coverage_observed = mean(part$covered_observed),
                 coverage_robust = mean(part$covered_robust))
    }))
  table[order(table$n_groups, table$parameter), ]
}
summary_table <- summarise_coverage(results)
print(summary_table, digits = 3, row.names = FALSE)
utils::write.csv(summary_table, "validation/covariate-slopes-coverage.csv",
                 row.names = FALSE)
