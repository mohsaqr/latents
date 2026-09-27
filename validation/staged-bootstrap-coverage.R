# Interval coverage for a staged fit (`fit_staged()`): the conditional Wald
# errors, which treat the first-stage measurement as known, against the
# two-stage bootstrap (`parameter_inference(method = "bootstrap")`), which
# refits both stages on every resample. Run from the package root:
#   Rscript validation/staged-bootstrap-coverage.R
#
# Data: 80 groups of eight rows in two latent group classes of equal share;
# profile 1 has probability 0.8 in group class 1 and 0.25 in group class 2;
# two Gaussian indicators with profile means (0, 0) and (1.6, 1.4), unit
# variances. The measurement model is correctly specified at the first stage,
# so the staged estimator is consistent and both intervals aim at the same
# values. Profiles are aligned by the first indicator's mean, group classes by
# their probability of profile 1.
pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260927)
replications <- 150L
bootstrap_iter <- 150L
truth <- c(high_profile_1 = 0.8, low_profile_1 = 0.25, high_share = 0.5)

replicate_once <- function(replication) {
  n_groups <- 80L
  size <- 8L
  n <- n_groups * size
  group_class <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size) + 1L
  profile <- 2L - stats::rbinom(n, 1L, c(0.8, 0.25)[group_class])
  data <- data.frame(g = rep(seq_len(n_groups), each = size),
                     y1 = stats::rnorm(n, c(0, 1.6)[profile]),
                     y2 = stats::rnorm(n, c(0, 1.4)[profile]))
  quiet <- function(expr) {
    withCallingHandlers(expr,
      latents_unconverged = function(condition) invokeRestart("muffleWarning"),
      latents_bootstrap_dropped = function(condition) invokeRestart("muffleWarning"))
  }
  fit <- quiet(fit_staged(data, c("y1", "y2"), "g", 2L, 2L, n_starts = 3, seed = 1))
  if (!fit$converged) return(NULL)
  wald <- tryCatch(parameter_inference(fit),
                   latents_boundary_fit = function(condition) NULL,
                   latents_singular_information = function(condition) NULL)
  if (is.null(wald)) return(NULL)
  bootstrap <- quiet(parameter_inference(fit, method = "bootstrap",
                                         iter = bootstrap_iter, n_starts = 2))
  # Profile 1 of the truth has the lower first-indicator mean.
  first <- if (fit$means[1L, 1L] < fit$means[2L, 1L]) "profile_1" else "profile_2"
  profile_rows <- function(table, group_class) {
    table[table$level == "profile" & table$outcome == first &
            table$term == group_class, ]
  }
  group_row <- function(table, group_class) {
    table[table$level == "group" & table$outcome == group_class, ]
  }
  shares <- vapply(c("group_class_1", "group_class_2"), function(group_class) {
    profile_rows(wald, group_class)$estimate
  }, numeric(1))
  high <- names(which.max(shares))
  low <- names(which.min(shares))
  covers <- function(table) {
    rows <- rbind(profile_rows(table, high), profile_rows(table, low),
                  group_row(table, high))
    rows$conf_low <= truth & truth <= rows$conf_high
  }
  standard_errors <- function(table) {
    rbind(profile_rows(table, high), profile_rows(table, low),
          group_row(table, high))$standard_error
  }
  data.frame(replication = replication, parameter = names(truth),
             covered_wald = covers(wald), covered_bootstrap = covers(bootstrap),
             se_wald = standard_errors(wald),
             se_bootstrap = standard_errors(bootstrap),
             estimate = rbind(profile_rows(wald, high), profile_rows(wald, low),
                              group_row(wald, high))$estimate)
}

# One row per parameter: coverage of both intervals and the spread they claim
# against the spread the estimates actually have.
summarise_coverage <- function(results) {
  do.call(rbind, lapply(split(results, results$parameter), function(part) {
    data.frame(parameter = part$parameter[1L], usable = nrow(part),
               empirical_sd = stats::sd(part$estimate),
               mean_se_wald = mean(part$se_wald),
               mean_se_bootstrap = mean(part$se_bootstrap),
               coverage_wald = mean(part$covered_wald),
               coverage_bootstrap = mean(part$covered_bootstrap))
  }))
}

results <- do.call(rbind, parallel::mclapply(
  seq_len(replications), replicate_once, mc.set.seed = TRUE,
  mc.cores = max(1L, parallel::detectCores() - 2L)))
summary_table <- summarise_coverage(results)
print(summary_table, digits = 3, row.names = FALSE)
utils::write.csv(summary_table, "validation/staged-bootstrap-coverage.csv",
                 row.names = FALSE)
