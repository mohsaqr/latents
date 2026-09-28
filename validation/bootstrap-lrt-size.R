# Size of the parametric bootstrap likelihood-ratio test on the two model
# kinds it newly accepts: fits with missing indicators (`missing = "fiml"`)
# and membership-covariate fits. Run from the package root:
#   Rscript validation/bootstrap-lrt-size.R
#
# In every data set the two-profile null is true, and it is tested against
# three profiles, so the p-values should be uniform: a test at level a should
# reject in about a proportion a of the data sets.
#   fiml:      50 groups of eight, two Gaussian indicators with profile means
#              (0, 0) and (2, 1.6), 15% of each indicator missing completely at
#              random; the bootstrap carries the observed pattern into every
#              replicate.
#   covariate: the same measurement on complete data, with profile membership
#              depending on a row covariate (logit slope 1.2); the bootstrap
#              holds the covariate at its observed values.
pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260928)
datasets <- 100L
bootstrap_iter <- 39L

replicate_once <- function(replication, arm) {
  n <- 400L
  data <- data.frame(g = rep(seq_len(50L), each = 8L), z = stats::rnorm(n))
  slope <- if (identical(arm, "covariate")) 1.2 else 0
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(slope * data$z))
  data$y1 <- stats::rnorm(n, c(0, 2)[state])
  data$y2 <- stats::rnorm(n, c(0, 1.6)[state])
  if (identical(arm, "fiml")) {
    data$y1[stats::runif(n) < 0.15] <- NA
    data$y2[stats::runif(n) < 0.15] <- NA
  }
  arguments <- list(data = data, vars = c("y1", "y2"), id = "g",
                    n_group_classes = 1L, n_starts = 3, seed = 1,
                    missing = if (identical(arm, "fiml")) "fiml" else "error")
  if (identical(arm, "covariate")) arguments$profile_covariates <- "z"
  quiet <- function(expr) {
    withCallingHandlers(expr, warning = function(condition) {
      invokeRestart("muffleWarning")
    })
  }
  null <- quiet(do.call(multilpa, c(arguments, list(n_profiles = 2L))))
  alternative <- quiet(do.call(multilpa, c(arguments, list(n_profiles = 3L))))
  test <- tryCatch(
    quiet(bootstrap_lrt(null, alternative, iter = bootstrap_iter, n_starts = 3,
                        max_iter = 3000)),
    error = function(condition) NULL)
  if (is.null(test)) {
    return(data.frame(arm = arm, replication = replication, p_value = NA_real_,
                      reason = "refused"))
  }
  data.frame(arm = arm, replication = replication,
             p_value = as.data.frame(test)$p_value,
             reason = if (is.na(as.data.frame(test)$p_value)) "failed replicates" else NA)
}

# One row per arm: how many data sets gave a p-value, and how often it fell
# at or below 0.05 and 0.10.
summarise_size <- function(results) {
  do.call(rbind, lapply(split(results, results$arm), function(part) {
    p <- part$p_value[!is.na(part$p_value)]
    data.frame(arm = part$arm[1L], datasets = nrow(part), usable = length(p),
               reject_05 = mean(p <= 0.05), reject_10 = mean(p <= 0.10),
               mean_p = mean(p))
  }))
}

results <- do.call(rbind, lapply(c("fiml", "covariate"), function(arm) {
  do.call(rbind, parallel::mclapply(seq_len(datasets), replicate_once, arm = arm,
                                    mc.set.seed = TRUE,
                                    mc.cores = max(1L, parallel::detectCores() - 2L)))
}))
summary_table <- summarise_size(results)
print(summary_table, digits = 3, row.names = FALSE)
utils::write.csv(summary_table, "validation/bootstrap-lrt-size.csv", row.names = FALSE)
