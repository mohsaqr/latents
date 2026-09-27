# Monte Carlo recovery and interval coverage for the two-level mixture
# regression, the one mixture_regression() model with no external implementation to
# compare against. Run from the package root:
#   Rscript validation/mixture-regression-recovery.R
#
# Data are simulated from a fitted-to-truth two-level model (150 groups of
# six rows, two regression classes, two group classes, a row-level and a
# group-level membership covariate), refitted, and every estimate compared
# with its generating value: bias, and the coverage of the nominal 95%
# Wald interval under the observed-information and robust covariances.
pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260927)
replications <- 200L

truth_fit <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                    n_group_classes = 2, membership = ~ sleep,
                    group_membership = ~ motivation, n_starts = 5, seed = 1,
                    tol = 1e-10)
truth <- coef(truth_fit)
seeds <- sample.int(.Machine$integer.max, replications)

one_replication <- function(seed) {
  simulated <- study_hours
  simulated$score <- simulate(truth_fit, nsim = 1, seed = seed)$sim_1
  fit <- withCallingHandlers(
    mixture_regression(score ~ hours, simulated, 2, id = "student", n_group_classes = 2,
           membership = ~ sleep, group_membership = ~ motivation,
           n_starts = 3, seed = seed, tol = 1e-10),
    latents_degenerate_start = function(w) invokeRestart("muffleWarning"))
  # Align labels to the truth: of the 2 x 2 relabellings, keep the one whose
  # estimates are closest to the generating values. Without this, a
  # replication whose classes came out in the other order would count as a
  # large error, which is label switching and not estimation error.
  relabellings <- expand.grid(classes = 1:2, groups = 1:2)
  candidates <- lapply(seq_len(nrow(relabellings)), function(r) {
    class_order <- if (relabellings$classes[r] == 1L) 1:2 else 2:1
    group_order <- if (relabellings$groups[r] == 1L) 1:2 else 2:1
    relabelled <- fit
    relabelled$params <- latents:::.mixture_permute(fit$spec, fit$params,
                                                   class_order, group_order)
    relabelled$expectation <- latents:::.mixture_expectation(
      fit$spec, relabelled$params)
    relabelled
  })
  distance <- vapply(candidates, function(candidate) {
    sum((coef(candidate) - truth)^2)
  }, numeric(1))
  fit <- candidates[[which.min(distance)]]
  fit$inference <- NULL
  estimate <- coef(fit)
  observed <- sqrt(diag(vcov(fit, type = "observed")))
  robust <- sqrt(diag(vcov(fit, type = "robust")))
  data.frame(replication = seed, parameter = names(truth),
             truth = unname(truth), estimate = unname(estimate),
             se_observed = unname(observed), se_robust = unname(robust),
             converged = fit$converged,
             relabelled = which.min(distance) != 1L)
}

rows <- do.call(rbind, parallel::mclapply(seeds, one_replication,
                                          mc.cores = 4L, mc.set.seed = TRUE))
z <- stats::qnorm(0.975)
rows$covered_observed <- abs(rows$estimate - rows$truth) <= z * rows$se_observed
rows$covered_robust <- abs(rows$estimate - rows$truth) <= z * rows$se_robust
summary_table <- do.call(rbind, lapply(split(rows, rows$parameter), function(p) {
  data.frame(parameter = p$parameter[1L], truth = p$truth[1L],
             mean_estimate = mean(p$estimate),
             bias = mean(p$estimate - p$truth),
             empirical_sd = stats::sd(p$estimate),
             mean_se_observed = mean(p$se_observed),
             coverage_observed = mean(p$covered_observed),
             coverage_robust = mean(p$covered_robust))
}))
summary_table <- summary_table[match(names(truth), summary_table$parameter), ]
rownames(summary_table) <- NULL
cat(sprintf("replications: %d, converged: %d, relabelled to the truth: %d\n",
            nrow(rows) / length(truth), sum(rows$converged) / length(truth),
            sum(rows$relabelled) / length(truth)))
print(summary_table, digits = 3)
# Monte Carlo standard error of a 95% coverage estimate from R replications.
cat(sprintf("coverage Monte Carlo SE: %.3f\n",
            sqrt(0.95 * 0.05 / replications)))
utils::write.csv(summary_table, file.path("validation", "mixture-regression-recovery.csv"),
                 row.names = FALSE)
