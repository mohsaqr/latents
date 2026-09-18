# Run from the project root: Rscript validation/synthetic-demo.R
source(file.path("R", "fit-ml-lpa.R"))
source(file.path("R", "gaussian-moments.R"))
source(file.path("R", "methods.R"))

set.seed(20260917)
n_groups <- 60L
group_sizes <- rep(c(16L, 20L, 24L), length.out = n_groups)
group_class <- rep(1:2, each = n_groups / 2L)
group_index <- rep(seq_len(n_groups), times = group_sizes)
profile_probability <- c(0.85, 0.15)[group_class[group_index]]
true_profile <- 1L + as.integer(runif(length(group_index)) > profile_probability)
synthetic_data <- data.frame(
  school_id = sprintf("school_%02d", group_index),
  reading = rnorm(length(group_index), c(-2, 2)[true_profile], 0.6),
  maths = rnorm(length(group_index), c(-1.5, 1.5)[true_profile], 0.8)
)
str(synthetic_data)
print(head(synthetic_data))
print(summary(synthetic_data))
print(vapply(synthetic_data, class, character(1)))
stopifnot(!anyNA(synthetic_data), !anyDuplicated(synthetic_data),
          nrow(synthetic_data) == sum(group_sizes))

elapsed <- system.time({
  fit <- fit_ml_lpa(synthetic_data, c("reading", "maths"), "school_id",
                    n_profiles = 2, n_group_classes = 2,
                    n_starts = 8, seed = 42)
})
str(fit, max.level = 1L)
print(dim(fit$subject_posteriors))
print(dim(fit$group_posteriors))
print(head(fit$subject_posteriors))
print(summary(fit))
print(elapsed)

# Align the arbitrary mixture labels before assessing recovery.
profile_order <- order(fit$means[, "reading"])
group_order <- order(fit$profile_probabilities[, profile_order[1L]], decreasing = TRUE)
individual_accuracy <- mean(match(fit$subject_profiles, profile_order) == true_profile)
group_accuracy <- mean(match(fit$group_classes, group_order) == group_class)
stopifnot(fit$converged, !fit$boundary,
          all(is.finite(fit$subject_posteriors)),
          all(is.finite(fit$group_posteriors)),
          max(abs(rowSums(fit$subject_posteriors) - 1)) < 1e-10,
          max(abs(rowSums(fit$group_posteriors) - 1)) < 1e-10,
          all(diff(fit$log_likelihood_history) >= -1e-7),
          individual_accuracy > 0.95, group_accuracy > 0.85)
cat(sprintf("Verified recovery: individual accuracy %.3f; group accuracy %.3f\n",
            individual_accuracy, group_accuracy))
dir.create("tmp", showWarnings = FALSE)
saveRDS(list(data = synthetic_data, true_profile = true_profile,
             true_group_class = group_class, fit = fit, elapsed = elapsed),
        file.path("tmp", "synthetic-demo.rds"))
