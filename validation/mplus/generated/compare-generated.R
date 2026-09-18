# Compare native R estimation with retained, independently fitted Mplus outputs.
# Run from the project root: Rscript validation/mplus/generated/compare-generated.R
source(file.path("R", "fit-ml-lpa.R"))
source(file.path("R", "gaussian-moments.R"))
source(file.path("R", "methods.R"))
artifact_dir <- file.path("validation", "mplus", "generated")
synthetic_data <- read.table(file.path(artifact_dir, "synthetic.dat"),
  col.names = c("y1", "y2", "clus", "id"))
str(synthetic_data)
print(head(synthetic_data))
print(summary(synthetic_data))
print(vapply(synthetic_data, class, character(1)))
stopifnot(nrow(synthetic_data) == 1200L, !anyNA(synthetic_data),
          !anyDuplicated(synthetic_data$id), length(unique(synthetic_data$clus)) == 60L)

# Decode the explicitly inspected TECH1 parameter ordering and align both levels.
comparisons <- lapply(c("varying", "equal"), function(variance_model) {
  stopifnot(is.character(variance_model), length(variance_model) == 1L)
  fit <- multilpa(synthetic_data, c("y1", "y2"), "clus", 2L, 2L,
    variance_model = variance_model, n_starts = 30L, max_iter = 5000L,
    tol = 1e-12, seed = 20260917)
  output <- readLines(file.path(artifact_dir, paste0(variance_model, ".out")))
  stopifnot(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", output)),
    any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", output)),
    !any(grepl("*** ERROR", output, fixed = TRUE)))
  values <- scan(file.path(artifact_dir, paste0(variance_model, "-results.dat")),
    quiet = TRUE)
  n_parameters <- if (variance_model == "varying") 11L else 9L
  stopifnot(length(values) == 2L * n_parameters + 7L,
    values[2L * n_parameters + 1L] == n_parameters,
    fit$n_parameters == n_parameters)
  mplus_means <- rbind(values[c(1, 2)], values[c(5, 6)])
  mplus_variances <- if (variance_model == "varying") {
    rbind(values[c(3, 4)], values[c(7, 8)])
  } else rbind(values[c(3, 4)], values[c(3, 4)])
  group_logit <- values[n_parameters - 2L]
  intercept <- values[n_parameters - 1L]
  slope <- values[n_parameters]
  profile_probability <- plogis(c(intercept + slope, intercept))
  mplus_profile <- cbind(profile_probability, 1 - profile_probability)
  mplus_group <- c(plogis(group_logit), 1 - plogis(group_logit))
  # Two profiles: match by the first indicator, whose means are distinct here.
  profile_order <- order(fit$means[, "y1"])[rank(mplus_means[, 1])]
  aligned_profile <- fit$profile_probabilities[, profile_order, drop = FALSE]
  group_order <- order(aligned_profile[, 1])[rank(mplus_profile[, 1])]
  mplus_saved <- read.table(file.path(artifact_dir,
    paste0(variance_model, "-posteriors.dat")), col.names = c("y1", "y2",
      "p11", "p12", "p21", "p22", "cb", "cw", "joint", "id", "clus"))
  mplus_saved <- mplus_saved[match(synthetic_data$id, mplus_saved$id), ]
  stopifnot(identical(mplus_saved$id, synthetic_data$id), !anyNA(mplus_saved),
    !anyDuplicated(mplus_saved$id))
  subject_posteriors <- cbind(mplus_saved$p11 + mplus_saved$p21,
                             mplus_saved$p12 + mplus_saved$p22)
  group_posteriors_by_subject <- cbind(mplus_saved$p11 + mplus_saved$p12,
                                      mplus_saved$p21 + mplus_saved$p22)
  first_group_rows <- match(unique(synthetic_data$clus), synthetic_data$clus)
  differences <- c(
    log_likelihood = abs(fit$log_likelihood - values[2L * n_parameters + 2L]),
    means = max(abs(fit$means[profile_order, ] - mplus_means)),
    variances = max(abs(fit$variances[profile_order, ] - mplus_variances)),
    profile_probabilities = max(abs(aligned_profile[group_order, ] - mplus_profile)),
    group_probabilities = max(abs(fit$group_probabilities[group_order] - mplus_group)),
    subject_posteriors = max(abs(fit$subject_posteriors[, profile_order] - subject_posteriors)),
    group_posteriors = max(abs(fit$group_posteriors[, group_order] -
      group_posteriors_by_subject[first_group_rows, ])),
    aic = abs(fit$aic - values[2L * n_parameters + 3L]),
    bic_individual = abs(fit$bic_individual - values[2L * n_parameters + 4L]))
  print(variance_model)
  print(differences, digits = 10)
  stopifnot(fit$converged, differences[c("log_likelihood", "aic", "bic_individual")] < 0.0001,
    differences[c("means", "variances", "profile_probabilities", "group_probabilities")] < 0.000001,
    differences[c("subject_posteriors", "group_posteriors")] <= 0.000001)
  reference <- list(
    data = synthetic_data,
    variance_model = variance_model,
    means = mplus_means,
    variances = mplus_variances,
    profile_probabilities = mplus_profile,
    group_probabilities = mplus_group,
    subject_posteriors = subject_posteriors,
    group_posteriors = group_posteriors_by_subject[first_group_rows, ],
    n_parameters = values[2L * n_parameters + 1L],
    log_likelihood = values[2L * n_parameters + 2L],
    aic = values[2L * n_parameters + 3L],
    bic_individual = values[2L * n_parameters + 4L],
    provenance = list(
      oracle = "Local Mplus9 Demo actual run on synthetic data",
      version = output[1],
      results_precision = "eight significant digits",
      posterior_precision = "twelve decimal places",
      parameter_order_source = "TECH1 in retained .out",
      date = "2026-09-17",
      md5 = tools::md5sum(file.path(artifact_dir, c("synthetic.dat",
        paste0(variance_model, c(".inp", ".out", "-results.dat", "-posteriors.dat")))))))
  saveRDS(reference, file.path("tests", "fixtures", "mplus",
    paste0("twolevel-synthetic-", variance_model, ".rds")))
  list(variance_model = variance_model, fit = fit, differences = differences,
    mplus_values = values, profile_order = profile_order, group_order = group_order)
})
names(comparisons) <- c("varying", "equal")
saveRDS(comparisons, file.path(artifact_dir, "comparison.rds"))
summary_table <- do.call(rbind, lapply(comparisons, function(comparison) {
  stopifnot(is.list(comparison), is.numeric(comparison$differences))
  data.frame(variance_model = comparison$variance_model,
    quantity = names(comparison$differences),
    absolute_difference = unname(comparison$differences))
}))
write.csv(summary_table, file.path(artifact_dir, "comparison.csv"), row.names = FALSE)
