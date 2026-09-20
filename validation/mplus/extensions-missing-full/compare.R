# Reproduce comparisons and compact fixtures from actual retained Mplus outputs.
# Run from the project root: Rscript validation/mplus/extensions-missing-full/compare.R
# Load the whole package rather than naming source files. `R/fit-ml-lpa.R`
# stopped existing when the sources were reorganised, and a comparison script
# should not depend on the internal file layout.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
artifact_dir <- file.path("validation", "mplus", "extensions-missing-full")
synthetic <- read.table(file.path(artifact_dir, "synthetic.dat"),
                        col.names = c("y1", "y2", "clus", "id"), na.strings = "-999")
str(synthetic)
print(head(synthetic))
print(summary(synthetic))
print(vapply(synthetic, class, character(1)))
stopifnot(nrow(synthetic) == 1200, !anyDuplicated(synthetic$id),
          sum(is.na(synthetic$y1)) == 240, sum(is.na(synthetic$y2)) == 138,
          !any(is.na(synthetic$y1) & is.na(synthetic$y2)))
comparisons <- lapply(c("varying", "equal"), function(variance_model) {
  stopifnot(is.character(variance_model), length(variance_model) == 1L)
  prefix <- if (variance_model == "varying") "full" else "equal"
  output <- readLines(file.path(artifact_dir, paste0(prefix, ".out")))
  stopifnot(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", output)),
            any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", output)),
            !any(grepl("*** ERROR", output, fixed = TRUE)))
  values <- scan(file.path(artifact_dir, paste0(prefix, "-results.dat")), quiet = TRUE)
  n_parameters <- if (variance_model == "varying") 13L else 10L
  stopifnot(length(values) == 2L * n_parameters + 7L,
            values[2L * n_parameters + 1L] == n_parameters)
  # Parameter order is inspected directly in TECH1, including off-diagonal theta.
  means <- rbind(values[1:2], values[6:7])
  covariance_1 <- values[c(3, 4, 4, 5)]
  covariance_2 <- if (variance_model == "varying") values[c(8, 9, 9, 10)] else covariance_1
  covariances <- array(c(covariance_1, covariance_2), c(2L, 2L, 2L))
  group_logit <- values[n_parameters - 2L]
  intercept <- values[n_parameters - 1L]
  slope <- values[n_parameters]
  profile_probability <- plogis(c(intercept + slope, intercept))
  profile_probabilities <- cbind(profile_probability, 1 - profile_probability)
  group_probabilities <- c(plogis(group_logit), 1 - plogis(group_logit))
  fit <- multilpa(synthetic, c("y1", "y2"), "clus", 2, 2,
                    covariance_model = "full", missing = "fiml",
                    variance_model = variance_model, n_starts = 20, seed = 20260917,
                    tol = 1e-13, max_iter = 10000)
  profile_order <- order(fit$means[, "y1"])[rank(means[, 1])]
  group_order <- order(fit$profile_probabilities[, profile_order[1L]])[
    rank(profile_probabilities[, 1])]
  saved <- read.table(file.path(artifact_dir, paste0(prefix, "-posteriors.dat")),
    na.strings = "*", col.names = c("y1", "y2", "p11", "p12", "p21", "p22",
                                    "cb", "cw", "joint", "id", "clus"))
  saved <- saved[match(synthetic$id, saved$id), ]
  stopifnot(identical(saved$id, synthetic$id), !anyDuplicated(saved$id))
  subject_posteriors <- cbind(saved$p11 + saved$p21, saved$p12 + saved$p22)
  group_posteriors <- cbind(saved$p11 + saved$p12, saved$p21 + saved$p22)
  group_posteriors <- group_posteriors[match(fit$group_values, saved$clus), ]
  differences <- c(log_likelihood = abs(fit$log_likelihood - values[2L * n_parameters + 2L]),
    means = max(abs(fit$means[profile_order, ] - means)),
    covariances = max(abs(fit$covariances[, , profile_order] - covariances)),
    profile_probabilities = max(abs(fit$profile_probabilities[group_order, profile_order] - profile_probabilities)),
    group_probabilities = max(abs(fit$group_probabilities[group_order] - group_probabilities)),
    subject_posteriors = max(abs(fit$subject_posteriors[, profile_order] - subject_posteriors)),
    group_posteriors = max(abs(fit$group_posteriors[, group_order] - group_posteriors)))
  print(variance_model)
  print(differences, digits = 12)
  stopifnot(fit$converged, !fit$boundary, all(differences < 1e-4), fit$n_parameters == n_parameters)
  reference <- list(data = synthetic, covariance_model = "full", missing = "fiml",
    variance_model = variance_model, means = means, covariances = covariances,
    profile_probabilities = profile_probabilities, group_probabilities = group_probabilities,
    subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
    n_parameters = n_parameters, mplus_standard_errors = values[n_parameters + seq_len(n_parameters)],
    log_likelihood = values[2L * n_parameters + 2L],
    aic = values[2L * n_parameters + 3L], bic_individual = values[2L * n_parameters + 4L],
    provenance = list(oracle = "Genuine Mplus 9 Demo: full covariance with missing indicators",
      version = output[1L], date = "2026-09-17", parameter_order_source = "Retained TECH1",
      results_precision = "eight significant digits", posterior_precision = "twelve decimal places",
      md5 = tools::md5sum(file.path(artifact_dir, c("synthetic.dat", paste0(prefix,
        c(".inp", ".out", "-results.dat", "-posteriors.dat")))))))
  saveRDS(reference, file.path("tests", "fixtures", "mplus",
                               paste0("twolevel-missing-full-", variance_model, ".rds")))
  list(variance_model = variance_model, differences = differences, fit = fit)
})
saveRDS(comparisons, file.path(artifact_dir, "comparison.rds"))
