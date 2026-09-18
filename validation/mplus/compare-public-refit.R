# Run from the project root after running the retained Mplus input files.
# This is a NEW matched LPA analysis on Example 10.4 data, not a reproduction
# of the original published CFA-mixture model (which this package cannot fit).
source(file.path("R", "fit-ml-lpa.R"))
source(file.path("R", "gaussian-moments.R"))
source(file.path("R", "methods.R"))
artifact_dir <- file.path("validation", "mplus", "public-refit")
public_dir <- file.path("validation", "mplus", "public")
indicator_names <- paste0("y", seq_len(5L))
data <- read.table(file.path(artifact_dir, "ex104-lpa.dat"),
                   col.names = c(indicator_names, "clus", "id"))
original <- read.table(file.path(public_dir, "ex10.4.dat"),
                       col.names = c(indicator_names, "generating_class", "clus"))
str(data)
print(head(data))
print(summary(data))
print(vapply(data, class, character(1)))
stopifnot(nrow(data) == 1000L, length(unique(data$clus)) == 110L,
          !anyNA(data), !anyDuplicated(data$id),
          isTRUE(all.equal(data[, c(indicator_names, "clus")],
                           original[, c(indicator_names, "clus")], tolerance = 0)))

comparisons <- lapply(c("varying", "equal"), function(variance_model) {
  stopifnot(is.character(variance_model), length(variance_model) == 1L)
  prefix <- paste0("public-", variance_model)
  output_file <- file.path(artifact_dir, paste0(prefix, ".out"))
  output <- readLines(output_file)
  stopifnot(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", output)),
            any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", output)),
            !any(grepl("*** ERROR", output, fixed = TRUE)))
  values <- scan(file.path(artifact_dir, paste0(prefix, "-results.dat")), quiet = TRUE)
  n_parameters <- if (variance_model == "varying") 23L else 18L
  stopifnot(length(values) == 2L * n_parameters + 7L,
            values[2L * n_parameters + 1L] == n_parameters)
  # TECH1 explicitly lists these parameter orders for the retained input files:
  # CW1 means, CW1 variances, CW2 means, optional CW2 variances, then three logits.
  parameter_names <- c(paste0("mu1_", indicator_names), paste0("var1_", indicator_names),
                        paste0("mu2_", indicator_names),
                        if (variance_model == "varying") paste0("var2_", indicator_names),
                        "group_logit", "profile_intercept", "profile_slope")
  stopifnot(length(parameter_names) == n_parameters)
  parameters <- setNames(values[seq_len(n_parameters)], parameter_names)
  mplus_means <- rbind(parameters[paste0("mu1_", indicator_names)],
                       parameters[paste0("mu2_", indicator_names)])
  mplus_variances <- rbind(parameters[paste0("var1_", indicator_names)],
    parameters[paste0(if (variance_model == "varying") "var2_" else "var1_", indicator_names)])
  profile_probability <- plogis(c(parameters["profile_intercept"] + parameters["profile_slope"],
                                  parameters["profile_intercept"]))
  mplus_profile <- cbind(profile_probability, 1 - profile_probability)
  mplus_group <- c(plogis(parameters["group_logit"]), 1 - plogis(parameters["group_logit"]))
  metrics <- setNames(tail(values, 7L),
                       c("n_parameters", "log_likelihood", "aic", "bic_individual",
                         "sabic", "joint_entropy", "condition_number"))
  saved <- read.table(file.path(artifact_dir, paste0(prefix, "-posteriors.dat")),
    col.names = c(indicator_names, "p11", "p12", "p21", "p22", "cb", "cw", "joint", "id", "clus"))
  saved <- saved[match(data$id, saved$id), , drop = FALSE]
  stopifnot(!anyNA(saved), !anyDuplicated(saved$id), identical(saved$id, data$id),
            identical(saved$clus, data$clus),
            max(abs(as.matrix(saved[indicator_names]) - as.matrix(data[indicator_names]))) < 1e-10)
  subject_posteriors <- cbind(saved$p11 + saved$p21, saved$p12 + saved$p22)
  group_by_subject <- cbind(saved$p11 + saved$p12, saved$p21 + saved$p22)
  first_group_rows <- match(unique(data$clus), data$clus)
  group_posteriors <- group_by_subject[first_group_rows, , drop = FALSE]
  stopifnot(max(abs(rowSums(subject_posteriors) - 1)) < 1e-10,
            max(abs(rowSums(group_posteriors) - 1)) < 1e-10,
            max(abs(group_by_subject - group_posteriors[match(data$clus, unique(data$clus)), ])) < 1e-10)
  # Only genuine Mplus values enter expected; native R fits never supply targets.
  reference <- list(data = data, indicators = indicator_names, variance_model = variance_model,
    expected = list(means = unname(mplus_means), variances = unname(mplus_variances),
      profile_probabilities = unname(mplus_profile), group_probabilities = unname(mplus_group),
      subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
      metrics = metrics),
    provenance = list(version = trimws(grep("^Mplus VERSION", output, value = TRUE)[1L]),
      data_url = "https://www.statmodel.com/usersguide/chap10/ex10.4.dat",
      original_model_url = "https://www.statmodel.com/usersguide/chap10/ex10.4.html",
      model = "New two-level diagonal Gaussian LPA; NOT original Example 10.4 CFA mixture",
      command = sprintf("mpdemo %s.inp %s.out", prefix, prefix),
      working_directory = artifact_dir,
      md5 = tools::md5sum(c(file.path(public_dir, "ex10.4.dat"),
        file.path(artifact_dir, c("ex104-lpa.dat", paste0(prefix, ".inp"), paste0(prefix, ".out"),
                                 paste0(prefix, "-results.dat"), paste0(prefix, "-posteriors.dat"))))),
      precision = "SAVEDATA RESULTS 8 significant digits; posterior FORMAT F20.12"))
  saveRDS(reference, file.path("tests", "fixtures", "mplus", paste0("twolevel-public-", variance_model, ".rds")))
  fit <- multilpa(data, indicator_names, "clus", 2L, 2L,
                    variance_model = variance_model, n_starts = 20L,
                    max_iter = 10000L, tol = 1e-14, seed = 20260917)
  profile_order <- order(fit$means[, "y1"])[rank(mplus_means[, 1L])]
  group_order <- order(fit$profile_probabilities[, profile_order[1L]])[rank(mplus_profile[, 1L])]
  differences <- c(
    log_likelihood = abs(fit$log_likelihood - unname(metrics["log_likelihood"])),
    aic = abs(fit$aic - unname(metrics["aic"])),
    bic_individual = abs(fit$bic_individual - unname(metrics["bic_individual"])),
    means = max(abs(fit$means[profile_order, ] - mplus_means)),
    variances = max(abs(fit$variances[profile_order, ] - mplus_variances)),
    profile_probabilities = max(abs(fit$profile_probabilities[group_order, profile_order] - mplus_profile)),
    group_probabilities = max(abs(fit$group_probabilities[group_order] - mplus_group)),
    subject_posteriors = max(abs(fit$subject_posteriors[, profile_order] - subject_posteriors)),
    group_posteriors = max(abs(fit$group_posteriors[, group_order] - group_posteriors)))
  print(variance_model)
  print(differences, digits = 10)
  stopifnot(fit$converged, !fit$boundary, fit$n_parameters == n_parameters,
    differences["log_likelihood"] <= 0.00005,
    all(differences[c("aic", "bic_individual")] <= 0.0005),
    all(differences[c("means", "variances", "profile_probabilities", "group_probabilities",
                      "subject_posteriors", "group_posteriors")] <= 0.0001))
  list(variance_model = variance_model, fit = fit, differences = differences,
        mplus_metrics = metrics, reference_path = paste0("twolevel-public-", variance_model, ".rds"))
})
names(comparisons) <- c("varying", "equal")
saveRDS(comparisons, file.path(artifact_dir, "comparison.rds"))
