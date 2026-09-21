# Compare native MLR robust standard errors, the MLR scaling correction factor,
# and the sample-size adjusted BIC against a retained genuine Mplus 9 run.
# Run from the project root: Rscript validation/mplus/robust/compare.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
artifact_dir <- file.path("validation", "mplus", "robust")

output <- readLines(file.path(artifact_dir, "robust.out"))
stopifnot(
  "Mplus did not terminate normally" =
    any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", output)),
  "Mplus did not replicate the best loglikelihood" =
    any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", output)),
  "Mplus reported an error" =
    !any(grepl("*** ERROR", output, fixed = TRUE)),
  "Mplus did not use the MLR estimator" =
    any(grepl("H0 Scaling Correction Factor", output)))

values <- scan(file.path(artifact_dir, "robust-results.dat"), quiet = TRUE)
n_parameters <- 11L
stopifnot("unexpected saved-results layout" =
            length(values) == 2L * n_parameters + 8L &&
            values[2L * n_parameters + 1L] == n_parameters)
mplus_standard_errors <- values[n_parameters + seq_len(n_parameters)]
mplus_log_likelihood <- values[2L * n_parameters + 2L]
mplus_scaling <- values[2L * n_parameters + 3L]
mplus_aic <- values[2L * n_parameters + 4L]
mplus_bic <- values[2L * n_parameters + 5L]
mplus_sabic <- values[2L * n_parameters + 6L]

reference <- readRDS(file.path("tests", "fixtures", "mplus",
                               "twolevel-synthetic-varying.rds"))
dat <- reference$data
str(dat)
print(head(dat))
print(vapply(dat, class, character(1)))

fit <- multilpa(dat, c("y1", "y2"), "clus", 2, 2, variance_model = "varying",
                  start = starting_values(reference), n_starts = 1,
                  max_iter = 5000, tol = 1e-13)
information <- parameter_inference(fit, dat, vcov_type = "robust")

# Mplus orders the measurement block by profile then indicator, reports
# variances rather than log variances, and parameterizes the mixing
# probabilities as a group logit, a within-class intercept at cb = 2, and the
# cb = 1 minus cb = 2 contrast.
q <- fit$n_parameters
n_measurement <- q - 3L
ordering <- c(1, 2, 5, 6, 3, 4, 7, 8)
measurement_derivatives <- c(rep(1, 4), as.vector(t(fit$variances)))
jacobian <- matrix(0, q, q)
jacobian[cbind(seq_len(n_measurement), ordering)] <- measurement_derivatives[ordering]
jacobian[q - 2L, q] <- 1
jacobian[q - 1L, q - 1L] <- 1
jacobian[q, c(q - 2L, q - 1L)] <- c(1, -1)
native_standard_errors <- sqrt(diag(jacobian %*% attr(information, "covariance_unconstrained") %*%
                                      t(jacobian)))

# `format = "long"` since 0.9.0: the default became one wide row, and the
# `convention` of a criterion that has none -- aic, kic, deviance -- became
# NA_character_ rather than "none". Both are matched explicitly here.
indices <- get_data(fit, "information_criteria", format = "long")
native <- function(criterion, convention) {
  matched <- indices$criterion == criterion &
    if (is.na(convention)) is.na(indices$convention) else
      !is.na(indices$convention) & indices$convention == convention
  stopifnot("exactly one criterion row must match" = sum(matched) == 1L)
  indices$value[matched]
}

differences <- c(
  log_likelihood = abs(fit$log_likelihood - mplus_log_likelihood),
  robust_standard_errors = max(abs(native_standard_errors - mplus_standard_errors)),
  scaling_correction = abs(attr(information, "scaling_correction") - mplus_scaling),
  aic = abs(native("aic", NA_character_) - mplus_aic),
  bic_individual = abs(native("bic", "individuals") - mplus_bic),
  sabic_individual = abs(native("sabic", "individuals") - mplus_sabic))
print(differences, digits = 12)

comparison <- data.frame(
  quantity = c(sprintf("robust_se_%d", seq_len(q)), "scaling_correction",
               "aic", "bic_individual", "sabic_individual"),
  native = c(native_standard_errors, attr(information, "scaling_correction"),
             native("aic", NA_character_), native("bic", "individuals"),
             native("sabic", "individuals")),
  mplus = c(mplus_standard_errors, mplus_scaling, mplus_aic, mplus_bic, mplus_sabic))
comparison$absolute_difference <- abs(comparison$native - comparison$mplus)
print(comparison, digits = 12, row.names = FALSE)

# Saved Mplus results carry eight significant digits, and the printed criteria
# three decimals; tolerances below account for that finite precision.
stopifnot(
  "robust standard errors disagree with Mplus" =
    differences[["robust_standard_errors"]] < 1e-6,
  "MLR scaling correction disagrees with Mplus" =
    differences[["scaling_correction"]] < 1e-6,
  "sample-size adjusted BIC disagrees with Mplus" =
    differences[["sabic_individual"]] < 1e-3)

saveRDS(list(data = dat, mplus_standard_errors = mplus_standard_errors,
             mplus_scaling = mplus_scaling, mplus_aic = mplus_aic,
             mplus_bic = mplus_bic, mplus_sabic = mplus_sabic,
             means = reference$means, variances = reference$variances,
             profile_probabilities = reference$profile_probabilities,
             group_probabilities = reference$group_probabilities),
        file.path("tests", "fixtures", "mplus", "twolevel-robust.rds"))
write.csv(comparison, file.path(artifact_dir, "comparison.csv"), row.names = FALSE)
cat("Robust comparison complete.\n")
