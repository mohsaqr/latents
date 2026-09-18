# Rebuild offline targets from genuine Mplus output, then refit R independently.
pkgload::load_all(".")
directory <- file.path("validation", "mplus", "covariates")
data <- read.table(file.path(directory, "covariates.dat"),
                   col.names = c("y1", "y2", "z", "w", "g", "id"))
str(data)
print(head(data))
print(summary(data))
print(vapply(data, class, character(1)))
stopifnot(!anyNA(data), !anyDuplicated(data$id), nrow(data) == 720L)
values <- scan(file.path(directory, "covariates-results.dat"), quiet = TRUE)
stopifnot(length(values) == 33L, values[27] == 13L)
means <- rbind(values[1:2], values[5:6])
variances <- rbind(values[3:4], values[7:8])
beta <- matrix(c(values[10] + values[13], values[10], values[12]), 3, 1)
gamma <- matrix(values[c(9, 11)], 2, 1)
saved <- read.table(file.path(directory, "covariates-posteriors.dat"),
                    col.names = c("y1", "y2", "z", "w", "p11", "p12", "p21", "p22",
                                  "cb", "cw", "joint", "id", "g"))
stopifnot(nrow(saved) == nrow(data), !anyNA(saved), !anyDuplicated(saved$id))
saved <- saved[match(data$id, saved$id), ]
subject_posteriors <- cbind(saved$p11 + saved$p21, saved$p12 + saved$p22)
group_posteriors <- cbind(saved$p11 + saved$p12, saved$p21 + saved$p22)[!duplicated(data$g), ]
profile_order <- order(means[, 1])
if (profile_order[1] == 2L) beta <- -beta
means <- means[profile_order, , drop = FALSE]
variances <- variances[profile_order, , drop = FALSE]
subject_posteriors <- subject_posteriors[, profile_order]
group_order <- order(beta[1:2, 1], decreasing = TRUE)
beta[1:2, ] <- beta[group_order, ]
if (group_order[1] == 2L) gamma <- -gamma
group_posteriors <- group_posteriors[, group_order]
target <- list(means = means, variances = variances, profile_coefficients = beta,
               group_coefficients = gamma, subject_posteriors = subject_posteriors,
               group_posteriors = group_posteriors, log_likelihood = values[28],
               aic = values[29], bic_individual = values[30], n_parameters = 13L)
fixture <- list(data = data, expected = target, version = "Mplus VERSION 9 DEMO (Mac)",
                source = "New matching two-level ML analysis on synthetic data",
                md5 = tools::md5sum(file.path(directory, c("covariates.dat", "covariates.inp", "covariates.out",
                                                          "covariates-results.dat", "covariates-posteriors.dat"))))
saveRDS(fixture, file.path("tests", "fixtures", "mplus", "twolevel-covariates.rds"))
fit <- fit_covariates(data, c("y1", "y2"), "g", 2, 2, "z", "w",
                             n_starts = 10, seed = 812, tol = 1e-12)
p <- order(fit$means[, 1])
beta_r <- unname(fit$profile_coefficients)
if (p[1] == 2L) beta_r <- -beta_r
h <- order(beta_r[1:2, 1], decreasing = TRUE)
beta_r[1:2, ] <- beta_r[h, ]
gamma_r <- unname(fit$group_coefficients) * if (h[1] == 2L) -1 else 1
comparison <- c(log_likelihood = abs(fit$log_likelihood - target$log_likelihood),
 means = max(abs(fit$means[p, ] - target$means)),
 variances = max(abs(fit$variances[p, ] - target$variances)),
 profile_coefficients = max(abs(beta_r - target$profile_coefficients)),
 group_coefficients = max(abs(gamma_r - target$group_coefficients)),
 subject_posteriors = max(abs(fit$subject_posteriors[, p] - target$subject_posteriors)),
 group_posteriors = max(abs(fit$group_posteriors[, h] - target$group_posteriors)))
print(comparison, digits = 12)
stopifnot(fit$converged, comparison[1] < 5e-5, all(comparison[-1] < 1e-5))
saveRDS(list(fit = fit, expected = target, comparison = comparison), file.path(directory, "comparison.rds"))
