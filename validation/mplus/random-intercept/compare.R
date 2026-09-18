# Run from the package root after the retained Mplus analyses.
pkgload::load_all(".")
directory <- file.path("validation", "mplus", "random-intercept")
data <- read.table(file.path(directory, "one-profile.dat"), col.names = c("y1", "y2", "group"))
str(data)
print(head(data))
print(summary(data))
stopifnot(!anyNA(data), all(vapply(data, is.numeric, logical(1))))
fit <- fit_multilpa_random_intercept(data, "y1", "group", 1L,
  n_starts = 2L, tol = 1e-11, seed = 983L)
reference <- scan(file.path(directory, "univariate-results.dat"), quiet = TRUE)
stopifnot(length(reference) == 21L, reference[7L] == 3L,
  any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", readLines(file.path(directory, "univariate.out")), fixed = TRUE)))
comparison <- data.frame(quantity = c("within_variance", "mean", "random_variance", "log_likelihood", "aic", "bic_individual"),
  native = c(fit$variances, fit$means, fit$random_sd^2, fit$log_likelihood, fit$aic, fit$bic_individual),
  mplus = reference[c(1L, 2L, 3L, 8L, 10L, 11L)])
comparison$absolute_difference <- abs(comparison$native - comparison$mplus)
print(comparison, digits = 12)
stopifnot(max(comparison$absolute_difference[seq_len(3L)]) < 1e-5,
  comparison$absolute_difference[4L] < 5e-5,
  max(comparison$absolute_difference[5:6]) < 1e-4)
write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
saveRDS(list(data = data[c("y1", "group")],
  mplus = list(parameters = reference[seq_len(3L)], n_parameters = reference[7L],
    log_likelihood = reference[8L], aic = reference[10L], bic = reference[11L]),
  mplus_version = "Mplus VERSION 9 DEMO (Mac)",
  source = file.path(directory, "univariate.out"),
  md5 = tools::md5sum(file.path(directory, c("one-profile.dat", "univariate.out", "univariate-results.dat")))),
  file.path("tests", "fixtures", "mplus", "random-intercept-one-profile.rds"))

# Diagnose the separate singular two-indicator Mplus specification without
# pretending it is the exact rank-one random-intercept model.
two_reference <- scan(file.path(directory, "one-profile-mixture-results.dat"), quiet = TRUE)
dense_log_likelihood <- function(between_ridge) {
  stopifnot(is.numeric(between_ridge), length(between_ridge) == 1L, between_ridge >= 0)
  sum(vapply(split(seq_len(nrow(data)), data$group), function(rows) {
    residual <- as.vector(sweep(as.matrix(data[rows, c("y1", "y2")]), 2L, two_reference[3:4], "-"))
    n <- length(rows)
    covariance <- diag(rep(two_reference[1:2], each = n)) + two_reference[5L] +
      kronecker(diag(2L) * between_ridge, matrix(1, n, n))
    -0.5 * (length(residual) * log(2 * pi) + as.numeric(determinant(covariance)$modulus) +
      sum(residual * solve(covariance, residual)))
  }, numeric(1)))
}
diagnostic <- c(exact_rank_one = dense_log_likelihood(0),
  between_ridge_1e_4 = dense_log_likelihood(1e-4), mplus = two_reference[12L])
print(diagnostic, digits = 14)
stopifnot(abs(diagnostic[2L] - diagnostic[3L]) < 5e-5,
  abs(diagnostic[1L] - diagnostic[3L]) > 1e-3)
