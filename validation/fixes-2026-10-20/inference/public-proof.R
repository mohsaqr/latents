pkgload::load_all(quiet = TRUE)
source("tests/testthat/helper-expected-warnings.R")
set.seed(792)
fit <- function(data) quietly(multilpa(data, "y", "g", 1, 1,
  missing = "fiml", n_starts = 1, seed = 3, tol = 1e-12))
data <- data.frame(g = 1, y = stats::rnorm(50))
padded <- rbind(data, data.frame(g = 2:96, y = NA_real_))
cat("One informative group plus 95 empty groups:\n")
invisible(lapply(c("robust", "opg"), function(type) {
  answer <- tryCatch(parameter_inference(fit(padded), vcov_type = type), error = identity)
  cat(type, paste(class(answer), collapse = "/"), conditionMessage(answer), "\n")
  stopifnot(inherits(answer, "latents_too_few_groups"))
}))
data <- data.frame(g = rep(1:5, each = 50), y = stats::rnorm(250))
padded <- rbind(data, data.frame(g = 6:100, y = NA_real_))
a <- parameter_inference(fit(data), vcov_type = "robust")
b <- parameter_inference(fit(padded), vcov_type = "robust")
cat("Five informative groups, before/after 95 empty groups:\n")
print(data.frame(estimate = c(a$estimate[1], b$estimate[1]),
                 standard_error = c(a$standard_error[1], b$standard_error[1])))
stopifnot(isTRUE(all.equal(attr(a, "covariance"), attr(b, "covariance"), tolerance = 1e-9)))

set.seed(802)
data <- data.frame(g = rep(1:20, each = 20), z = stats::rnorm(400), w = stats::rnorm(400))
data$y <- 0.8 + 0.7 * data$z - 0.4 * data$w + stats::rnorm(400)
direct_fit <- function(data) quietly(multilpa(data, "y", "g", 1, 1,
  direct_effects = "y ~ z + w", n_starts = 1, seed = 3, tol = 1e-12))
base <- direct_fit(data)
shifted_data <- transform(data, z = z + 2026, w = w - 1000)
shifted <- direct_fit(shifted_data)
cat("Direct effects after predictor translation z + 2026, w - 1000:\n")
cat("Log likelihood:", format(c(base$log_likelihood, shifted$log_likelihood), digits = 15), "\n")
invisible(lapply(c("observed", "robust", "opg"), function(type) {
  a <- parameter_inference(base, vcov_type = type)
  b <- parameter_inference(shifted, vcov_type = type)
  cat(type, "\n")
  print(data.frame(parameter = a$term[c(1, 3, 4)], estimate = b$estimate[c(1, 3, 4)],
                   original_se = a$standard_error[c(1, 3, 4)], shifted_se = b$standard_error[c(1, 3, 4)]))
}))
reference <- stats::lm(y ~ z + w, shifted_data)
reference_covariance <- summary(reference)$cov.unscaled * mean(stats::residuals(reference)^2)
covariance <- attr(parameter_inference(shifted), "covariance")[c(1, 3, 4), c(1, 3, 4)]
cat("Maximum relative lm normal-MLE covariance discrepancy:",
    max(abs(covariance - reference_covariance)) / max(abs(reference_covariance)), "\n")
stopifnot(isTRUE(all.equal(unname(covariance), unname(reference_covariance), tolerance = 2e-7)))
