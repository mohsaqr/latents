Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(quiet = TRUE)
capture <- function(expr) tryCatch(expr, error = function(e) list(error = conditionMessage(e), class = class(e)))
set.seed(210)
n <- 120L
d <- data.frame(g = seq_len(n), a = rnorm(n), b = rnorm(n), w = rep(c(1, 10), each = n / 2))
d$b[seq_len(n / 2)] <- d$a[seq_len(n / 2)] + rnorm(n / 2, sd = 0.1)
d$b[n / 2 + seq_len(n / 2)] <- -d$a[n / 2 + seq_len(n / 2)] + rnorm(n / 2, sd = 0.1)
fit <- multilpa(d, c("a", "b"), "g", n_profiles = 1, n_group_classes = 1,
                weights = "w", n_starts = 1, seed = 1)
cat("WEIGHTED RESIDUALS\n")
actual <- get_results(fit, "residuals")
print(actual)
r <- sweep(as.matrix(d[c("a", "b")]), 2, fit$means[1, ], "-")
scatter <- crossprod(r * d$w, r)
weighted_correlation <- scatter[1, 2] / sqrt(prod(diag(scatter)))
cat("Independent weighted residual correlation:", weighted_correlation, "\n")
cat("Converged:", fit$converged, "sampling weights:", range(fit$sampling_weights), "\n")
cat("ONE CLASS IDENTICAL IMPUTATIONS\n")
pooled <- capture(pool_imputations(list(d, d), c("a", "b"), "g", n_profiles = 1,
                                   n_group_classes = 1, n_starts = 1, seed = 1))
if (inherits(pooled, "latents_pooled")) {
  print(get_results(pooled))
  print(summary(pooled))
} else print(pooled)
cat("RUBIN ZERO WITHIN VARIANCE\n")
tab <- data.frame(level = "group", outcome = "group_1", term = "(Intercept)", parameter = "probability",
                  estimate = 1, standard_error = 0, statistic = NA_real_)
print(latents:::.multilpa_rubin(list(tab, tab), 0.95))
cat("DESCRIPTIVES MATRIX COLUMN\n")
mat <- data.frame(g = rep(1:3, each = 2)); mat$a <- matrix(1:12, ncol = 2)
print(capture(descriptives(mat, vars = "a", id = "g")))
cat("BOOTSTRAP DATA OVERRIDE\n")
set.seed(719)
original <- data.frame(g = rep(seq_len(30), each = 4), a = rnorm(120), b = rnorm(120))
one <- multilpa(original, c("a", "b"), "g", n_profiles = 1,
               n_group_classes = 1, n_starts = 1, seed = 1)
wrong <- original
wrong$a <- wrong$a * 10
wrong$b <- wrong$b * 10
wald_wrong <- capture(parameter_inference(one, data = wrong))
print(wald_wrong)
correct_boot <- capture(parameter_inference(one, method = "bootstrap", iter = 12,
                                            n_starts = 1, seed = 1))
wrong_boot <- capture(parameter_inference(one, data = wrong, method = "bootstrap", iter = 12,
                                          n_starts = 1, seed = 1))
if (is.data.frame(correct_boot) && is.data.frame(wrong_boot)) {
  print(data.frame(parameter = correct_boot$parameter,
                   original_point = correct_boot$estimate,
                   wrong_point = wrong_boot$estimate,
                   correct_se = correct_boot$standard_error,
                   wrong_se = wrong_boot$standard_error,
                   se_ratio = wrong_boot$standard_error / correct_boot$standard_error,
                   wrong_ci_low = wrong_boot$conf_low, wrong_ci_high = wrong_boot$conf_high))
} else {
  print(correct_boot)
  print(wrong_boot)
}
cat("R3STEP PREDICTOR UNITS\n")
set.seed(606)
n <- 400L
membership <- rnorm(n)
state <- rbinom(n, 1, plogis(1.5 * membership))
steps <- data.frame(g = seq_len(n), a = rnorm(n, 6 * state, 0.8),
                    b = rnorm(n, 5 * state, 0.8), z = membership)
first_stage <- multilpa(steps, c("a", "b"), "g", 2, 1, n_starts = 1, seed = 1)
standard <- capture(r3step(first_stage, steps, "z", vcov_type = "observed"))
scaled <- steps
scaled$z <- scaled$z * 1e-9
tiny <- capture(r3step(first_stage, scaled, "z", vcov_type = "observed"))
print(standard)
print(tiny)
