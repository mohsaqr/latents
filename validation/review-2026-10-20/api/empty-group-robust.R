library(latents, lib.loc = "validation/review-2026-10-20/rlib")
set.seed(17)
d <- data.frame(g = rep(1:5, each = 30), y = rnorm(150))
padded <- rbind(d, data.frame(g = 6:100, y = NA_real_))
fit <- function(data) suppressMessages(multilpa(data, "y", "g", n_profiles = 1,
  n_group_classes = 1, missing = "fiml", n_starts = 1, seed = 3))
a <- fit(d)
b <- fit(padded)
ta <- parameter_inference(a, vcov_type = "robust")
tb <- parameter_inference(b, vcov_type = "robust")
print(ta)
print(tb)
cat("Informative groups:", a$n_informative_groups, b$n_informative_groups, "\n")
cat("Log likelihood:", a$log_likelihood, b$log_likelihood, "\n")
ratios <- tb$standard_error[1:2] / ta$standard_error[1:2]
cat("SE ratios:", ratios, "\n")
cat("CR1 count ratio:", sqrt((100 / 99) / (5 / 4)), "\n")
stopifnot(abs(a$log_likelihood - b$log_likelihood) < 1e-9,
          max(abs(ta$estimate - tb$estimate)) < 1e-9,
          all(abs(ratios - sqrt((100 / 99) / (5 / 4))) < 1e-8))
