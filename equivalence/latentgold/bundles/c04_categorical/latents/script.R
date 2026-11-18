library(latents)
data <- utils::read.delim("../data/c04.dat", na.strings = ".")
fit <- multilpa(data, sprintf("c%d", 1:4), "g", n_profiles = 2, n_group_classes = 2, categorical = sprintf("c%d",
       1:4), n_starts = 50L, max_iter = 200000L, tol = 1e-14, seed = 1L)
print(fit)
summary(fit)
result <- data.frame(quantity = "log likelihood", value = as.numeric(logLik(fit)))
