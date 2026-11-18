library(latents)
data <- utils::read.delim("../data/c01.dat", na.strings = ".")
fit <- multilpa(data, c("y1", "y2", "y3"), "g", n_profiles = 2, n_group_classes = 2, variance_model = "equal",
       n_starts = 50L, max_iter = 200000L, tol = 1e-14, seed = 1L)
print(fit)
summary(fit)
result <- data.frame(quantity = "log likelihood", value = as.numeric(logLik(fit)))
