library(latents)
data <- utils::read.delim("../data/c05.dat", na.strings = ".")
fit <- multilpa(data, c("y1", "y2", "c1", "c2"), "g", n_profiles = 2, n_group_classes = 2,
       categorical = c("c1", "c2"), n_starts = 50L, max_iter = 200000L, tol = 1e-14,
       seed = 1L)
print(fit)
summary(fit)
result <- data.frame(quantity = "log likelihood", value = as.numeric(logLik(fit)))
