library(latents)
data <- utils::read.delim("../data/c10.dat", na.strings = ".")
fit <- lta(data, c("y1", "y2"), "g", n_profiles = 2, time = "t", n_group_classes = 2L, n_starts = 50L,
       max_iter = 200000L, tol = 1e-14, seed = 1L)
print(fit)
summary(fit)
result <- data.frame(quantity = "log likelihood", value = as.numeric(logLik(fit)))
