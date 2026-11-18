library(latents)
data <- utils::read.delim("../data/c12.dat", na.strings = ".")
fit <- multilpa(data, c("browse", "lectures", "forum_read", "forum_post", "attendance"),
       "g", n_profiles = 2, n_group_classes = 2, n_starts = 50L, max_iter = 200000L,
       tol = 1e-14, seed = 1L)
print(fit)
summary(fit)
result <- data.frame(quantity = "log likelihood", value = as.numeric(logLik(fit)))
