pkgload::load_all(".", quiet = TRUE)
set.seed(21)
class <- sample(1:2, 300, TRUE)
data <- data.frame(a = stats::rnorm(300, c(-1.5, 1.5)[class]),
                   b = stats::rnorm(300, c(-1.5, 1.5)[class]))
fit <- function(d) enumerate_lpa(d, c("a", "b"), n_profiles = 1:2,
  model = "EEI", missing = "fiml", n_starts = 5, seed = 1, tol = 1e-12)
original <- get_results(fit(data), "lmr")
padded <- get_results(fit(rbind(data, data.frame(a = rep(NA_real_, 95),
  b = rep(NA_real_, 95)))), "lmr")
print(original)
print(padded)
cat("Differences:\n")
print(padded[c("n_units", "lr", "adjusted_lr", "vlmr_p", "almr_p")] -
        original[c("n_units", "lr", "adjusted_lr", "vlmr_p", "almr_p")])
