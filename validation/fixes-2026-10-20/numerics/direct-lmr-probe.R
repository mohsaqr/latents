pkgload::load_all(quiet = TRUE)
set.seed(77)
n <- 240
class <- sample(1:2, n, TRUE)
z <- stats::rnorm(n, 1, 2)
data <- data.frame(z = z, a = stats::rnorm(n, c(-2, 2)[class]) + 0.7 * z,
                   b = stats::rnorm(n, c(-2, 2)[class]))
fit <- function(data, size = 1) suppressWarnings(enumerate_lpa(data, c("a", "b"),
  n_profiles = 1:2, model = "EEI", direct_effects = "a ~ z", n_starts = 3,
  seed = 1, tol = 1e-12, min_variance = 1e-6 * size^2))
base <- get_results(fit(data), "lmr")
shifted <- get_results(fit(transform(data, z = z + 1e6)), "lmr")
scaled <- get_results(fit(transform(data, a = a * 1e8, b = b * 1e8, z = z + 1e6),
                          1e8), "lmr")
print(base)
print(shifted)
print(scaled)
stopifnot(all(is.na(c(base$note, shifted$note, scaled$note))))
stopifnot(max(abs(base$reference_mean - c(shifted$reference_mean, scaled$reference_mean))) < 1e-5)
