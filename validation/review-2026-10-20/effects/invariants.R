pkgload::load_all(quiet = TRUE)

# Independent one-profile Gaussian reference for a nonadjacent, disjoint
# covariance pattern. Covariance uses MLE divisor n, not stats::cov's n - 1.
set.seed(615)
n <- 600
raw <- matrix(rnorm(n * 4), n, 4)
raw[, 4] <- .6 * raw[, 1] + .8 * raw[, 4]
raw[, 3] <- -.4 * raw[, 2] + sqrt(.84) * raw[, 3]
d <- as.data.frame(raw)
names(d) <- paste0("y", seq_len(4))
f <- lpa(d, names(d), 1, local_dependence = c("y1:y4", "y2:y3"),
         n_starts = 1, seed = 1, tol = 1e-12)
expected <- cov(raw) * (n - 1) / n
expected[!f$covariance_pattern] <- 0
stopifnot(max(abs(f$means[1, ] - colMeans(raw))) < 1e-10,
          max(abs(f$covariances[, , 1] - expected)) < 1e-10,
          f$n_parameters == 10)
print(c(pattern_max_covariance_error = max(abs(f$covariances[, , 1] - expected))))
print(subset(parameter_inference(f), parameter == "covariance"))

# Missing Gaussian responses should reduce to a complete-case regression for
# the affected indicator at one profile, with no contribution by missing y1.
set.seed(616)
d <- data.frame(z = rnorm(n), w = rnorm(n))
d$y1 <- 1 + .4 * d$z + .7 * d$w + rnorm(n)
d$y2 <- rnorm(n)
d$y1[sample.int(n, 120)] <- NA_real_
f <- lpa(d, c("y1", "y2"), 1, direct_effects = "y1 ~ z + w",
         missing = "fiml", n_starts = 1, seed = 1, tol = 1e-12)
reference <- lm(y1 ~ z + w, data = d)
coefficients <- c(f$means[1, "y1"], unname(f$direct_effects))
print(c(missing_effect_max_error = max(abs(coefficients - coef(reference)))))
stopifnot(max(abs(coefficients - coef(reference))) < 1e-5,
          abs(f$variances[1, "y1"] - sum(resid(reference)^2) / nobs(reference)) < 1e-5)
print(subset(parameter_inference(f), parameter == "direct_effect"))

# Sampling weights affect the direct-effect least-squares update too.
set.seed(617)
d$weight <- runif(n, .2, 2)
f <- lpa(d, c("y1", "y2"), 1, direct_effects = "y1 ~ z + w", weights = "weight",
         missing = "fiml", n_starts = 1, seed = 1, tol = 1e-12)
reference <- lm(y1 ~ z + w, data = d, weights = weight)
coefficients <- c(f$means[1, "y1"], unname(f$direct_effects))
print(c(weighted_missing_effect_max_error = max(abs(coefficients - coef(reference)))))
stopifnot(max(abs(coefficients - coef(reference))) < 1e-5)
cat("Independent covariance, missing-response and weighted regression references passed.\n")
