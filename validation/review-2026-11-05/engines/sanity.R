pkgload::load_all(".", quiet = TRUE)
options(width = 110)
set.seed(110507)
n <- 240
state <- sample(1:2, n, TRUE)
dd <- data.frame(a = rnorm(n, c(-2, 2)[state], .7),
                 b = rnorm(n, c(1, -1)[state], .9))
fit <- lpa(dd, c("a", "b"), 2, model = "VVI", seed = 46, n_starts = 3, tol = 1e-10)
print(class(fit)); print(names(fit)); print(head(fit$subject_posterior))
weighted_densities <- vapply(seq_len(2), function(k) {
  fit$profile_probabilities[1, k] * dnorm(dd$a, fit$means[k, "a"], sqrt(fit$variances[k, "a"])) *
    dnorm(dd$b, fit$means[k, "b"], sqrt(fit$variances[k, "b"]))
}, numeric(n))
ll_hand <- sum(log(rowSums(weighted_densities)))
cat("VVI hand likelihood error:", ll_hand - fit$log_likelihood, "\n")
tiny <- transform(dd, a = a * 1e-6, b = b * 1e-6)
fit_tiny <- lpa(tiny, c("a", "b"), 2, model = "VVI", seed = 46, n_starts = 3, tol = 1e-10)
jacobian_error <- fit_tiny$log_likelihood - fit$log_likelihood + 2 * n * log(1e-6)
cat("Tiny-unit likelihood Jacobian error:", jacobian_error, "\n")
reg <- transform(dd, y = 1 + .7 * a + rnorm(n, sd = .8))
fr <- mixture_regression(y ~ a, reg, 1, seed = 22)
reference <- lm(y ~ a, reg)
cat("Gaussian regression LL error vs lm ML:", fr$log_likelihood - as.numeric(logLik(reference)), "\n")
print(coef(fr)); print(coef(reference))
stopifnot(abs(ll_hand - fit$log_likelihood) < 1e-8,
          abs(jacobian_error) < 1e-5,
          abs(fr$log_likelihood - as.numeric(logLik(reference))) < 1e-8)
