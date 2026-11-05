suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# One-class mixture regressions against stats::glm, MASS, hand likelihoods, by family.
n <- 300
d <- data.frame(x = rnorm(n), z = rbinom(n, 1, .4), e = runif(n, 1, 3))
d$y <- 1 + .5 * d$x - .7 * d$z + rnorm(n, sd = 1.3)
d$cnt <- rpois(n, exp(.3 + .4 * d$x + log(d$e)))
d$nb <- rnbinom(n, size = 2, mu = exp(.5 + .4 * d$x))
d$s <- rbinom(n, 8, plogis(-.2 + .6 * d$x)); d$f <- 8 - d$s
d$bin <- rbinom(n, 1, plogis(.2 + .6 * d$x))
d$o <- cut(d$x * .8 + rlogis(n), c(-Inf, -1, .3, 1.5, Inf), labels = c("a", "b", "c", "d"), ordered_result = TRUE)
show <- function(label, a, b) cat(sprintf("%-38s latents %.8f  ref %.8f  diff %.2e\n", label, a, b, a - b))
cmp <- function(label, fit, ref_coef, ref_se, ref_ll, ref_k) {
  tb <- as.data.frame(fit)
  show(paste(label, "LL"), fit$log_likelihood, ref_ll)
  cat(sprintf("  k latents %d ref %d; max|coef diff| %.2e; max|se diff| %.2e  max rel se %.2e\n", fit$n_parameters, ref_k,
    max(abs(tb$estimate - ref_coef)), max(abs(tb$std_error - ref_se)), max(abs(tb$std_error / ref_se - 1))))
}
# Gaussian: lm's ML sigma; SE: sigma_ML^2 (X'X)^-1
f <- mixture_regression(y ~ x + z, d, 1, seed = 1, tol = 1e-12)
l <- lm(y ~ x + z, d); X <- model.matrix(l); s2 <- sum(resid(l)^2) / n
cmp("gaussian", f, coef(l), sqrt(diag(s2 * solve(crossprod(X)))), as.numeric(logLik(l)), 4L)
cat("  sigma", get_results(f, "classes")$sigma, sqrt(s2), " se(sigma)", get_results(f, "classes")$sigma_std_error, sqrt(s2 / (2 * n)), "\n")
ft <- get_results(f, "fit"); cat(sprintf("  BIC %.6f stats %.6f; AIC %.6f stats %.6f; sabic %.6f hand %.6f\n", ft$bic, BIC(l), ft$aic, AIC(l), ft$sabic, -2*logLik(l)+4*log((n+2)/24)))
# Poisson with offset
f <- mixture_regression(cnt ~ x + offset(log(e)), d, 1, family = "poisson", seed = 1, tol = 1e-12)
g <- glm(cnt ~ x + offset(log(e)), poisson, d)
cmp("poisson+offset", f, coef(g), sqrt(diag(vcov(g))), as.numeric(logLik(g)), 2L)
tb <- as.data.frame(f); cat("  exp_estimate", tb$exp_estimate, " exp(coef)", exp(coef(g)), "\n")
# Binomial counts and Bernoulli
f <- mixture_regression(cbind(s, f) ~ x, d, 1, family = "binomial", seed = 1, tol = 1e-12)
g <- glm(cbind(s, f) ~ x, binomial, d)
cmp("binomial counts", f, coef(g), sqrt(diag(vcov(g))), as.numeric(logLik(g)), 2L)
f <- mixture_regression(bin ~ x, d, 1, family = "binomial", seed = 1, tol = 1e-12)
g <- glm(bin ~ x, binomial, d)
cmp("bernoulli", f, coef(g), sqrt(diag(vcov(g))), as.numeric(logLik(g)), 2L)
# NB vs MASS::glm.nb
f <- mixture_regression(nb ~ x, d, 1, family = "negative_binomial", seed = 1, tol = 1e-12)
g <- MASS::glm.nb(nb ~ x, d, control = glm.control(epsilon = 1e-12, maxit = 200))
cmp("negbin (glm.nb SE conditional)", f, coef(g), sqrt(diag(vcov(g))), as.numeric(logLik(g)), 3L)
cl <- get_results(f, "classes"); cat(sprintf("  dispersion %.8f vs 1/theta %.8f; se %.6f vs SE.theta/theta^2 %.6f\n", cl$dispersion, 1 / g$theta, cl$dispersion_std_error, g$SE.theta / g$theta^2))
# hand NB ll
mu <- exp(cbind(1, d$x) %*% as.data.frame(f)$estimate); a <- cl$dispersion
show("negbin hand", f$log_likelihood, sum(lgamma(d$nb + 1/a) - lgamma(1/a) - lgamma(d$nb + 1) + (1/a) * log((1/a)/(1/a + mu)) + d$nb * log(mu/(1/a + mu))))
# Ordinal vs polr
f <- mixture_regression(o ~ x, d, 1, family = "ordinal", seed = 1, tol = 1e-12)
g <- MASS::polr(o ~ x, d, Hess = TRUE, control = list(reltol = 1e-14))
tb <- as.data.frame(f); print(tb[c("term", "estimate", "std_error", "exp_estimate", "p_adjusted")])
print(summary(g)$coefficients[, 1:2])
cat(sprintf("  LL %.8f polr %.8f  k %d\n", f$log_likelihood, as.numeric(logLik(g)), f$n_parameters))
# Gaussian no-intercept and common with factor
d$g <- factor(sample(c("u", "v", "w"), n, TRUE))
f <- mixture_regression(y ~ 0 + g + x, d, 1, seed = 1, tol = 1e-12); l <- lm(y ~ 0 + g + x, d)
cat("no intercept diff", max(abs(as.data.frame(f)$estimate - coef(l))), "\n")
print(as.data.frame(f)[c("term", "p_value", "p_adjusted")]); print(p.adjust(summary(l)$coefficients[, 4], "BH"))
