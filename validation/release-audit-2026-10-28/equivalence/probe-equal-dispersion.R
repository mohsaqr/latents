pkgload::load_all('.', quiet = TRUE)
set.seed(2810)
n <- 400L
x <- rnorm(n)
d <- data.frame(x = x, y = rnbinom(n, mu = exp(0.5 + 0.3 * x), size = 3))
fit <- mixture_regression(y ~ x, d, n_classes = 1L, family = 'negative_binomial', variance = 'equal', n_starts = 1L, seed = 1L)
str(d)
print(fit$converged)
print(fit$params$sigma2)
print(coef(fit))
print(as.data.frame(get_results(fit, 'classes')))
stopifnot(fit$converged, all(is.na(get_results(fit, 'classes')$dispersion)), all(is.na(get_results(fit, 'classes')$dispersion_std_error)))
# Counterfactual: drop incidental input names only inside the packer's runtime wrapper.
original <- latents:::.mixture_pack
fixed <- function(spec, params) {
  params$sigma2 <- unname(params$sigma2)
  original(spec, params)
}
assignInNamespace('.mixture_pack', fixed, ns = 'latents')
corrected <- mixture_regression(y ~ x, d, n_classes = 1L, family = 'negative_binomial', variance = 'equal', n_starts = 1L, seed = 1L)
print(coef(corrected))
print(as.data.frame(get_results(corrected, 'classes')))
stopifnot(isTRUE(all.equal(unname(coef(fit)), unname(coef(corrected)), tolerance = 1e-12)), all(is.finite(get_results(corrected, 'classes')$dispersion)))
