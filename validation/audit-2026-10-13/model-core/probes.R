Sys.setenv(NOT_CRAN = 'true')
pkgload::load_all('.', quiet = TRUE)
ns <- asNamespace('latents')
f <- function(name) get(name, ns)
results <- list()
check_growth <- function(amplitude, noise = c(-1, 0, 1), sigma2 = 1) {
  y <- amplitude + noise
  n <- length(y)
  spec <- list(x = matrix(1, n, 1), z = matrix(numeric(), n, 0),
               random_design = matrix(1, n, 1), y = y, offset = rep(0, n),
               n = n, group_index = rep(1L, n))
  statistics <- f('.growth_statistics')(spec)
  g <- amplitude^2
  got <- f('.growth_class_terms')(statistics, 0, matrix(sqrt(g), 1, 1), sigma2)
  quadratic <- sum((y - mean(y))^2) / sigma2 + n * mean(y)^2 / (sigma2 + n * g)
  reference <- -0.5 * (n * log(2*pi) + (n-1) * log(sigma2) + log(sigma2+n*g) + quadratic)
  data.frame(amplitude, actual = got$log_density, reference, difference = got$log_density-reference,
             inferred_quadratic = -2*got$log_density - n*log(2*pi) - (n-1)*log(sigma2) - log(sigma2+n*g), reference_quadratic = quadratic)
}
results$growth_cancellation <- do.call(rbind, lapply(c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8), check_growth))
print(results$growth_cancellation)
# Overflow-resistant sampling weights, including exact zero weights.
results$weights <- f('.latents_sampling_weights')(data.frame(w=c(0,1e308,1e308)), 'w', 1:3, 3)
print(results$weights)
# Invalid counts and designs should receive classed input refusals.
run_condition <- function(expr) tryCatch({force(expr); 'accepted'}, error = function(e) paste(class(e), collapse=';'))
set.seed(100)
d <- data.frame(y=rnorm(50), x=rnorm(50), id=rep(1:10, each=5))
results$invalid <- list(nonfinite_response=run_condition(mixture_regression(y~x, transform(d,y=Inf),1,vcov_type='none')),
                       nonfinite_design=run_condition(mixture_regression(y~I(exp(x)),transform(d,x=1e308),1,vcov_type='none')))
print(results$invalid)
saveRDS(results, 'validation/audit-2026-10-13/model-core/probes.rds')
