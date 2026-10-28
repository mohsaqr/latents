pkgload::load_all(quiet = TRUE)
set.seed(551)
n <- 400
d <- data.frame(z = rnorm(n), w = rnorm(n))
d$y1 <- 1 + .4 * d$z + .7 * d$w + rnorm(n)
d$y2 <- rnorm(n)
attempt <- function(code) tryCatch(force(code), error = function(e) e)
fit <- function(data) suppressWarnings(lpa(data, c("y1", "y2"), 1,
  direct_effects = "y1 ~ z + w", n_starts = 1, seed = 1, tol = 1e-12))
base <- fit(d)
base_inf <- parameter_inference(base)
cat("Covariate shifts, unchanged maximum likelihood, and inference result:\n")
shifts <- c(0, 10, 100, 300, 1000, 2026)
shift_results <- do.call(rbind, lapply(shifts, function(shift) {
  changed <- transform(d, z = z + shift)
  fitted <- fit(changed)
  inference <- attempt(parameter_inference(fitted))
  reference <- lm(y1 ~ z + w, data = changed)
  reference_se <- sqrt(diag(vcov(reference)) * reference$df.residual / n)
  data.frame(shift = shift, likelihood_delta = fitted$log_likelihood - base$log_likelihood,
    effect_delta = max(abs(fitted$direct_effects - base$direct_effects)),
    slope_reference_se = reference_se[["z"]],
    reference_intercept_se = reference_se[["(Intercept)"]],
    inference = if (inherits(inference, "error")) class(inference)[1L] else "success")
}))
print(shift_results, digits = 12)
stopifnot(all(abs(shift_results$likelihood_delta) < 1e-8),
          shift_results$inference[1L] == "success",
          shift_results$inference[length(shifts)] == "latents_singular_information")
cat("Base direct-effect SEs versus independent lm MLE reference:\n")
print(subset(base_inf, parameter == "direct_effect", select = c(term, standard_error)))
ref <- lm(y1 ~ z + w, data = d)
print(sqrt(diag(vcov(ref)) * ref$df.residual / n))

cat("\nMultiple-covariate units: fit status, likelihood and coefficients:\n")
scales <- c(1, 1e4, 1e6, 1e7, 1e8)
unit_results <- do.call(rbind, lapply(scales, function(scale) {
  fitted <- attempt(fit(transform(d, w = w * scale)))
  data.frame(scale = scale,
    likelihood = if (inherits(fitted, "error")) NA_real_ else fitted$log_likelihood,
    result = if (inherits(fitted, "error")) conditionMessage(fitted) else "success")
}))
print(unit_results, digits = 12)

set.seed(22)
items <- as.data.frame(matrix(rbinom(400 * 4, 1, .5), 400, 4))
names(items) <- c("a", "b", "a:b", "c")
colliding <- lca(items, names(items), 1, local_dependence = "a:b", n_starts = 1)
manual <- items
manual$joint <- interaction(manual$a, manual$b, sep = ":", lex.order = TRUE)
correct <- lca(manual, c("joint", "a:b", "c"), 1, n_starts = 1)
cat("\nExisting indicator a:b overwritten by generated a:b:\n")
print(data.frame(implementation = c("current local_dependence", "collision-free joint reference"),
  n_parameters = c(colliding$n_parameters, correct$n_parameters),
  log_likelihood = c(colliding$log_likelihood, correct$log_likelihood)))
stopifnot(colliding$n_parameters == 4, correct$n_parameters == 5,
          !identical(as.character(items[["a:b"]]),
                     as.character(get_results(colliding, "data")[["a:b"]])))

cat("\nAdditional combination smoke checks:\n")
set.seed(111)
mixed <- transform(d, a = factor(rbinom(n, 1, .5)), b = factor(rbinom(n, 1, .5)),
                   count = rpois(n, 2), ordinal = ordered(sample(1:3, n, TRUE)))
extended <- lpa(mixed, c("y1", "y2", "a", "b", "count", "ordinal"), 1,
  categorical = c("a", "b"), count = "count", ordinal = "ordinal",
  local_dependence = "a:b", direct_effects = "y1 ~ z + w", n_starts = 1, seed = 1,
  tol = 1e-12)
stopifnot(isTRUE(extended$converged))
stopifnot(isTRUE(all.equal(predict(extended, type = "posterior"),
                         predict(extended, mixed, type = "posterior"))))
extended_inf <- parameter_inference(extended)
print(subset(extended_inf, parameter == "direct_effect"))
cat("Combination fit, original-data prediction and direct-effect inference passed.\n")
