pkgload::load_all(quiet = TRUE)
set.seed(551)
n <- 400
d <- data.frame(z = rnorm(n), w = rnorm(n))
d$y1 <- 1 + .4 * d$z + .7 * d$w + rnorm(n)
d$y2 <- rnorm(n)
attempt <- function(code) tryCatch(suppressWarnings(force(code)), error = function(e) e)
show <- function(label, x) {
  cat("\n", label, "\n", sep = "")
  if (inherits(x, "error")) { cat(class(x), conditionMessage(x), "\n") } else {
    print(x$direct_effects)
    print(c(log_likelihood = x$log_likelihood, converged = x$converged))
  }
}
base <- attempt(lpa(d, c("y1", "y2"), 1, direct_effects = "y1 ~ z + w", n_starts = 1, seed = 1, tol = 1e-12))
show("base", base)
large <- attempt(lpa(transform(d, w = w * 1e4), c("y1", "y2"), 1, direct_effects = "y1 ~ z + w", n_starts = 1, seed = 1, tol = 1e-12))
show("one covariate in large units", large)
small <- attempt(lpa(transform(d, w = w * 1e-4), c("y1", "y2"), 1, direct_effects = "y1 ~ z + w", n_starts = 1, seed = 1, tol = 1e-12))
show("one covariate in small units", small)
shifted <- attempt(lpa(transform(d, z = z + 1e5), c("y1", "y2"), 1, direct_effects = "y1 ~ z + w", n_starts = 1, seed = 1, tol = 1e-12))
show("shifted", shifted)
cat("\nBase inference\n")
print(attempt(parameter_inference(base)))
cat("\nShifted inference\n")
print(attempt(parameter_inference(shifted)))

set.seed(22)
items <- as.data.frame(matrix(rbinom(400 * 4, 1, .5), 400, 4))
names(items) <- c("a", "b", "a:b", "c")
plain <- attempt(lca(items, names(items), 1, n_starts = 1))
collision <- attempt(lca(items, names(items), 1, local_dependence = "a:b", n_starts = 1))
cat("\nName collision plain parameter count, fit parameter count, original vars and fitted vars\n")
print(plain$n_parameters)
if (inherits(collision, "error")) print(collision) else {
  print(collision$n_parameters)
  print(collision$original_vars)
  print(collision$vars)
}
