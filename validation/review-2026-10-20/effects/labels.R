pkgload::load_all(quiet = TRUE)
attempt <- function(code) tryCatch(force(code), error = function(e) e)

# Ordinary categorical encoding protects distinct numeric values whose printed
# labels coincide. Joining a pair first bypasses that protection.
d <- data.frame(a = rep(c(.3, .1 + .2), 100),
                b = rep(c(0, 0, 1, 1), 50), c = rep(0:1, each = 100))
plain <- attempt(lca(d, names(d), 1, n_starts = 1))
paired <- lca(d, names(d), 1, local_dependence = "a:b", n_starts = 1)
reference_data <- transform(d, a = match(a, unique(a)))
reference <- lca(reference_data, names(reference_data), 1,
                 local_dependence = "a:b", n_starts = 1)
cat("Distinct native values:", length(unique(d$a)), "\n")
cat("Unpaired fit:", conditionMessage(plain), "\n")
print(data.frame(implementation = c("current pair", "distinct numeric categories"),
  parameters = c(paired$n_parameters, reference$n_parameters),
  likelihood = c(paired$log_likelihood, reference$log_likelihood)))
print(paired$categorical_levels)
stopifnot(inherits(plain, "latents_bad_categorical"), paired$n_parameters == 2,
          reference$n_parameters == 4)
