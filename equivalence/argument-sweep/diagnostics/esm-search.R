# Does a longer search reach the optimum poLCA found on the incomplete ESM items?
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("equivalence/argument-sweep/helpers.R")
source("equivalence/argument-sweep/datasets.R")
data <- sweep_punch(sweep_categorical()$esm, 0.1)
codes <- as.data.frame(lapply(data, function(v) as.integer(factor(v))))
formula <- stats::as.formula(sprintf("cbind(%s) ~ 1", paste(names(data), collapse = ", ")))
reference <- sweep_seeded(31L, function() {
  suppressWarnings(poLCA::poLCA(formula, codes, nclass = 3L, maxiter = 50000L, tol = 1e-12,
                                na.rm = FALSE, verbose = FALSE, calc.se = FALSE, nrep = 15L))
})
options(digits = 12, width = 200)
cat("poLCA best of 15:", reference$llik, "\n")
print(sort(round(reference$attempts, 4), decreasing = TRUE))
searches <- sweep_grid(n_starts = c(10L, 30L, 100L), search = c("kmeans", "emEM"))
print(do.call(rbind, lapply(searches, function(o) {
  fit <- suppressWarnings(lca(data, names(data), n_classes = 3L, missing = "fiml", seed = 21L,
                              n_starts = o$n_starts, search = o$search, tol = 1e-10,
                              max_iter = 20000L))
  starts <- sweep_table(fit, "starts")
  data.frame(n_starts = o$n_starts, search = o$search, log_likelihood = fit$log_likelihood,
             gap_to_polca = reference$llik - fit$log_likelihood,
             distinct_optima = length(unique(round(starts$log_likelihood, 3))),
             replicated = fit$n_best_replicated, boundary = fit$boundary)
})))
