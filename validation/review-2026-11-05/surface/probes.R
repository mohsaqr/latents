Sys.setenv(NOT_CRAN = "true", VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")
pkgload::load_all(".", quiet = TRUE)
warnings_seen <- character()
capture_warnings <- function(expr) withCallingHandlers(expr, warning = function(w) {
  warnings_seen <<- c(warnings_seen, paste(class(w)[1L], conditionMessage(w)))
  invokeRestart("muffleWarning")
})
cat("Working version:", as.character(utils::packageVersion("latents")), "\n")
set.seed(1106)
n <- 400L
z <- sample(1:2, n, replace = TRUE)
d <- data.frame(a = ifelse(z == 1L, "x", "y"),
                b = ifelse(z == 1L, "u", "v"),
                c = ifelse(runif(n) < c(.85, .2)[z], "r", "s"),
                d = ifelse(runif(n) < c(.9, .3)[z], "p", "q"),
                e = ifelse(runif(n) < c(.75, .2)[z], "m", "n"))
str(d)
fit <- capture_warnings(lca(d, names(d), n_classes = 2, n_starts = 6, seed = 1105,
                            tol = 1e-10, max_iter = 3000))
joined <- capture_warnings(lca(d, names(d), n_classes = 2, n_starts = 6, seed = 1105,
                               local_dependence = "c:d", tol = 1e-10, max_iter = 3000))
cat("\nBOUNDARY PROBE\n")
print(fit$response_probabilities$a)
cat("min probability", min(unlist(fit$response_probabilities)),
    "floor", .multilpa_probability_floor(fit), "boundary", fit$boundary, "\n")
checks <- get_checks(fit)
print(checks[checks$check == "estimates at a bound", ])
warnings_seen <- character()
lrt <- capture_warnings(compare_models(independent = fit, joined = joined, test = "lrt"))
print(as.data.frame(lrt)[, c("model", "lr_statistic", "lr_df", "lr_p_value")])
cat("LRT warnings:", paste(warnings_seen, collapse = " | "), "\n")
infer <- tryCatch(parameter_inference(fit), error = identity)
cat("Inference result:", class(infer)[1L], if (inherits(infer, "error")) conditionMessage(infer), "\n")
stopifnot(!fit$boundary, any(unlist(fit$response_probabilities) <= .multilpa_probability_floor(fit)),
          all(checks$status[checks$check == "estimates at a bound"] == "ok"),
          inherits(infer, "latents_boundary_fit"), !length(warnings_seen))

cat("\nCATEGORICAL BLRT NAMES PROBE\n")
set.seed(1116)
id <- rep(seq_len(45L), each = 5L)
g <- sample(1:2, 45L, replace = TRUE)
k <- 1L + rbinom(length(id), 1L, c(.15, .85)[g[id]])
pair_data <- data.frame(id = id, a = ifelse(runif(length(id)) < c(.2,.8)[k], "x", "y"),
                        b = ifelse(runif(length(id)) < c(.7,.3)[k], "u", "v"),
                        y = rnorm(length(id), c(-1.5,1.5)[k]))
str(pair_data)
make_pair <- function(groups) capture_warnings(multilpa(pair_data, c("a","b","y"), "id",
  n_profiles = 2, n_group_classes = groups, categorical = c("a","b"),
  local_dependence = "a:b", n_starts = 3, seed = 1105, max_iter = 3000, tol = 1e-8))
pair_null <- make_pair(1)
pair_alt <- make_pair(2)
cat("Converged", pair_null$converged, pair_alt$converged, "vars", pair_null$vars, "\n")
simulation <- tryCatch(.multilpa_simulate(pair_null), error = identity)
cat("Direct simulation:", class(simulation)[1L], if(inherits(simulation,"error")) conditionMessage(simulation), "\n")
blrt <- capture_warnings(bootstrap_lrt(pair_null, pair_alt, iter = 3, seed = 1105))
print(as.data.frame(blrt))
print(get_results(blrt, "replicates"))
stopifnot(inherits(simulation, "error"), as.data.frame(blrt)$n_valid == 0L,
          all(grepl("undefined columns selected", get_results(blrt, "replicates")$error)))

cat("\nCROSS LEVEL PREVALENCE AND ZERO COUNTS\n")
set.seed(1127)
groups <- seq_len(60L)
type <- rep(1:2, length.out = length(groups))
cross_data <- do.call(rbind, lapply(groups, function(j) {
  n_j <- c(4L, 35L)[type[j]]
  profile <- sample(1:2, n_j, replace = TRUE, prob = if(type[j] == 1) c(.85,.15) else c(.25,.75))
  data.frame(team = j, y1 = rnorm(n_j, c(-1.5,1.5)[profile], .6),
             y2 = rnorm(n_j, c(-1,1)[profile], .6))
}))
str(cross_data)
make_cross <- function(family, data = cross_data, weights = NULL) capture_warnings(
  multilpa(data, c("y1","y2"), "team", n_profiles = 2, n_group_classes = 2,
            family = family, n_starts = 3, seed = 1105, weights = weights, tol = 1e-9))
warnings_seen <- character()
full <- make_cross("full_cross_level")
cat("Full warning classes:", paste(warnings_seen, collapse = " | "), "\n")
restricted <- make_cross("restricted_cross_level")
full_prevalence <- unique(subset(get_results(full,"profiles"), level == "individual")[c("profile","prevalence")])
restricted_prevalence <- unique(subset(get_results(restricted,"profiles"), level == "individual")[c("profile","prevalence")])
cat("Full table:\n"); print(full_prevalence)
cat("Full posterior means:", colMeans(full$subject_posteriors), "\n")
cat("Restricted table:\n"); print(restricted_prevalence)
cat("Restricted posterior means:", colMeans(restricted$subject_posteriors), "\n")
stopifnot(max(abs(full_prevalence$prevalence - colMeans(full$subject_posteriors))) > .1,
          max(abs(restricted_prevalence$prevalence - colMeans(restricted$subject_posteriors))) < 1e-4)
cross_data$w <- ifelse(cross_data$team <= 6L, 0, 1)
weighted <- make_cross("full_cross_level", weights = "w")
subset_fit <- make_cross("full_cross_level", data = subset(cross_data, w > 0))
cat("Zero-weight fit table:\n"); print(get_results(weighted, "fit"))
cat("Subset fit table:\n"); print(get_results(subset_fit, "fit"))
cat("nobs", nobs(weighted), "LL difference", weighted$log_likelihood - subset_fit$log_likelihood,
    "BIC difference", weighted$bic - subset_fit$bic, "\n")
stopifnot(get_results(weighted,"fit")$n_groups == 60L, nobs(weighted) == 54L,
          abs(weighted$log_likelihood - subset_fit$log_likelihood) < .01)

cat("\nLATEST FIX AND WILDCARD TESTS\n")
results <- lapply(c("test-review-1105-fixes.R", "test-wildcards.R", "test-plot-docs.R"), function(file) {
  testthat::test_file(file.path("tests/testthat", file), reporter = "summary", stop_on_failure = FALSE)
})
totals <- do.call(rbind, lapply(results, as.data.frame))
print(aggregate(cbind(failed,error,warning,skipped,passed) ~ 1, totals, sum))
saveRDS(results, "validation/review-2026-11-05/surface/test-results.rds")
cat("\nSession:\n"); print(sessionInfo())
