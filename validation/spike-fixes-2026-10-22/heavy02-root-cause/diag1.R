# heavy_02: what the two maxima look like, and how varied the k-means starts are.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
sys.source(file.path(root, "reference.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
vars <- case$spec$vars
str(case$data)
str(case$spec)
str(case$truth)
fit_with <- \(search) suppressWarnings(lpa(case$data, vars, 2, model = "VVI", n_starts = 10,
  seed = 1, tol = 1e-10, max_iter = 5000, direct_effects = "y3 ~ z", search = search))
low <- fit_with("kmeans"); high <- fit_with("emEM")
cat("\nlow ", low$log_likelihood, " high ", high$log_likelihood, "\n")
show <- \(fit) { print(get_results(fit, "estimates")) }
cat("\n--- low (kmeans) ---\n"); print(names(low)); print(low$means); print(low$variances); print(low$profile_probabilities)
cat("\n--- high (emEM) ---\n"); print(high$means); print(high$variances); print(high$profile_probabilities)
saveRDS(list(low = low, high = high), "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/heavy/fits.rds")
