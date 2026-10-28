# heavy_02 with the direct effect removed (y3 - b z at each solution's b):
# latents kmeans/emEM versus mclust's model-based hierarchical start.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(mclust))
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
vars <- case$spec$vars
fits <- readRDS("/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/heavy/fits.rds")
b_values <- vapply(fits, \(fit) unname(fit$direct_effects), numeric(1))
print(b_values)
rows <- lapply(names(b_values), function(name) {
  adjusted <- transform(case$data, y3 = y3 - b_values[[name]] * z)
  fit_latents <- \(search) suppressWarnings(lpa(adjusted, vars, 2, model = "VVI", n_starts = 10,
    seed = 1, tol = 1e-10, max_iter = 5000, search = search))$log_likelihood
  mclust_ll <- \(use) Mclust(subset(adjusted, select = vars), G = 2, modelNames = "VVI",
                             initialization = list(hcPairs = hc(subset(adjusted, select = vars),
                               modelName = "VVV", use = use)), verbose = FALSE)$loglik
  data.frame(b_from = name, latents_kmeans = fit_latents("kmeans"), latents_emEM = fit_latents("emEM"),
             mclust_hc_svd = mclust_ll("SVD"), mclust_hc_vars = mclust_ll("VARS"),
             mclust_hc_std = mclust_ll("STD"))
})
print(do.call(rbind, rows), digits = 10)
