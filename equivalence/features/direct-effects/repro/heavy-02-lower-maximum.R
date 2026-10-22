# heavy_02: Latent GOLD and Mplus both reach log likelihood -3098.4868; latents
# stops at -3106.9860 from every start, with or without the direct effect.
# From the project root (after make-kit.R; reads the retained LG output):
#   Rscript equivalence/features/direct-effects/repro/heavy-02-lower-maximum.R
pkgload::load_all(".", quiet = TRUE)
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
sys.source(file.path(root, "reference.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
vars <- case$spec$vars
fits <- lapply(1:3, \(seed) suppressWarnings(lpa(case$data, vars, 2, model = "VVI",
                                                 n_starts = 200, seed = seed, tol = 1e-10,
                                                 max_iter = 5000, direct_effects = "y3 ~ z")))
latents_values <- vapply(fits, `[[`, numeric(1), "log_likelihood")
reached <- vapply(fits, \(fit) sum(abs(fit$starts$log_likelihood - max(latents_values)) < 1e-4),
                  numeric(1))
# Latent GOLD's estimates, evaluated by the independent likelihood (reference.R).
estimates <- de_latentgold("heavy_02", case, root)$estimates
by_class <- \(kind) rbind(estimates[sprintf("%s:%s:1", kind, vars)],
                          estimates[sprintf("%s:%s:2", kind, vars)])
point <- list(means = by_class("mean"), variances = by_class("variance"),
              effects = case$truth$effects,
              proportions = unname(estimates[c("proportion:1", "proportion:2")]))
point$effects$b <- unname(estimates[["effect:y3:z"]])
colnames(point$means) <- colnames(point$variances) <- vars
stationary <- Reduce(\(p, i) ref_m_step(ref_posteriors(p, case$data, vars), p, case$data, vars,
                                        "VVI"), seq_len(500), point)
print(data.frame(
  quantity = c(sprintf("latents, 200 starts, seed %d", 1:3),
               "independent likelihood at Latent GOLD's estimates",
               "after 500 further EM steps from there (a fixed point)",
               "Mplus (retained output)"),
  log_likelihood = c(latents_values, ref_log_likelihood(point, case$data, vars),
                     ref_log_likelihood(stationary, case$data, vars),
                     de_mplus("heavy_02", case, root)$log_likelihood),
  starts_reaching_it = c(reached, NA, NA, NA)))
print(point$variances)
