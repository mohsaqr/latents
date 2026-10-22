# FINDING 1: EM declares convergence short of the maximum (twopair_mod_05).
# Run from the project root:
#   Rscript equivalence/features/joint-items/findings/squarem-early-stop.R
#
# twopair_mod_05: 3 classes, two binary pairs, 800 rows. Latent GOLD
# (Newton-Raphson finish, tolerance 1e-12) prints -2916.7570, Mplus -2916.757,
# and the reference EM polished from latents' own estimates reaches
# -2916.757011. latents, asked for tol = 1e-10, reports converged = TRUE at
# -2916.757139 (kit fit: 15 starts x 2 seeds, max_iter 5000), 1.28e-4 below
# the maximum (4.4e-8 relative, 440 times the requested tolerance). Its
# posteriors are up to 7.8e-3 from Latent GOLD's and Mplus's (which agree with
# each other to 5e-5), its item probabilities 3.3e-3, class sizes 9.9e-4.
# Output of this script (2026-10-20):
#   squarem, tol 1e-10: -2916.757202, 1431 iterations, converged
#   squarem, tol 1e-13: -2916.757012, 1650 iterations, converged
#   plain EM, tol 1e-10: -2916.757907, 8398 iterations, converged
#   plain EM, tol 1e-13: -2916.757014, 44438 iterations, converged
#   the pairs joined by hand, squarem, tol 1e-10: -2916.757202 (same as above)
# Diagnosis: the stopping rule is the gain per step or SQUAREM cycle,
# |gain| <= tol (1 + |LL|) (R/kernel-em.R). Along a flat direction EM's gain
# is tiny long before the maximum (remaining gap ~ gain / (1 - rate)), so the
# rule is met early. Not specific to local_dependence (the hand-joined model is
# identical); poLCA, also gain-based, stops 3.7e-5 short on these data. A
# score check or a Newton finish (as Latent GOLD does; LEARNINGS.md
# 2026-10-19 "EM near the optimum") would close it. Over the other 53 kit
# cases latents' gap was at most 5.3e-10 relative, probabilities 8.7e-5,
# posteriors 4.1e-4. In the reference study (140 datasets, 8 starts) 7 fits
# end 1.8e-6 to 1.04e-5 below the reference maximum (checked for the three
# largest: the same mode, reached by polishing latents' estimates), e.g. twopair_mod
# seed 50086: latents -1833.0852513 (converged, 189 iterations), the same
# estimates polished -1833.0852409, the reference EM's best -1833.0852409.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
case <- readRDS(file.path("equivalence", "features", "joint-items", "kit", "cases.rds"))$twopair_mod_05
fit <- function(...) {
  withCallingHandlers(lca(case$data, case$spec$vars, 3, local_dependence = case$spec$pairs,
                          n_starts = 15, seed = 1, max_iter = 100000, ...),
                      message = \(m) invokeRestart("muffleMessage"))
}
runs <- list(squarem_1e10 = fit(tol = 1e-10), squarem_1e13 = fit(tol = 1e-13),
             none_1e10 = fit(tol = 1e-10, acceleration = "none"),
             none_1e13 = fit(tol = 1e-13, acceleration = "none"))
print(data.frame(run = names(runs),
                 log_likelihood = sprintf("%.6f", vapply(runs, `[[`, 0, "log_likelihood")),
                 iterations = vapply(runs, `[[`, 0, "iterations"),
                 converged = vapply(runs, \(f) isTRUE(f$converged), TRUE)))
# The same model fitted by hand (the two pairs joined beforehand): not specific
# to local_dependence.
by_hand <- transform(case$data, pq1 = paste(p1, p2), pq2 = paste(q1, q2))
hand <- withCallingHandlers(lca(by_hand, c("pq1", "pq2", "u1", "u2"), 3, n_starts = 15, seed = 1,
                                tol = 1e-10, max_iter = 100000),
                            message = \(m) invokeRestart("muffleMessage"))
cat(sprintf("joined by hand, squarem, tol 1e-10: %.6f (%d iterations)\n",
            hand$log_likelihood, hand$iterations))
out <- readLines(file.path("equivalence", "features", "joint-items", "returned", "twopair_mod_05.out"))
cat("Latent GOLD prints -2916.7570 (4 decimals); Mplus:", trimws(grep("H0 Value", out, value = TRUE)), "\n")
