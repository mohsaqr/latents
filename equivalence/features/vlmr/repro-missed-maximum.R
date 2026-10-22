# Reproduction: random starts miss a regular maximum that Latent GOLD finds.
# Run from the project root: Rscript equivalence/features/vlmr/repro-missed-maximum.R
#
# eee_k2_04 (simulab, n = 500, 3 indicators, EEE), 3 classes. Latent GOLD
# (200 start sets) reports LL -2220.9820 with class shares 0.51 / 0.13 / 0.35
# and no variance or share near a bound. The independent likelihood of
# reference.R evaluated at Latent GOLD's printed estimates gives the same
# value (-2220.9822 at the rounded estimates). latents reaches only -2222.2043
# from 60, 150 or 600 random starts (seed 11),
# but converges to -2220.9822 when started from Latent GOLD's estimates, so
# the maximum is regular and reachable: the random starts never land in its
# basin. The 2-vs-3 LR is then 2.44 too small and the 3-vs-4 LR wrong too.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
here <- file.path("equivalence", "features", "vlmr")
sys.source(file.path(here, "reference.R"), envir = environment())
case <- readRDS(file.path(here, "kit", "cases.rds"))[["eee_k2_04"]]
vars <- case$spec$vars
x <- as.matrix(case$data[vars])

shares <- exp(c(0, -1.3505, -0.3808))
shares <- shares / sum(shares)
means <- cbind(0.2750 + c(0, -1.1211, 1.9046), -0.0186 + c(0, 0.1954, 1.7324),
               0.1819 + c(0, -0.6036, 1.3416))
covariance <- matrix(c(0.7891, 0.4521, 0.1374, 0.4521, 1.3563, 0.3755,
                       0.1374, 0.3755, 0.6511), 3)
latent_gold <- list(model = "EEE", k = 3L, proportions = shares, means = means,
                    covariance = covariance)
cat(sprintf("reference LL at Latent GOLD's printed estimates: %.4f\n",
            sum(ref_unit_loglik(ref_pack(latent_gold), x, "EEE", 3L))))

random <- vapply(c(60L, 150L, 600L), function(n_starts) {
  suppressWarnings(lpa(case$data, vars, n_profiles = 3, model = "EEE", n_starts = n_starts,
                       seed = 11, tol = 1e-12, max_iter = 5000))$log_likelihood
}, numeric(1))
cat(sprintf("latents from %d random starts: %.4f\n", c(60L, 150L, 600L), random), sep = "")

start <- list(means = means, variances = matrix(diag(covariance), 3, 3, byrow = TRUE),
              covariances = array(covariance, c(3, 3, 3)),
              profile_probabilities = matrix(shares, 1), group_probabilities = 1)
started <- suppressWarnings(lpa(case$data, vars, n_profiles = 3, model = "EEE", start = start,
                                n_starts = 1, seed = 1, tol = 1e-12, max_iter = 5000))
cat(sprintf("latents started at Latent GOLD's estimates: %.4f (converged %s)\n",
            started$log_likelihood, started$converged))
