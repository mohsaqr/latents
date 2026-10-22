# Finding 1: three_step(method = "ml") returns a local maximum of the
# step-three likelihood, without a warning, when that likelihood has several.
#
# From the project root:
#   Rscript equivalence/features/ml-three-step/repro-local-maximum.R
#
# sweep_105 of the kit: three profiles (simulab LPA, n = 1000, entropy 0.64)
# and a heavy-tailed outcome (t with 3 df). latents runs EM from the
# proportional moments (Latent GOLD's non-random start) and stops at the
# first maximum; the independent reference (25 starts) and Latent GOLD with
# 100 random start sets reach a higher one.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
sys.source(file.path("equivalence", "features", "ml-three-step", "helpers.R"),
           envir = environment())
sys.source(file.path("equivalence", "features", "ml-three-step", "reference.R"),
           envir = environment())
case <- readRDS(file.path(ml3_kit(), "cases.rds"))$sweep_105
fit <- suppressWarnings(ml3_fit(case))
ours <- three_step(fit, case$data, "y", method = "ml", vcov_type = "independent")
print(ours)
pieces <- .multilpa_level_assignments(fit, "individuals")
inner <- .multilpa_step_ml_fit(pieces, case$data$y, "varying")
reference <- ref_step3(unname(fit$subject_posteriors), case$data$y, "varying")
cat(sprintf("latents log likelihood   %.4f (converged: %s)\n", inner$log_likelihood,
            inner$converged))
cat(sprintf("reference log likelihood %.4f (its starts: %s)\n", reference$log_likelihood,
            paste(sprintf("%.3f", reference$all_log_likelihoods), collapse = ", ")))
random <- ml3_lg_result("x_sweep_105", "vr", 3L, "varying")
default <- ml3_lg_result("sweep_105", "vr", 3L, "varying")
if (!is.null(random)) {
  cat(sprintf("Latent GOLD default start %.4f, 100 random start sets %.4f\n",
              default$log_likelihood, random$log_likelihood))
}
print(data.frame(class = 1:3, latents = ours$estimate, reference = reference$classes$estimate,
                 reference_se = reference$classes$standard_error,
                 shift_in_se = (ours$estimate - reference$classes$estimate) /
                   reference$classes$standard_error))
