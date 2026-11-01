# Compare latents' three_step(method = "ml", step_one = ) with Latent GOLD's
# step-three listings and its own step-one classification-error files. Run
# from the project root after run-kit.sh; writes comparison.csv.
#
# Rows:
# - step-one error covariance: the covariance of the ~err logits e from the
#   independent likelihood here (J V J') against Latent GOLD's own step one
#   (its Step1Jacobian, Step1Hessian), matched by label. Both are invariant to
#   how the step-one model is parameterized. Latent GOLD reads the logits
#   with 16 digits; the tolerance is relative, 1e-3 of the largest entry, for
#   two optimizers' maxima and numerical Jacobians.
# - step-three class means / probabilities and robust standard errors,
#   latents (vcov_type = "independent") against Latent GOLD's robust ones,
#   with the fixed parameters read as plain values (step_one = "fixed") and
#   with their covariance (step_one = "estimated"), at half a unit in the
#   printed place plus 1e-4 of the value for the step-one optimizer.
# - the overall Wald test of the class effect, the same two ways.
pkgload::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-step-one")
source(file.path("equivalence", "latentgold", "compare-functions.R"))
source(file.path(base, "common.R"))
comparison <- step_comparison(file.path(base, "returned"),
                              readRDS(file.path(base, "kit", "independent.rds")))
utils::write.csv(comparison, file.path(base, "comparison.csv"), row.names = FALSE)
print(table(comparison$status, comparison$step_one))
print(subset(comparison, status != "agree"), row.names = FALSE)
