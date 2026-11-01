# Step-one uncertainty in three_step() against Latent GOLD 6.1. Run from the
# project root:
#   Rscript equivalence/latentgold-step-one/make-kit.R
#   sh equivalence/latentgold-step-one/run-kit.sh
#   Rscript equivalence/latentgold-step-one/compare.R
#
# For each case (common.R): latents fits step one; an independent likelihood
# gives the step-one covariance and the Jacobian of the fixed step-three
# parameters (error and size logits); Latent GOLD fits step three with them
# fixed, once reading their covariance (variant "estimated") and once as
# plain values ("fixed"), with model-based and robust standard errors, and
# with the case gradients too ("second", Latent GOLD's second-order correction). Latent
# GOLD also fits step one itself and writes its own classification-error file.
# kit/independent.rds keeps the independent quantities for compare.R;
# returned/INPUTS.md5 fingerprints the kit Latent GOLD was given.
pkgload::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-step-one")
source(file.path(base, "common.R"))
kit_dir <- file.path(base, "kit")
dir.create(kit_dir, recursive = TRUE, showWarnings = FALSE)
independent <- lapply(stats::setNames(step_cases, step_cases), function(case) {
  made <- step_case(case)
  model <- step_one_model(made)
  estimate <- step_one_estimate(made, model)
  modal <- made$fit$subject_profiles
  # The values latents itself holds fixed; the Jacobian and covariance are the
  # independent likelihood's.
  values <- step_values(made$fit$subject_posteriors, modal, made$k)
  jacobian <- numDeriv::jacobian(function(theta) {
    step_fixed_values(theta, model, modal, made$k)
  }, estimate$theta)
  utils::write.table(step_data_file(made), file.path(kit_dir, sprintf("%s.dat", case)),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(cbind(id = seq_len(nrow(made$data)),
                           subset(made$data, select = made$indicators)),
                     file.path(kit_dir, sprintf("%s_step1.dat", case)),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  gradients <- numDeriv::jacobian(model$unit_log_likelihood, estimate$theta)
  writeLines(step_fixed_file(values, jacobian, estimate$covariance),
             file.path(kit_dir, sprintf("%s_fixed.txt", case)))
  writeLines(step_fixed_file(values, jacobian, estimate$covariance, gradients),
             file.path(kit_dir, sprintf("%s_fixed_gradients.txt", case)))
  writeLines(step_one_syntax(made), file.path(kit_dir, sprintf("%s_step1.lgs", case)))
  lapply(step_outcomes, function(outcome) {
    lapply(c("estimated", "second", "fixed"), function(variant) {
      lapply(c(FALSE, TRUE), function(robust) {
        file <- sprintf("%s_%s_%s%s.lgs", case, outcome, variant, if (robust) "_robust" else "")
        writeLines(step_three_syntax(case, outcome, variant, robust, made$k, values),
                   file.path(kit_dir, file))
      })
    })
  })
  list(values = values, jacobian = jacobian, covariance = estimate$covariance,
       gradients = gradients,
       log_likelihood = estimate$log_likelihood, latents_log_likelihood = made$fit$log_likelihood)
})
saveRDS(independent, file.path(kit_dir, "independent.rds"), version = 2)
invisible(lapply(names(independent), function(case) {
  cat(sprintf("%s: independent LL %.8f, latents LL %.8f, %d fixed values, %d step-one parameters\n",
              case, independent[[case]]$log_likelihood, independent[[case]]$latents_log_likelihood,
              length(independent[[case]]$values), ncol(independent[[case]]$jacobian)))
}))
files <- sort(list.files(kit_dir, pattern = "[.](dat|lgs|txt)$", full.names = TRUE),
              method = "radix")
writeLines(sprintf("%s  %s", unname(tools::md5sum(files)), basename(files)),
           file.path(kit_dir, "MANIFEST.md5"))
cat("kit written:", length(files), "files\n")
