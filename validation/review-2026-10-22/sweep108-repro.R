# Run from the repository root with single-threaded BLAS.
pkgload::load_all(quiet = TRUE)
source("equivalence/features/ml-three-step/reference.R")
case <- readRDS("equivalence/features/ml-three-step/kit/cases.rds")$sweep_108
pieces <- list(posteriors = case$posteriors, modal = case$modal,
               n_classes = case$spec$k)
current <- latents:::.multilpa_step_ml_fit(pieces, case$data$y, "equal")
reference <- ref_step3(case$posteriors, case$data$y, "equal")
se <- sqrt(colSums(latents:::.multilpa_step_ml_influence(
  current, case$data$y, "equal")^2) * nrow(case$data) / (nrow(case$data) - 1))
options(digits = 14)
cat("Known documented residual: HEAD 12e0145; not a novel regression.\n")
print(c(current_ll = current$log_likelihood,
        reference_ll = reference$log_likelihood,
        ll_gap = reference$log_likelihood - current$log_likelihood,
        current_converged = current$converged,
        current_at_floor = current$at_floor,
        reference_score = reference$max_score))
print(data.frame(class = seq_along(current$means), current = current$means,
                 reference = reference$classes$estimate,
                 current_robust_se = se,
                 reference_robust_se = reference$classes$standard_error,
                 reference_model_se = reference$classes$standard_error_model))
saveRDS(list(commit = "12e0145", current = current, reference = reference,
             current_robust_se = se),
        "/private/tmp/latents-review-20261022/sweep108-results.rds")
