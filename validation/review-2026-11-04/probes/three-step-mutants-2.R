source("/Users/mohammedsaqr/Documents/Github/latents/validation/review-2026-11-04/probes/three-step-mutants-common.R")
bg <- c("test-latents-controls.R","test-review-1015-covariates.R","test-audit-1101-transitions-three-step.R")
vapply(bg, function(f){ t0<-proc.time()[[3]]; run_files(f); proc.time()[[3]]-t0 },0) |> print()
mutate("by_group intercept not back-shifted", ".multilpa_r3step_group_model",
  "(shift %*% slope_inverse)", "(0 * shift %*% slope_inverse)", bg)
mutate("by_group robust meat doubled", ".multilpa_r3step_by_group_class",
  ".multilpa_cross_product(scores)", "2 * .multilpa_cross_product(scores)", bg)
mutate("by_group family includes group logit", ".multilpa_r3step_group_frame",
  "tested <- labels$level == \"individuals\" & labels$term %in% terms", "tested <- !is.na(raw)", bg)
