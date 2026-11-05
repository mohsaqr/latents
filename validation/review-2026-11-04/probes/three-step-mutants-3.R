source("/Users/mohammedsaqr/Documents/Github/latents/validation/review-2026-11-04/probes/three-step-mutants-common.R")
src <- deparse(get(".lta_three_step_inputs",ns),width.cutoff=500L); print(grep("errors",src,value=TRUE))
lt <- c("test-lta-three-step.R","test-audit-1101-transitions-three-step.R")
mutate("lta: D transposed in emissions", ".lta_three_step_inputs",
  "t(errors)[modal, , drop = FALSE]", "errors[modal, , drop = FALSE]", lt)
mutate("lta: chain posteriors, not occasion", ".lta_occasion_posteriors",
  "joint/rowSums(joint)", "posteriors", lt)
mutate("lta: robust meat doubled", ".lta_three_step_fit",
  ".inference_sandwich(inverse, crossprod(scores), nrow(scores))", ".inference_sandwich(inverse, 2 * crossprod(scores), nrow(scores))", lt)
mutate("lta: adjust includes intercepts", ".lta_three_step",
  "tested <- fitted$labels$term != \"(Intercept)\"", "tested <- rep(TRUE, length(p_value))", lt)
so <- "test-three-step-step-one.R"
mutate("step-one ML: sizes not updated", ".multilpa_three_step_ml",
  "prior <- colMeans(updated$posteriors)", "prior <- colMeans(pieces$posteriors)", so)
mutate("step-one: V1 scaled by 0 (no correction)", ".multilpa_covariance_root",
  "decomposition$vectors %*% diag(sqrt(pmax(values, 0)), length(values))", "0 * decomposition$vectors", so)
