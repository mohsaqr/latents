source("/Users/mohammedsaqr/Documents/Github/latents/validation/review-2026-11-04/probes/three-step-mutants-common.R")
setwd("/Users/mohammedsaqr/Documents/Github/latents")
assign("root", "/Users/mohammedsaqr/Documents/Github/latents/tests/equivalence", envir=globalenv())
f <- "test-lta-three-step.R"
mutate("lta: robust meat doubled", ".lta_three_step_fit",
  ".inference_sandwich(inverse, crossprod(scores), nrow(scores))", ".inference_sandwich(inverse, 2 * crossprod(scores), nrow(scores))", f)
mutate("lta: adjust includes intercepts", ".lta_three_step",
  "tested <- fitted$labels$term != \"(Intercept)\"", "tested <- rep(TRUE, length(p_value))", f)
mutate("ML sizes uniform instead of avg posterior", ".multilpa_step_ml_fit",
  "prior <- colMeans(pieces$posteriors)", "prior <- rep(1/pieces$n_classes, pieces$n_classes)", "test-three-step-ml.R")
mutate("r3step meat not clustered", "r3step",
  "rowsum(score_rows(fitted$par), pieces$group_index, reorder = FALSE)", "score_rows(fitted$par)", "test-tidy-three-step.R")
assign("root", "/Users/mohammedsaqr/Documents/Github/latents/equivalence", envir=globalenv())
mutate("r3step meat not clustered", "r3step",
  "rowsum(score_rows(fitted$par), pieces$group_index, reorder = FALSE)", "score_rows(fitted$par)", "test-jstats-mplus.R")
