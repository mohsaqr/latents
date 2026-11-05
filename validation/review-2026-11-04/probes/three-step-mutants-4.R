source("/Users/mohammedsaqr/Documents/Github/latents/validation/review-2026-11-04/probes/three-step-mutants-common.R")
root <- "/Users/mohammedsaqr/Documents/Github/latents/equivalence"; assign("root", root, envir=globalenv())
setwd("/Users/mohammedsaqr/Documents/Github/latents")
t0 <- proc.time()[[3]]
mutate("categorical overall: no 1/p delta", ".multilpa_three_step_categorical",
  "sweep(clustered[, columns, drop = FALSE], 2L, p, \"/\")", "clustered[, columns, drop = FALSE]",
  "test-latentgold-categorical-distal.R")
mutate("categorical CI half without 1/p(1-p)", ".multilpa_three_step_categorical",
  "half <- quantile * errors/(estimates * (1 - estimates))", "half <- quantile * errors",
  "test-latentgold-categorical-distal.R")
mutate("ML sizes uniform instead of avg posterior", ".multilpa_step_ml_fit",
  "prior <- colMeans(pieces$posteriors)", "prior <- rep(1/pieces$n_classes, pieces$n_classes)",
  "test-features-ml-three-step.R")
mutate("r3step meat not clustered", "r3step",
  "rowsum(score_rows(fitted$par), pieces$group_index, reorder = FALSE)", "score_rows(fitted$par)",
  "test-jstats-mplus.R")
cat("elapsed",proc.time()[[3]]-t0,"\n")
