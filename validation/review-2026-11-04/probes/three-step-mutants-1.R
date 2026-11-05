source("/Users/mohammedsaqr/Documents/Github/latents/validation/review-2026-11-04/probes/three-step-mutants-common.R")
fast <- c("test-three-step.R","test-three-step-clusters.R","test-three-step-ml.R","test-three-step-categorical.R","test-review-1101-three-step-categorical.R","test-audit-1101-transitions-three-step.R")
vapply(fast, function(f){ t0<-proc.time()[[3]]; run_files(f); proc.time()[[3]]-t0 },0) |> print()
mutate("pair SE assumes independence", ".multilpa_step_pairs",
  "(clustered[, pairs[2L, ], drop = FALSE] - clustered[, pairs[1L, ], drop = FALSE])^2",
  "clustered[, pairs[2L, ], drop = FALSE]^2 + clustered[, pairs[1L, ], drop = FALSE]^2", fast)
mutate("r3step meat not clustered", "r3step",
  "rowsum(score_rows(fitted$par), pieces$group_index, reorder = FALSE)", "score_rows(fitted$par)",
  c("test-three-step.R","test-three-step-clusters.R","test-audit-1101-transitions-three-step.R"))
mutate("categorical overall: no 1/p delta", ".multilpa_three_step_categorical",
  "sweep(clustered[, columns, drop = FALSE], 2L, p, \"/\")", "clustered[, columns, drop = FALSE]",
  c("test-three-step-categorical.R","test-review-1101-three-step-categorical.R","test-audit-1101-transitions-three-step.R"))
mutate("categorical CI half without 1/p(1-p)", ".multilpa_three_step_categorical",
  "half <- quantile * errors/(estimates * (1 - estimates))", "half <- quantile * errors",
  c("test-three-step-categorical.R","test-review-1101-three-step-categorical.R"))
mutate("ML sizes uniform instead of avg posterior", ".multilpa_step_ml_fit",
  "prior <- colMeans(pieces$posteriors)", "prior <- rep(1/pieces$n_classes, pieces$n_classes)",
  c("test-three-step-ml.R","test-audit-1101-transitions-three-step.R"))
