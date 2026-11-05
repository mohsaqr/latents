# Accuracy and rigor review, 2026-11-04: summary

Eight read-only reviewers, one report each in this folder. 0 CRITICAL, 19 MAJOR, 50 MINOR, plus EDGE items.
All 19 MAJOR findings were re-run by the parent on new seeds (probes/verify-1.R to verify-3.R; own hand or finite-difference checks where marked). MINOR and EDGE findings were not re-run.

| # | Finding | Where | Parent's reproduction | Fix changes numbers? |
|---|---|---|---|---|
| 1 | Two-step LCA SE leaves held response logits free | R/covariate-inference.R:790 | SE 0.185 vs 0.118 by own finite differences | yes |
| 2 | compare_models() penalises a staged fit for second-stage parameters only | R/compare-models.R:418 | k 3 vs 23, weight 1.000; with 23, BIC 16987.42 vs joint 16986.89 | yes, that comparison |
| 3 | Cross-level BIC selects group classes and the full form on data with no group structure | R/cross-level.R:437; vignettes/articles/additive.Rmd:280 | 6 of 6: full -84 to -132, restricted -35 to -67; ordinary family +8 to +9 | no if documented |
| 4 | Unused factor levels in membership, group_membership, random count as parameters | R/mixture-regression.R:891, :690 | k 9 vs 8, BIC +4.09, all SEs NA, equal LL | yes, affected fits |
| 5 | Near-empty origin counted and reported as estimated | R/lta-engine.R:1370; R/transitions.R:972 | row 1e-10/1/1e-10, expected moves 3e-6, estimated TRUE, k 38, no warning | yes, affected fits |
| 6 | Weights with grand centring or scaling use unweighted offsets | R/centering.R:192, :116 | means differ from duplicated rows by 0.68; 0 without the transform | yes, or document |
| 7 | Person centring turns a singly observed value into an exact zero | R/centering.R:315 | spike on 1 of 12 seeds (BIC off by 376); k3 preferred 10 of 12 vs 6 of 12 | warning no; repair yes |
| 8 | Column named `class` as id or cluster corrupts assignment tables | R/mixture-regression-methods.R:283; R/growth-mixture-methods.R:586 | ids 501... show as 2, 1, 2 | no estimates |
| 9 | Boundary flag and "estimates at a bound" miss response probabilities | R/diagnostics-checks.R:68; R/compare-lrt.R:48 | prob 1e-10, flag FALSE, check "ok", Wald refused | status and warnings |
| 10 | Unidentified transition specifications fit silently | R/transitions.R:872 | T=2 stayer: LL equal to 6 decimals, k 17 vs 15, stayer share 0.45, no warning | warning no |
| 11 | Unidentified growth fits report success | R/mixture-regression.R:318 | exact ridge confirmed by own likelihood, converged TRUE, no warning. Part (b), binary without predictors, not re-run | warning no |
| 12 | step_one = "estimated" on a staged fit omits measurement uncertainty | R/three-step.R:593 | ratio 1.003, 1.004 vs joint 1.058, 1.076 | yes, or refusal |
| 13 | Staged bootstrap with a foreign measurement | R/bootstrap-inference.R:137 | estimate outside own interval 8 of 14, 3 of 14, control 0 of 14 | refusal |
| 14 | No boundary guard on r3step chain or pooled fits | R/lta-three-step.R:264; R/three-step.R:1125 | missing guard confirmed; robust p .0005-.0007 beside observed p .12-.13 on 2 of 11 fits. Reviewer's extreme slopes did not recur; chain not re-run | NA on such fits |
| 15 | Cross-level `prevalence` means different things in the two forms | R/cross-level.R:452 | 0.515 vs 0.359; truth 0.478 per average group, 0.64 per individual | owner's choice |
| 16 | p_adjusted is raw in get_results(), BH in get_residuals() | R/get-data.R:859 | equal to p_value in one, not in the other. workflow-lca.Rmd filter count not re-run | yes, or NA |
| 17 | OPG SEs under grand transforms break invariance | R/inference.R:1165; R/covariate-inference.R:595 | 0.081 to 0.113; observed and robust unchanged | yes, narrow |
| 18 | bootstrap_lrt() fails every replicate for categorical joint pairs | R/enumeration.R:1102 | 0 of 3 valid, "undefined columns selected" | NA becomes a number |
| 19 | Transition bootstrap crashes on one longer person | R/lta-bootstrap.R:184, :197 | unclassed error, homogeneous engine completes | no |

Reviewer's suggested fixes are in each report and were not applied. One was contradicted: a 1e-6 tolerance for #5 would miss the row at 2.9e-6.

Not reviewed: regression BLRT, growth bootstrap, multilevel robust covariance, parts of R/kernel-structures.R and R/multilpa-pool.R, plot rendering; articles were read and partly run, not rendered.

Full suite: 14847 expectations, 1 failed (stale test-plot-polish.R:146, wording removed in 1.1.14), 10 errors (socket workers loading the installed 0.95.1; those files pass against a temporary install of the source), 16 skipped.
