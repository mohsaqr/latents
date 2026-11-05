# Review 2026-11-04: `three-step`

Probes: `validation/review-2026-11-04/probes/three-step-*.R`. Every number below came from a probe run in this session on R/ as it stands (1.1.15, working tree), with cores = 1 through the two-slot wrapper. Latent GOLD itself was not run; its retained outputs were used only through `equivalence/test-latentgold-*.R` in the mutation checks.

## 1. What was reviewed

| File | Depth |
|---|---|
| R/three-step.R (1522 lines) | Read in full. D, BCH weights, weighted means, sandwich, thin-class rule, pairs, step-one closure, r3step pooled, by_group_class, refusals. Probed. |
| R/three-step-ml.R (553) | Read in full. Likelihood, scores, Hessian, polish, influence, correction transform. Probed. Start-selection logic read, and checked only on nine bounded-outcome data sets. |
| R/three-step-categorical.R (422) | Read in full. BCH/modal shares, EM, information, Jacobian, overall Wald on logits, pairs. Probed. |
| R/lta-three-step.R (332) | Read in full. Step-two posteriors, basis, likelihood, scores, fit, table. Brute-force and independent ML. |
| Helpers used by these | `.inference_sandwich`, `.inference_score_information`, `.multilpa_observed_hessian`, `.multilpa_design_basis`, `.mixture_invert`, `.lta_designs`, `.lta_log_transition`, `.multilpa_sequence_layout`, `.multilpa_sequence_emission`, read only as far as these paths use them. |

Read first, not repeated: HANDOFF, LEARNINGS (step-one, Latent GOLD, ML starts), validation/audit-2026-11-01/transitions-three-step.md, the verification matrix (`three_step() / r3step()` is `++++`).

NOT reviewed: the trajectory inputs (`.trajectory_three_step_inputs`, `.trajectory_r3step_inputs`, growth-mixture area); the coverage studies (not re-run); `drop_empty = FALSE` paths; V1 charts for ordinal and count indicators (VVV and binary LCA checked here; the 2026-11-01 audit checked the base point only for those kinds); `.multilpa_step_ml_fit`'s behaviour above K = 3.

## 2. Findings, most severe first

### F1. MAJOR, CONFIRMED. A membership or transition logit heading to infinity is reported as a finite estimate with a small sandwich SE and a significant p, no warning

R/lta-three-step.R:264-281 and :316 (chain), R/three-step.R:1125-1136 (pooled r3step, same absence at the by_group_class path :1339-1348).

One-step covariate fits refuse this: `latents_extreme_coefficients`, "a standard error at a coefficient heading to infinity is meaningless" (R/covariates.R:792-801). No three-step path has any such check. The only guards are `.multilpa_observed_hessian` (refuses a singular information) and `.mixture_invert`; a nearly flat direction passes both. The default `vcov_type = "robust"` makes it worse: the scores of a logit at its boundary vanish, so the sandwich collapses while the observed-information SE explodes. The two SEs differ by up to two orders of magnitude and nothing tells the reader.

Evidence (same estimate, robust vs observed; no warning in any row):

| Case | Row | Estimate | Robust SE, p | Observed SE, p |
|---|---|---|---|---|
| Chain, K = 3, n = 80, T = 4, seed 6 | move 1 to 3, intercept | -14.39 | 4.83, p = 0.0029 | 663.3, p = 0.98 |
| same | move 1 to 3, slope z | 0.553 | 0.751 | 197.4 |
| r3step level = "groups", G = 25, seed 16 | slope z | -46.07 | 18.6, p = 0.013 | 76.7, p = 0.55 |
| same, seed 17 | slope z | -22.76 | 8.25, p = 0.0058 | 23.0, p = 0.32 |
| same, seed 13 | slope z | 14.17 | 3.75, p = 0.00016 | 9.29, p = 0.13 |
| Pooled r3step, n = 400, covariate built to nearly separate the assignment (noise sd 0.20) | slope | 367.8 | 14.8, p = 5e-136 | 36971, p = 0.99 |
| same, noise sd 0.10 | slope | 323.3 | 53.2, p = 1e-9 | 5.6e26, p = 1 |

In the 14 small-chain data sets: 4 refused `latents_no_converge` (cause inspected on one n = 150 case only, `three-step-chain-refusal.R`: BFGS hit maxit with |score| 0.001 and coefficients near +-200, and more iterations reach code 0 at -234), 1 flagged singular with NA SEs, 2 returned |estimate| > 10 silently. So the outcome on the same kind of data is a refusal, a warning or a silent number depending on where BFGS happens to stop. The group rows have no constructed link between covariate and assignment: a true slope of 4 on 25 groups. On the last two rows the one-step fit with the same covariate warned `latents_extreme_coefficients`; r3step gave no condition. At noise sd 0.30 and 0.45 (slopes 13.3 and 6.1) robust and observed SEs agree within 11% and neither path warns, so the problem is confined to the boundary. A data set with the covariate nearly exact on another draw was refused `latents_singular_information`, as it should be.

Who is affected: K >= 3 chains with a rare transition cell, and small numbers of groups or persons with a strong covariate.

Smallest fix: after the fit, flag any logit whose fitted probability is below, say, 1e-5 on every row of its design (a covariate-free -14.39 has p = 5.6e-7, under the existing `limit = 20`, so that limit alone would not catch this), warn with the existing class, and report NA SE and p for that block, or give the observed SE beside the sandwich. Would change reported numbers: yes, only on such fits (NA instead of a number). Owner decides the threshold.

Probes: three-step-chain-separation.R, three-step-separation-pooled.R, three-step-separation-groups.R, three-step-chain-refusal.R.

### F2. MAJOR, CONFIRMED. `step_one = "estimated"` on a staged or measurement-held fit silently leaves out the measurement uncertainty

R/three-step.R:593 and :601-608 (V1 from `vcov()` of a fit whose free coordinates exclude the held measurement); docs R/three-step.R:289-293 and :351-361.

`fit_staged()` returns class `multilpa`, so it passes every check. Its V1 has only the structural coordinates (3 by 3 here, against 15 by 15 for the same model estimated jointly). The classification-error uncertainty the option exists to add comes almost wholly from the means and variances. The attribute still says `step_one = "estimated"`. `?fit_staged` says Wald errors omit the first stage, but `?three_step` promises the option adds "the uncertainty of the classification-error probabilities" for "profile fits ... when their Wald inference is available".

Evidence, entropy about 0.6, G = 100, m = 10, ratio SE(estimated)/SE(fixed):

| Fit | BCH | ML |
|---|---|---|
| Joint `multilpa()` | 1.064, 1.084 | 1.071, 1.078 |
| `fit_staged()` | 1.003, 1.004 | 1.003, 1.004 |

Who is affected: anyone combining `fit_staged()` (or a held `measurement`) with the correction. No test covers the combination (grep of test-three-step*.R).

Smallest fix: refuse with `latents_unsupported_three_step` when the fit holds measurement, or propagate the first stage's V1, or at least warn and document. Would change reported numbers: refusal removes a table; propagation raises those SEs by the gap above.

Probe: three-step-staged-stepone.R. It reports on one data set; the ratio of the gap will vary with separation.

### F3. MINOR, CONFIRMED. Categorical pairs: non-BH corrections count one hypothesis as several

R/three-step-categorical.R:355. `p.adjust` runs over every pair-by-category row. With a binary outcome the two rows of a pair are one test (`estimate` opposite in sign, same SE, same p). BH is unaffected, because duplicating every p leaves it unchanged, so the default is correct. Holm, Hochberg, Hommel, Bonferroni and BY are not: the family is doubled. With C > 2 categories each pair has C - 1 free differences, so counting C is still too many.

Evidence: 3 classes, binary outcome, `p_adjusted` from a factor against the same outcome coded 0/1 (three-step-supplement.R): BH 0.003199, 0.024695, 0.350354, identical. Holm 0.0064, 0.0659, 0.7007 against 0.0032, 0.0329, 0.3504. Hochberg and Hommel 0.0053, 0.0494, 0.3504 against 0.0032, 0.0329, 0.3504. Bonferroni 0.0064, 0.0988, 1 against 0.0032, 0.0494, 1. BY 0.0078, 0.0605, 0.8584 against 0.0059, 0.0453, 0.6423. Under Holm, Bonferroni and BY the middle pair crosses 0.05 by coding alone. A second data set in three-step-misc.R gives the same pattern.

Smallest fix: for C = 2 adjust over one category per pair and copy to its twin; for C > 2 document that the family counts every category. Changes `p_adjusted` under non-default `adjust` only. Probes: three-step-supplement.R, three-step-misc.R section 1.

### F4. MINOR, CONFIRMED (the mismatch); size varies. V1 ignores nesting while step three clusters

R/three-step.R:593. `vcov(x)` is the observed information. With one group class on nested data that likelihood treats observations as independent, while the same call's step three is cluster-robust by default. The correction therefore mixes a clustered V3 with an unclustered V1. It is not stated on the help page.

Evidence, G = 60, m = 15, one group class, between-group shifts in indicators and membership, SE by class:

| Seed | Fixed | Estimated, observed V1 | Estimated, robust V1 |
|---|---|---|---|
| 1 | 0.080, 0.125 | 0.127, 0.183 | 0.193, 0.266 |
| 2 | 0.099, 0.090 | 0.111, 0.103 | 0.111, 0.104 |
| 3 | 0.082, 0.084 | 0.086, 0.089 | 0.085, 0.089 |
| 4 | 0.080, 0.079 | 0.081, 0.081 | 0.082, 0.081 |
| 5 | 0.086, 0.083 | 0.087, 0.085 | 0.088, 0.086 |

Mostly under 1.5%; 50% on seed 1. Whether the robust V1 is the right reference on that seed is not settled by one data set; a small coverage cell would settle it. Smallest fix: use `vcov_type = "robust"` for V1 when `vcov_type = "cluster"` and the fit has more than one observation per group, or say on the page that V1 is model-based. Any Latent GOLD step-one kit case with nesting and one group class would move under a code change; I did not check whether the kit has one. Probe: three-step-v1-nesting.R.

### F5. MINOR, CONFIRMED. Modal overall test returns NA with no condition when a category is absent from a class

R/three-step-categorical.R:108. Any p exactly 0 or 1 returns an NA statistic, df and p silently. BCH outside the range warns and ML on the bound warns, but `p == 0` under `method = "modal"` passes `outside` (strict inequality) and has no message. The pairs on the same data report p = 6.8e-64 with class 1's SE of zero treated as known. Evidence: `three_step(fit, d, "rare", method = "modal", contrast = "overall")` gives `statistic = NA, df = 1, p_value = NA`, no warning. Smallest fix: warn `latents_boundary` there. Changes no number. Probe: three-step-misc.R section 2.

### F6. MINOR, CONFIRMED. Help page misstates the `variance` attribute

R/three-step.R (return section: "`variance` on a numeric ML result"); R/three-step-ml.R:543 sets it only for `contrast = "none"`. Evidence: `attr(., "variance")` is "varying" for none, NULL for pairs and overall. Same page: `adjust` "ignored for contrast = none" omits overall. Fix: set the attribute in the early returns, or narrow the sentence. Changes no number. Probe: three-step-v1-nesting.R, last line.

### F7. MINOR, CONFIRMED by output and derivation. The three-step chain's structure is undisclosed for general-engine step-one fits

R/lta-three-step.R:17-99, docs R/three-step.R "Transition fits". Two facts the page does not state:
(a) step three is always a homogeneous first-order chain. A step-one fit with `transitions = "occasion"` returns one intercept per origin and destination (table in three-step-lta-endtoend.R); `order = 2` follows the same code path (`.lta_designs(..., "homogeneous")`, not separately run).
(b) under `measurement = "occasion"`, one pooled D is applied to every occasion, though D then really differs by occasion, which breaks the fixed-error assumption. Under invariant measurement the pooled shares make D the same at every occasion by construction, so this concerns only this option.

Evidence for (b), separation falling over four occasions: D diagonal by occasion 1.00, 0.985, 0.915, 0.80; pooled 0.925. Slopes with pooled D 0.768, 0.836; hand ML with occasion D 0.787, 0.862, about 0.1 SE. A small effect on this case. Fix: say so on the page, or refuse `measurement = "occasion"` there. A code change would move only those fits. Probe: three-step-lta-occasion-measure.R.

### EDGE, listed for completeness

- `method = "modal", step_one = "estimated"` refuses on a fit without Wald inference (prior), though the page says modal is unaffected (R/three-step.R:593 precedes :599).
- Step-one SEs change by about 0.1% when one indicator is multiplied by 1e-4 and another by 1e4 (0.18955 vs 0.18942 at tol 1e-13); fixed SEs are identical to 7 digits and the likelihoods agree. Not traced further. three-step-mutants-5.R.
- A stayer `lta()` appears to carry more than one group class internally, so `r3step(initial_covariates =)` gives the "one group class" message rather than the stayer one; still refused correctly.
- A one-occasion-only data set reaches `moves[, 1L]` on NULL (R/lta-three-step.R:238-244), from reading only; not run.

## 3. Verified sound

| Item | How | Result |
|---|---|---|
| Three-step chain likelihood | Brute-force path sum in plain base R, K = 3, including a missing emission and unequal lengths | Per-person difference 8.9e-16 |
| Transition logits use staying as reference; initial against last | Independent softmax coding against `.lta_log_transition` | Identical to 6 decimals |
| Chain analytic scores | numDeriv of brute force; per person | 1.6e-9; 1.1e-9 |
| Chain end to end | Independent base-R forward ML from hand posteriors, modal and D, K = 2 | LL equal; estimates within 2.9e-7; observed and robust SEs equal to 6 digits |
| Occasion posteriors (`w_k f_k(y_it)`) and D | By hand, both engines | 9.4e-16; 5.6e-16; D diff 0 |
| Chain equivariance and row order | Units and origin of covariates; refit on shuffled rows | Slope ratios exact; LL diff 1e-13; estimates to 7 digits |
| Chain multiplicity | BH on covariate rows only | Equal to `p.adjust` by hand |
| r3step pooled, K = 3 | Independent likelihood; clustered sandwich by numDeriv | Estimates to 6 digits; SEs equal to 7 digits |
| r3step reference class | Label swap reverses sign, equal SE; attribute names last class | Exact |
| r3step group level with unsorted rows and character ids | Independent optim on `p$units`-keyed covariate; refit on shuffled rows | Same to 7 digits. `values[!duplicated()]` is safe because group ids number by first appearance |
| by_group_class | Independent two-level likelihood, numDeriv Hessian and CR1 sandwich | LL equal; SE within 1e-5; gradient vs FD 1.4e-7 (K = 3) |
| by_group_class transform | x2 times 100 plus 1e5 | LL diff 1e-13; slope ratio 0.01, SE ratio 0.01 |
| by_group_class with one group class | Against pooled | Same to 1e-5; intercept row is named `group_class_1` |
| Single-start optima | 10 random starts each, chain K = 3 (5 cases) and by_group K = 3, H = 2 or 3 (6 cases) | None above the reported LL by 1e-3 |
| Step-one first-order term | Independent VVV Cholesky chart and binary LCA, own likelihood, Hessian, Jacobian (BCH, proportional, ML) | SE and pair SE equal to 5 or 6 digits (ML: re-solved by EM at 1e-13) |
| Step-one with V1 = 0 | Patched closure | SE difference exactly 0, numeric and categorical, all methods |
| Sign of the step-one Jacobian | Derivation: dtheta3 = I33^-1 dS3 (R/three-step-ml.R:511-518, R/three-step-categorical.R:66-79) is right. | Correct; note the sign cannot affect any reported number, since only J V1 J' enters and no cross term |
| Omitted cross term | Stacked M-estimator on one data set per separation | Adds 9-22% of the increment; SE +1 to 2%. Bootstrap SD with step one refitted 0.069, 0.074 (entropy 0.57) against estimated 0.076, 0.071 and fixed 0.067, 0.063; 0.089, 0.090 (entropy 0.42) against 0.100, 0.089 and 0.077, 0.071. Single bootstraps of 300, not a coverage claim |
| ML numeric Hessian and scores | numDeriv at outcome units 1, 1e-4, 1e6, both variance models | Gradient 1e-7; Hessian 1.8e-7; identical across units |
| ML start selection | Nine bounded 1-5 outcomes, against a hand EM from the proportional start | Same solution every time |
| Categorical ML | Closed-form information vs numDeriv; sandwich SE; overall Wald on logits by hand, df 4 | 1.7e-7; 3.1e-8; 52.57795 both |
| r3step Hessian at step 1e-4 | numDeriv | 1.6e-7; gradient at BFGS stop 2e-7 |
| Outcome unit invariance with correction | 1e4 scale and 1e6 shift | SE ratio exactly 1e4 (BCH, ML) |
| Binary BCH pairs against 0/1 numeric BCH | Same data | Estimate, SE and p identical |
| Refusals | 16 fit kinds by 8 calls (three-step-refusals.R), plus direct effects (three-step-chain-refusal.R) | Weights, noise, covariates, additive, cross-level, occasion transitions, stayer all refused with their class. Count, ordinal, prior and direct effects run with sensible output. Only slip: F2 |

Line-by-line (matches the mathematics; earlier audit also ran these numerically): D[t, s] rows by true class; BCH row of inv(D) by assigned class; ratio influence w(y - mu)/T with raw negative weights; absolute weights only in the Kish thin-class count; CR1 applied before pseudo-rows; logit-scale interval half-width q SE / (p(1 - p)); pseudo-rows carrying covariance into pairs and overall.

## 4. Test gaps

Mutants patched into the namespace at run time, nothing written to R/; three-step-mutants-*.R.

| Formula broken | tests/testthat | Equivalence |
|---|---|---|
| Chain: D transposed in emissions | MISSED (test-lta-three-step, audit-1101) | Caught, test-latentgold-step3-markov |
| Chain: chain posteriors instead of occasion posteriors | MISSED | Caught, same |
| Chain: robust meat doubled | MISSED | MISSED: test-latentgold-step3-markov and tests/equivalence/test-lta-three-step |
| Chain: BH family includes intercepts | MISSED | MISSED: the same two files |
| Pooled r3step: meat not clustered | MISSED (test-three-step, clusters, audit-1101) | MISSED by tests/equivalence/test-tidy-three-step; test-jstats-mplus did not finish inside the 5-minute cap, undetermined |
| Categorical: logit CI half-width without 1/p(1-p) | MISSED (categorical, review-1101) | MISSED, test-latentgold-categorical-distal |
| Categorical overall: no 1/p delta | MISSED | Caught, same |
| Numeric ML: class sizes uniform rather than average posteriors | MISSED (ml, audit-1101) | Caught, tests/equivalence/test-three-step-ml |
| Pair SE assuming independence | Caught only by a categorical test; no numeric one in those files | not run |
| by_group_class: back-shift, robust meat, family | Caught, test-latents-controls | n/a |
| Step-one: sizes not updated; V1 zero | Caught, test-three-step-step-one | n/a |

Scope of this table: only the files named were run, not equivalence/run.R or the harness, whose `equivalence/harness/three-step-reference.R` may catch further mutants. "MISSED" means by those files. Mutant 4's two later lines were cut by the 5-minute limit, so they were rerun in three-step-mutants-6.R.

The chain's robust SE and its BH family: my probe (three-step-lta-endtoend.R) checks both by hand and they are right today; nothing in the named test files would notice if they broke.

Other gaps:
- `three_step(step_one = "estimated")` on `fit_staged()` is untested (F2).
- No test asserts any behaviour at a diverging logit in r3step (F1). `test-three-step-clusters` tests only the group count.
- The sign of the step-one Jacobian cannot be tested through the SEs.
