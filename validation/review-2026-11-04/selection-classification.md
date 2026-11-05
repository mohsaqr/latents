# Review 2026-11-04: selection and classification

Package version 1.1.15, branch `features`. Probes are in `validation/review-2026-11-04/probes/selection-classification-*.R`. Nothing in R/, man/ or tests/ was changed. The mutation runs edited throwaway copies in the scratchpad.

Counts: 0 CRITICAL, 2 MAJOR, 7 MINOR, 4 EDGE.

## 1. What was reviewed

| File | Depth |
|---|---|
| diagnostics.R | Read in full. Criteria, entropy, AvePP and OCC hand-computed. |
| classification-entropy.R | Read in full. EN recomputed independently on plain and covariate two-level fits. |
| compare-lrt.R | Read in full. About 45 real pairs run, including transition fits. |
| compare-models.R | Read in full. Schwarz weights recomputed. |
| lmr.R | Read in full. One K=1 v K=2 case rebuilt with numDeriv, and a weighted one. `.latents_wynn` was read only, not separately exercised. |
| enumeration.R | The bootstrap LRT, simulation and criteria parts were read and run. `get_candidate` and the print helpers were read lightly. |
| selection.R | Deltas, knees and neighbours read and checked. The plot functions were not reviewed. |
| bivariate-residuals.R | Read in full. Gaussian, categorical and pooled statistics hand-computed. |
| diagnostics-checks.R, report.R | Read in full. Check wording compared with the numbers on real fits. |
| tree.R | Flow, share and matching calculations read and checked. Drawing geometry and colours were not reviewed. |
| align-labels.R, alignment-service.R | Read in full. Permutation invariance run on 10 specifications. |
| sensitivity.R | Profile and transition paths read; profile path run once. Growth path read only. |
| candidates.R, enumeration-type.R, grid-views.R | Partly read, not independently probed. I relied on the 2026-11-01 audit for these. |

Not reviewed: `.additive_bootstrap_lrt`, `.lta_simulate` internals, plot rendering, parallel code. I did not rerun the BLRT or LRT size studies; I only read their CSVs.

## 2. Findings

### F1. MAJOR, CONFIRMED: "estimates at a bound" reports "no / ok" while a parameter is on its bound

Where: `R/diagnostics-checks.R:68-75` reads `x$boundary`; the same flag drives `R/compare-lrt.R:48`, `R/enumeration.R:294` (grid column) and the trust view in `R/selection.R`.

What is wrong: `x$boundary` comes from `.multilpa_covariance_boundary()` (`R/gaussian-moments.R:175`), which looks only at Gaussian variances and eigenvalues. General `lta()` fits add transition and initial probabilities (`R/lta-engine.R:1516`). Response probabilities are not looked at, and neither are negative-binomial dispersions. The check's own text says "a variance or probability sits at its lower bound", and `?diagnostics` promises `"check"` for active parameter bounds.

Evidence (`selection-classification-bound-lrt.R`, 2-class LCA, n = 300):
- Class 1, item `a`: probabilities 1.0000000 and 3.855e-09. `any(p <= .multilpa_probability_floor(fit))` is TRUE; `fit$boundary` is FALSE.
- `diagnostics()` prints "estimates at a bound: no, ok" and "No check needs attention."
- `parameter_inference()` on the same fit raises `latents_boundary_fit`: "Wald inference is unavailable for bound-active categorical response probabilities".
- `compare_models(free, joint, test = "lrt")`, with both fits at the floor, returns LR 2.231, df 2, p 0.328 and raises no warning. The `latents_boundary` warning exists to warn about exactly this reference problem.

Further evidence:
- Bundled data (`selection-classification-bound-realistic.R`): `enumerate_lca(student_esm, ...)` with 4 items. `classes2` and `classes4` have probabilities at the floor and `boundary` FALSE, while every LMR row is `NA` with a note about bound-active probabilities. The grid contradicts itself. The 4- and 5-class fits on these 4 items also warn that they are not identified; `classes2` does not.
- Same flag on a simulated LCA grid: `classes2` and `classes3` at the floor, `boundary` FALSE.
- Homogeneous `lta()` (`selection-classification-checks.R`): response probabilities 1 and 0, `boundary` FALSE.
- Negative binomial: dispersion 1e-08, fit warns "This is a boundary fit", `boundary` FALSE, check says "no, ok". The `NB equal < NB varying` LRT (`selection-classification-lrt.R`) raised no LRT warning on such a fit (LR 0.0026, df 1, p 0.959).

Who is affected: anyone fitting LCA or mixed models, where a response probability at 0 or 1 is common. They get a clean diagnostic and a chi-square p-value with no caution.

Smallest fix: add a separate predicate covering response probabilities, using `.multilpa_probability_floor()`, and count parameters at zero. Use it in the check, the LRT warning and the trust view. Do not widen `fit$boundary` itself: `bootstrap_lrt()` refuses fits with `isTRUE(model$boundary)`, so widening it would newly refuse LCA fits.

Does the fix change reported numbers: no estimate, criterion or p-value changes. Statuses and warnings change; if the grid `boundary` column is also changed, its values change. Owner's decision.

### F2. MAJOR, CONFIRMED: `bootstrap_lrt()` fails every replicate for joint item pairs and non-syntactic categorical names

Where: `R/enumeration.R:1102`, `as.data.frame(drawn, stringsAsFactors = FALSE)`.

What is wrong: `as.data.frame()` applies `check.names`, so `a:b` becomes `a.b`, and `result[, object$vars]` then fails with "undefined columns selected". `.latents_blrt_replicates()` catches this as a failed replicate.

Evidence (`selection-classification-blrt-joint.R`):
- Categorical `local_dependence = "a:b"`, 2x1 against 2x2: statistic 25.72, `p_value` NA, 0 of 3 valid, each replicate's error "undefined columns selected".
- LCA with items named `item one` and `item-two`: statistic 14.71, `p_value` NA, 0 of 3 valid.
- Continuous non-syntactic names work (`score 1`, `score 2`: statistic 14.08, p 0.25, 3 valid). Only categorical blocks fail.

The warning reads "a refit failed, did not converge or reversed the likelihood", which points the analyst at estimation rather than at a bug.

Who is affected: every categorical `local_dependence` fit, whose joint name always contains `:`, and any categorical indicator with a non-syntactic name. `?bootstrap_lrt` lists local dependence among what the models share, so this is a supported use. `.lta_simulate` was not tested for the same problem.

Smallest fix: add `check.names = FALSE` at that line. Tested in memory only (`selection-classification-more.R`):
- 5 of 5 replicates valid, p 0.1667.
- Simulated joint margins 0.182, 0.299, 0.331, 0.187 against model 0.180, 0.305, 0.330, 0.185, over 40 simulated data sets.
- Half-observed joint labels had zero simulated frequency.
- A `missing = "fiml"` joint fit also ran: 3 of 3 valid.

Does the fix change reported numbers: only `NA` p-values become numbers.

### F3. MINOR, CONFIRMED: a regular nested pair of transition models is refused

Where: `R/lta-compare.R:225`, the extra `rank(a$variance_model, ...)` condition in `.lta_nested()`.

EII stores `variance_model = "varying"`, while EEI and EEE store `"equal"`, so the condition fails although the letter rule passes. Of 61 pairs nested by the letters, 59 agree. The two refused are EII within EEI and EII within EEE (`selection-classification-lta-nesting.R`). `compare_models(test = "lrt")` raises `latents_bad_nesting` and `bootstrap_lrt()` raises `latents_bad_argument`, both saying the null is not nested. The same pair in `lpa()` is accepted (LR 0.799, df 2).

Smallest fix: drop that condition, since the letter rule already covers it. No existing number changes; two refusals become results.

### F4. MINOR, CONFIRMED: the local-independence detail says "largest" but names the smallest adjusted p

Where: `R/diagnostics-checks.R:388-389`.

Evidence (`selection-classification-more.R`): the detail reads "largest: y1 x y2 in profile_1, residual 0.240", while the table's largest residual is 0.351 in profile_2 (n 69). Profile_1 has the smaller adjusted p (8.5e-11 against 2.9e-03).

Smallest fix: change the wording to "smallest adjusted p", or pick by `abs(residual)`. No number changes.

### F5. MINOR, CONFIRMED: "modal vs estimated proportions" compares unlike quantities under sampling weights

Where: `R/diagnostics-checks.R:290`. The modal share is an unweighted count, the estimated share is weighted.

Evidence (`selection-classification-class.R`): perfectly separated classes (relative entropy 1.000, AvePP 1) with weights 1 and 5 give "largest gap 0.333". The gap is entirely due to the weights, not classification.

Smallest fix: weight the modal share in this one check. This changes that displayed value on weighted fits only. Status stays `info`.

### F6. MINOR, CONFIRMED by reading: the BLRT print misstates what makes a replicate invalid

Where: `R/enumeration.R:1732-1734` says invalid replicates include a spurious class "reaching a boundary". By `R/enumeration.R:1672-1673`, validity reads only convergence and the reversal window, so a boundary refit stays valid. The summary print at the same file states it correctly. Text fix only.

### F7. MINOR, CONFIRMED: general `lta()` fits have no residuals, and the docs do not list them

Where: `?get_residuals` and `?diagnostics` list regression, growth, additive, cross-level and noise fits as lacking the table. A homogeneous `lta()` fit returns it. The same data with `transitions = "occasion"` raises `latents_bad_argument`, and `diagnostics()` says "not defined for this model family" (`selection-classification-more.R`). Doc fix; no number changes.

### F8. MINOR, CONFIRMED by reading: the `compare_models()` print mixes conventions under plain labels

Where: `R/compare-models.R:447-449`. "BIC" is `bic_groups` and "ICL" is `icl_individual`, in one row of a two-level comparison, with nothing in the print saying so. The two penalties use different `n`. Label fix; no number changes.

### F9. MINOR, CONFIRMED: ordinal within nominal is refused, as documented

A nominal indicator contains the adjacent-category ordinal model with scores 1 to C, and the restriction is interior, so this is a regular nested pair. It is refused because the coding differs (`selection-classification-lrt.R`). The help says coding must agree, so this is conservative rather than wrong. Listed for the owner to decide; no action needed if the rule is intended.

### EDGE

- `drop_empty = FALSE`: `multilpa()`'s help says the criteria "count the informative rows either way". In the probe, individual `n` stays 115 of 120 rows, but group `n` is 30 against 29, and the entropy table uses 120 rows (entropy sum 40.118 against 36.839). By reading `.multilpa_complete_entropy()`, ICL's EN follows those same rows; I did not print the ICL difference.
- A one-indicator two-level LPA reports residuals as "not defined for this model family"; the true reason is that there is no pair (`R/report.R`, default `residual_reason`).
- `sensitivity()`: `best` and `optimum` include unconverged refits (`R/sensitivity.R:201`), while its plot counts only converged ones. Read only, not run.
- PLAUSIBLE: the help calls the categorical statistic "Latent GOLD's bivariate residual times its df". From my recollection LG's BVR is pooled, which would match `by = "overall"` and not the default `by = "profile"`. An LG listing would settle it.

### Checked and not a defect

- A correctly specified homogeneous LTA was flagged (r = -0.242, adjusted p 0.014). The true-state sample correlation in that data set was -0.248, so it was chance. Over 25 data sets the rate of minimum p below 0.05 was 0.08, the same as the LPA control (`selection-classification-lta-bvr.R`). This was a small check, not a size study.
- The categorical BVR is conservative, as documented. The repo's `validation/bivariate-residual-size.csv` shows per-test rates of 0.0007 to 0.0068 at the 5% level and mean statistics of 0.17 to 0.52 on df 1 (not rerun). It keeps power against strong dependence: one data set with real within-class dependence gave statistic 137.8 and was flagged (`selection-classification-more.R`).

## 3. Verified sound

**Criteria** (`selection-classification-ic.R`). All 14 columns equal a base-R hand computation to 8 decimals on a VVI 2x2 fit (log likelihood -777.9456, q = 11, n 40 and 240). Single-level conventions coincide.

**Complete-data entropy.** The hand value equals the package's on the same fit through those criteria, and on a covariate two-level fit (33.12952 both ways, log likelihood -777.935). The derivation matches the docstring.

**Entropy.** Relative entropy by hand: 0.8407129 for individuals, 0.6905905 for groups, both equal.

**Weights against duplicated groups.** Log likelihood scales by 82/40 (-782.283 against -1603.680). Relative entropy 0.8014155 both ways, groups 0.6337489 and 0.6337515. AvePP 0.9457514 and 0.9457515.

**AvePP, OCC, modal counts** (`selection-classification-class.R`). Equal to hand values at both levels. `subject_profiles` and `group_classes` equal `max.col(first)`.

**Error matrix and BCH.** The matrix equals `P(assigned | true)` by hand. BCH weights differ from `D^-1[assigned, class]` by 0. Summed weights equal posterior totals (119.9964, 180.0036). The direction also follows from the derivation.

**Bivariate residuals.** Gaussian r, Fisher z and p, categorical Pearson (0.01232495), df, and the pooled statistic (0.09931699) all equal hand values.

**Tree.** Total mass is 1 and margins match shares within 4.4e-16, including weights with FIML. Mass equals the posterior cross-product.

**Deltas, knees, Schwarz.** Deltas read correctly, Schwarz weights differ by 0. The knee on a rising curve returned `NA` by its documented rule (hand distances 0, 0.653, 0.474, 0.228, 0 were not reported); the distance formula itself is calibrated in `test-review-1026-selection.R:171`.

**LRT nesting** (`selection-classification-lrt.R`).
- Correctly accepted, with hand-checked df where stated: covariates (hand LR 0.5254386, df 1), slopes, direct effects, pair covariances, joint items (df 2), and transition pairs (occasion transitions df 4, occasion measurement df 24, second order df 4, covariate df 2).
- Correctly refused: reversed order, direct effects on other indicators, direct effect against membership covariate, VVI pair against EEE, Poisson within negative binomial (boundary), cross-level fits, different centering, spherical within EVE.
- I found no non-nested or boundary-count pair that gets through. The letter rule and the singular list hold by derivation at generic nulls. A probability or dispersion on its bound does get through, without a warning (F1).

**Bootstrap LRT** (`selection-classification-tests.R`).
- p = (1 + b) / (V + 1) by hand: 0.05 at 19 of 19 valid.
- Refits use each model's own controls (read).
- The refit rule reproduces q and log likelihood (difference 0) on VEI, EVE, pairs, direct effects, ordinal, negative binomial, grand centering and grand scaling.
- Negative-binomial simulated moments: mean 2.87 against 2.93, variance 11.99 against 12.05 (`selection-classification-blrt-joint.R`).

**VLMR and aLMR** (`selection-classification-tests.R`). W rebuilt from numDeriv scores and Hessians, independently of the package. Eigenvalues 1.027, 0.592, 0.125, -0.0048, 0, 0. Mean 1.739854 and sd 1.68602 equal the package. Monte Carlo p 0.38562 against 0.3846056, aLMR 0.42077 against 0.4199949, adjusted LR 1.507476 equal. Weighted case with B on squared weights: mean 2.00593 and sd 1.968854 equal. Imhof at chi-square(3), q = 5: 0.1717971 equals `pchisq`.

**Label alignment** (`selection-classification-align.R`). Log likelihood difference 0 on 7 plain specifications including ordinal and negative binomial. Likelihood reproduced on 3 covariate specifications including slopes by group class. Modal labels, criteria, classification, BVR and BCH are invariant. A person-shuffled refit aligned to the reference. Noise class stays 0.

## 4. Test gaps

Method: one-line mutations on a throwaway copy, run against the named test files only. A survivor below means it survived those files, not the whole suite. A wider rerun of survivors against more files had produced no results when this report was written, so those files are unconfirmed.

Survived the files run:

| Mutation | Files run, all passing |
|---|---|
| Schwarz `exp(-d/2)` to `exp(-d)` | test-enumeration, test-compare-lrt, test-review-1026-selection. `bic_weight` also appears in test-audit-1101-additive-compare, test-growth-mixture, test-review-1026-codex, which were not confirmed. |
| Fisher `sqrt(n-3)` to `sqrt(n)` | test-bivariate-residuals, test-review-1026-selection and -r2 |
| Tree flows ignoring weights | test-tree, test-review-1026-selection and -r2. test-tree's weighted test (line 218) checks shares, not weighted mass. |
| Simulation ignoring the group class | test-bootstrap-lrt-extended, test-enumeration, test-direct-effects |
| Pair-complete total in the categorical BVR expected table | test-bivariate-residuals |

Caught: SABIC, CAIC, AWE, CLC (only test-review-1026-selection-r2), OCC (only test-review-1027-diagnostics), relative entropy, EN group part, BCH direction, error matrix, aLMR `d`, Vuong sign, BLRT plus-one, knee (only test-review-1026-selection), categorical df, LRT df, logit rebasing.

KIC survived the files I ran, but `test-release-items.R:31` asserts its hand value. That is from reading; the mutant was not run against that file.

Missing tests tied to findings:
- No `bootstrap_lrt()` test on a categorical joint pair or a non-syntactic categorical name (F2).
- No test that the bound check or the LRT warning fires on a response probability at its floor (F1).
