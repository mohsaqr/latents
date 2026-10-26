# Round 2 mathematical review: selection and three-step

Reviewer: round-2 agent, area `selection`, 2026-10-26. Everything here was run
under `pkgload::load_all()` with `NOT_CRAN=true OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1` on the working tree as of this date. The
reproduction scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-selection/`
(written below as `<s>/`). No file under `R/`, `tests/`, `equivalence/` or
`man/` was edited.

## Scope (files)

R/enumeration.R (grid criteria, minima, AvePP summary, BLRT simulation and
p-value), R/enumeration-type.R, R/candidates.R, R/compare-models.R,
R/selection.R (deltas, knee), R/diagnostics.R (criteria, entropy, AvePP, OCC,
level weights), R/descriptives.R (ICC), R/lmr.R (Vuong W, Imhof, aLMR),
R/three-step.R (error matrix, BCH, proportional, modal, r3step pooled and by
group class), R/three-step-ml.R (ML step three), R/tree.R (class-split flows
and matching), R/grid-views.R. I also read the helpers these call in
R/noise.R (informative rows and groups, group sample size) and R/weights.R.

Round 1 (`area-reports/selection.md`) is not repeated here: F1 (BFGS
tryCatch) and F2 (row names) are fixed; F3 (soft against hard ICL), F4 and F5
are still as reported. F5 was fixed in `.multilpa_r3step_frame()` but not in
its by-group-class twin (S-6 below).

## Findings

| Id | Severity | File:line | What the code does | What it should do (reference) | Reproduction | Measured effect | Proposed fix |
|---|---|---|---|---|---|---|---|
| S-1 | **high** (single-level fits); medium (two-level) | R/diagnostics.R:74-93 | The group-convention `icl`, `awe` and `clc` charge only the entropy of the **group-class** posteriors. A fit with one group class, which includes every single-level `lpa()`/`lca()`, has group entropy 0. So `icl_groups` equals `bic_groups`, `awe_groups` is the AWE penalty alone and `clc_groups` equals the deviance. On a single-level fit the two conventions have the same `n` and the same units, yet `icl_groups` and `icl_individual` differ. `summary(grid)$best` then reports `clc_groups` as choosing the largest candidate, because the deviance always falls. | ICL, AWE and CLC charge the entropy of the complete-data classification (Biernacki, Celeux & Govaert 2000; Banfield & Raftery 1993; Biernacki & Govaert 1997). For a two-level model the latent variables are the group class C_j **and** the profiles X_ij, so the entropy is H(C_j \| y) + sum_h P(C_j = h \| y) sum_i H(X_ij \| C_j = h, y_ij). This is the individual EN when there is one group class. On a single-level fit the group and individual versions must coincide. | `<s>/ml_criteria.R`, `<s>/group_entropy_grid.R`, `<s>/lpa_grid_best.R` | Single-level `enumerate_lpa(1:5)`: `icl_groups` = BIC exactly (4036.796 for 2 profiles, against `icl_individual` 4326.109); `clc_groups` falls monotonically from 4158.853 to 3957.975, and `best` lists `clc groups -> profiles5`, so "4 distinct candidates are minimal under some criterion" counts this degenerate column. Two-level, 2 profiles x 2 group classes, 60 groups: EN(groups) 8.91, marginal EN(individuals) 120.57, joint complete-data EN 123.90. The `icl_groups` penalty is 230.0 too small. The plot (`.multilpa_distinct_criteria`, R/plot.R:449) already hides the group versions when every candidate has one group class, but the tables, `get_results(fit, "information_criteria")` and `summary()$best` do not. | Charge the joint entropy above under the groups convention. It is computable from the stored fit: the stored individual posteriors equal sum_h P(h \| y) P(k \| h, y) to 4.4e-16, so the conditionals are available. With one group class it reduces to the individual EN, which makes the single-level conventions identical. If that is not wanted, return `NA` for the three group-convention columns whenever `n_group_classes == 1`, as the documentation already does for an undefined group entropy, and leave them out of `best`. |
| S-2 | medium | R/compare-models.R:541-542 (`.latents_comparison_signature`) | Categorical items enter the data signature as character labels under `categorical`, and ordinal items as ordered factors under `extra`. Two fits of the same items, one nominal and one ordinal (`ordinal =`), are therefore refused as "fitted to different data". | The adjacent-category ordinal model is nested in the nominal one on identical data. Comparing them by BIC is the standard check (Latent GOLD offers both). Both likelihoods are multinomial probabilities of the same responses. | `<s>/compare_ordinal.R` | `lca(d, items, 2)`: logLik -2570.105, q = 25. With `ordinal = items`: logLik -2571.382, q = 17. `compare_models()` raises `latents_incomparable_models`. | Canonicalize every discrete item (nominal or ordinal) to `as.character()` of its labels, item by item, in one signature slot. Keep count and Gaussian items separate, since a density and a probability are not comparable. |
| S-3 | medium | R/diagnostics.R:281-300 (`.multilpa_classification_summary`) | Under sampling weights `odds_correct_classification` divides the odds of the **unweighted** AvePP by the odds of the **weighted** class share. The ratio mixes two populations. | OCC = [AvePP/(1 - AvePP)] / [pi/(1 - pi)] (Nagin 2005) compares classification odds with the prior odds of the **same** population. Use both weighted (the population convention the shares already use) or both unweighted. | `<s>/weighted_entropy.R` | n = 400, class 1 over-represented by weight 4: class 1 OCC reported **2.37**. Unweighted throughout it is 5.80, weighted throughout 3.37. The reported value is below the conventional 5 cut-off that the documentation quotes, while the sample-consistent value is above it. | Compute AvePP with the row weights (`weighted.mean` within each modal class) whenever the shares are weighted. Alternatively report OCC with the unweighted share. |
| S-4 | low/medium | R/diagnostics.R:74-76, 16-22; R/three-step.R:100-112 | Under sampling weights the log likelihood is the weighted pseudo log likelihood. The entropy charge of ICL, AWE and CLC, relative entropy, AvePP and the classification-error matrix (`get_results(fit, "classification_errors")`, `"bch_weights"`) are unweighted. `three_step()`/`r3step()` refuse weighted fits, so only the displayed tables are affected. | The entropy term of ICL, AWE and CLC is part of the complete-data **pseudo** log likelihood, sum_i w_i sum_k z_ik log(pi_k f_k), so it should carry the same w_i (Skinner 1989 pseudo-likelihood). The error matrix of the population is D[t,s] = sum w p_it 1(W=s) / sum w p_it. I did not establish what Mplus or Latent GOLD print, so the choice should be documented whichever way it is made. | `<s>/weighted_entropy.R` | EN unweighted 115.37 against weighted 93.24, so `icl_individual` moves 2800.456 to 2756.208 (44.2). Relative entropy 0.584 against 0.664. D[1,2] 0.102 against 0.043, D[2,1] 0.155 against 0.284. | Weight the entropy charge by the unit weights (they already sum to the number of units, so the scale is unchanged). Weight D, or document that D and AvePP describe the sample, not the population. |
| S-5 | low | R/diagnostics.R:31, 80; R/get-data.R:518 | `kic = -2L + 3(q + 1)` with q the total parameter count, citing Cavanaugh (1999). | Cavanaugh's general form is -2L + 3k with k the dimension of theta. The "+1" appears only in his regression form, where it counts sigma^2 beside k regressors. The (q + 1) form is the tidyLPA convention. | algebra | Every KIC is 3 larger than Cavanaugh's general form. A constant offset changes no ranking or difference. | Either use 3q, or keep 3(q + 1) and say in the documentation that it follows tidyLPA's convention. |
| S-6 | low | R/three-step.R:1177 (`.multilpa_r3step_group_frame`) | `sqrt(pmax(diag(covariance), 0))` turns a negative variance into SE 0, then statistic and p-value `NA`, with no warning. This is round-1 F5, which was fixed in `.multilpa_r3step_frame` but not in this twin used by `by_group_class = TRUE`. | As in `.multilpa_r3step_frame`: `NA` SE plus `latents_singular_information` warning. | by reading. Not triggered: on a boundary fixture (`<s>/r3step_check.R` with its first settings, G = 80, m = 8) a group-class intercept ran to 17.2 (independent optimum 35.4, same log likelihood -400.1597, a flat direction) and was reported with SE 242.7 and no warning, but the diagonal stayed positive. | 0 on the fixtures tried | Reuse the `.multilpa_r3step_frame` rule. Consider also a boundary warning when an intercept's \|estimate\| exceeds ~15 with a flat likelihood (Mplus prints such parameters as fixed). |
| S-7 | low (convention) | R/diagnostics.R:85-93; R/get-data.R:521 | (Round-1 F3, still present, and widened.) ICL and AWE use the **soft** entropy, while the cited Biernacki et al. (2000) ICL and Banfield & Raftery (1993) AWE use the hard (MAP) classification likelihood. The documented formulas are the soft ones, so the code matches its documentation but not its references. | Cite McLachlan & Peel (2000, section 6.10.3) for the soft ICL-BIC, or switch to the hard entropy -sum_i log t_{i, MAP}. | round 1 | round 1: soft -1891.17 against hard -1854.13 (mclust sign) | Correct the reference, or compute the hard form. |

## Verified correct

| Algorithm | Check | Max discrepancy |
|---|---|---|
| AIC, KIC, BIC, SABIC, CAIC, AWE, ICL, CLC under both conventions on a **two-level** fit with 2 group classes (n = 600 rows, 60 groups). Groups convention uses n = groups and group EN; individuals convention uses n = rows and individual EN, as documented. | hand formulas (`<s>/ml_criteria.R`) | 0 |
| Consistency of stored posteriors: individual posterior = sum_h P(h \| y_j) P(k \| h, y_ij), from the stored means, variances and profile mix | independent recomputation | 4.4e-16 |
| Zero-weight units: no effect on log likelihood, any criterion, entropy or classification table (ML-1 holds for this area) | refit with 6 extra zero-weight groups (`<s>/zero_weights_ic.R`) | 0 (all columns identical) |
| r3step pooled (Vermunt 2010; Asparouhov & Muthen 2014): estimates, observed-information SE, group-clustered CR1 sandwich. Covariate at origin 50 with SD 10, plus a group-level covariate. | independent likelihood sum log sum_t pi_t(z) D[t, W], BFGS -> Nelder-Mead -> BFGS, numDeriv Hessian and per-unit gradients (`<s>/r3step_check.R`) | est 4.3e-7; SE rel 2.2e-8 (observed), 4.0e-8 (robust) |
| r3step `by_group_class = TRUE`: two-level likelihood sum_j log sum_h eta_h prod_i sum_k pi_hk(z) D[k, a_i], the alpha/beta/gamma packing and the back-transform from the orthonormal chart | independent likelihood, three optimisers, numDeriv Hessian and per-group scores | log likelihood identical (-947.2649); est 4.2e-7; SE rel 6.9e-7 (observed), 1.4e-6 (robust), at the numDeriv limit |
| Classification-error matrix at `level = "groups"` (Vermunt 2010) | hand D[t,s] = sum_j h_jt 1(W_j = s) / sum_j h_jt (`<s>/bch_groups.R`) | 0 |
| BCH at `level = "groups"` (Bolck, Croon & Hagenaars 2004): weights D^-1[W_j, ], means, CR1 SEs over groups, pairwise contrast SE | hand | est 1.3e-15, SE 2.8e-17, contrast SE 5.6e-17 |
| ML step three at `level = "groups"` (pi fixed at mean posterior, D fixed) | independent likelihood + numDeriv sandwich with G/(G-1) | means 3.0e-8, variances rel 1.3e-7, SE rel 6.1e-9 |
| ICC(1) constancy: exact constants (0.1, 0.3, 0.7, 1/3, 2/3, 0.2, 123.456, 1e8+0.1) in 7 unbalanced groups give `NA`; group-constant values give exactly 1 | direct call | exact |
| Kneedle / chord distance (Satopaa et al. 2011) | read: the vertical distance below the first-to-last chord in min-max coordinates is proportional to the perpendicular distance, so the argmax is the same as Kneedle's difference curve for a falling curve. Round-1 numeric check of 1.4e-14 stands. | by derivation |
| BLRT p-value (1 + #{T* >= T})/(B + 1), withheld when any replicate fails; Monte Carlo SE sqrt(p(1 - p)/(B + 1)) | read (R/simulation-service.R:211) | exact |
| BLRT simulation: group class drawn from `group_probabilities`, profile from the class's mix, observed missingness mask carried for FIML, covariates held for covariate fits | read | n/a |
| aLMR = LR / (1 + 1/(d log n)) with d the parameter difference (Lo, Mendell & Rubin 2001, eq. 15; Mplus TECH11) | read; `equivalence/test-mplus-lmr.R` arithmetic | exact |
| VLMR/aLMR (Vuong 1989 Thm 3.3; Lo, Mendell & Rubin 2001) for a **bivariate** VVI LPA, 2 against 3 profiles, n = 300 (round 1 checked a univariate one) | independent W = [-B1 A1^-1, -B10 A0^-1; B10' A1^-1, B0 A0^-1] from numDeriv per-unit scores and Hessians in a (logit, mean, log-sd) chart, then CompQuadForm::imhof (`<s>/vlmr_bivariate.R`) | LR 42.064125, VLMR p 3.6375e-06, aLMR p 5.7018e-06, reference mean 5.697498, sd 3.709558: identical in every printed digit. The log likelihood at the package estimates equals the package's (difference 0) |
| Tree flows: sum_i w_i p_i(a) p_i(b) / sum w over the same informative units for every row; exact permutation matching up to 8 classes | read | n/a |
| Unit alignment for `level = "groups"`: groups are numbered in order of first appearance (`.multilpa_prepare_groups`), which is the order of `group_posteriors`, `split(outcome, group_index)` and the r3step group design (`!duplicated(group_index)`) | read | n/a |

## Not checked (and why)

- **VLMR for a categorical LCA.** `<s>/vlmr_lca.R` (independent Vuong W in a
  logit chart) ran on two designs (n = 400, 2-class truth; n = 1500, 3-class
  truth). Both gave a 3-class fit with a bound-active response probability,
  for which the package returns `NA` with the note "Wald inference is
  unavailable for bound-active categorical response probabilities", so no
  categorical p-value was compared. The independent unconstrained optimum was
  7.2e-5 above the package's log likelihood, the cost of the 1e-10 bound. The
  Gaussian check above and round 1's univariate one are the VLMR evidence.
- **What Mplus and Latent GOLD print for entropy, AvePP and the error matrix
  under sampling weights.** Neither program was run, so S-4 is argued from the
  pseudo-likelihood and not from a reference output.
- **The multilevel three-step assumption.** The error matrix used by
  `r3step(by_group_class = TRUE)` and the BCH and ML methods assumes
  P(W \| X, C) = P(W \| X). Under a two-level model, W_ij is computed from
  posteriors that use the whole group's data, so this is an approximation.
  I did not measure its bias by simulation.
- **BLRT reference distribution and size.** Not simulated: only the p-value
  arithmetic and the simulation code were read.
- **mixture_regression's ICL and BIC inside `compare_models()`.** These are
  computed in the regression files, which belong to another area.
- **Plots** (selection views, tree drawing): no visual output was changed or
  rendered.

## Fix status (batch 1: S-2, S-5, S-6)

| Id | Change | Test (tests/testthat/test-review-1026-selection-r2.R) | Fails on the old code? |
|---|---|---|---|
| S-2 | `R/compare-models.R`: the data signature holds one character-label vector per discrete item (nominal or ordinal), in item-name order, in a `discrete` slot, and counts in a `count` slot (`.latents_comparison_responses()`). | "nominal and ordinal fits of the same items are compared"; guards: a changed response is still refused; a Gaussian against an ordinal treatment is refused | Yes: the first test errors on the old code (`latents_incomparable_models`). The two guards pass on both versions, as intended. |
| S-5 | `R/get-data.R` (Information criteria section and references) and `vignettes/articles/model-selection.Rmd`: KIC `-2L + 3(q + 1)` is described as Cavanaugh (1999) in the form tidyLPA computes (Rosenberg et al. 2018), with the constant-3 difference from `-2L + 3q` stated. The number is unchanged. | none (documentation only) | n/a |
| S-6 | `R/three-step.R`: new `.multilpa_r3step_standard_errors()`, used by both `.multilpa_r3step_frame()` and `.multilpa_r3step_group_frame()`: NA plus `latents_singular_information` for a negative variance. It also no longer leaks base R's "NaNs produced" warning. | "a negative variance in the group-class r3step table is NA, with a warning" | Yes: 3 failures on the old code. |

The "fails on old code" checks were run on a copy of R/ in the scratchpad with only these two edits reverted.
S-1, S-3, S-4 and S-7 wait for the diagnostics files.

## Fix status (batch 2: S-1, S-3, S-4, S-7)

| Id | Change | Effect re-measured with the reproduction scripts |
|---|---|---|
| S-1 | New `R/classification-entropy.R`: `.latents_row_conditional_entropy()` (each row's E_C[H(X_i \| C, y)] from the E-step's `joint`), `.multilpa_complete_entropy()`. `.multilpa_finalise()` (`R/multilpa-starts.R`) and the covariate finalise (`R/covariates.R`) store `classification_entropy_rows`. `.multilpa_information_criteria()` (`R/diagnostics.R`) charges the complete-data entropy under both conventions. It falls back to the row posterior entropy with one group class, and to `NA` for older fits with several group classes. Transition fits keep the per-level entropies (option B). The definitions now name the entropy. | Single-level `enumerate_lpa(1:5)`: `icl`, `awe` and `clc` are identical under both conventions, and `best` names one candidate per criterion (CLC → profiles1, no longer profiles5). Two-level fit: `icl_groups` 4209.4871 equals the independent joint-entropy value. |
| S-3 | `.multilpa_average_posterior_matrix(p, weights)` (weighted AvePP) feeds the classification table, the average-posterior table, the enumeration and comparison AvePP summaries, and the quality columns of non-profile grids. | Weighted fixture: class 1 OCC 3.3667, against 2.3746 before; equal to the weighted-throughout value. |
| S-4 | Under weights the `"entropy"` table (sum and relative entropy) is weighted, through `.multilpa_relative_entropy(p, weights)`. So are the grid's and `compare_models()`'s `profile_entropy`/`group_entropy`, through `.multilpa_level_relative_entropy()`. `.multilpa_level_assignments()` carries `weights` and `.multilpa_error_matrix()` weights D, so the BCH weights follow. Modal counts stay counts. | See the tests below. |
| S-7 | Documented convention: soft entropy (McLachlan & Peel 2000, section 6.10.3), not the hard form; stated in `?get_results` and in the printed definitions. | n/a |

Documentation changed: the *Information criteria* and *Classification quality* sections of `R/get-data.R`, and its `"classification"`, `"average_posteriors"`, `"classification_errors"`, `"bch_weights"` and `"entropy"` entries; `vignettes/articles/model-selection.Rmd` ("Two conventions"). man/ was not regenerated.

The new tests in `test-review-1026-selection-r2.R` (8 added, 12 in total, 40 expectations) all fail on the code before batch 2 (a copy of R/ with the diagnostics and error-matrix edits reverted). The existing expectations that encoded the per-level convention were updated: `test-diagnostics.R` (ICL/AWE) and the EN-1 test in `test-review-1026-multilevel-r2.R`.

Not changed, and needing a decision:
- `.multilpa_distinct_criteria()` (R/plot.R:441) still drops the `_groups` panels of AWE, ICL and CLC on grids with one group class. Its comment ("ICL over groups equals BIC over groups") is no longer true, but changing it alters plots, so it was left for a rendered review.
- The average-posterior plot (`R/plot-gg.R:1205`) is unweighted.
- `R/lta-compare.R:515`, the transition grid's entropy column, is unweighted (LTA file).
