# Round 2 review: multilevel (two-level LPA/LCA, membership covariates, weights, search, cross-level)

Date: 2026-10-26. Reviewer: round-2 agent, area **multilevel**. No file under `R/`, `tests/`,
`equivalence/` or `man/` was edited. Every number below was produced in this session with
`NOT_CRAN=true OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript`, under
`pkgload::load_all()`, `cores = 1`. The scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-multilevel/`
(below, `<S>/` stands for that folder). `common.R` holds the shared simulator.

## Scope (files)

Reviewed: `R/multilpa.R`, `R/multilpa-starts.R`, `R/multilpa-pool.R`, `R/start-pool.R`,
`R/lpa.R`, `R/multilca.R`, `R/covariates.R`, `R/covariate-inference.R`, `R/cross-level.R`,
`R/weights.R`, `R/search.R`.

I also read the helpers they call: `R/kernel-structures.R` (nested and logit E-steps,
`.multilpa_weigh`), `R/kernel-em.R` (`.latents_em`, `.latents_run_starts`, the emEM screen),
`R/noise.R` (informative rows and groups, the group sample size), `R/diagnostics.R` (criteria,
entropy, classification tables), `R/centering.R`, `R/inference-service.R` (sandwich, CR1, OPG)
and `R/additive.R` (`.additive_prepare`, used by cross-level).

Round 1 is already fixed and was not re-reported: ML-1 (zero-weight normalization in
`multilpa()`) and ML-2 (constant-indicator check over zero-weight rows). I re-verified both:
they hold for `drop_empty = TRUE`. Findings ZW-3, ZW-4 and CL-3 below are what ML-1 did not
reach.

## Findings

| id | sev. | file:line | what the code does | what it should do | reproduction | measured effect | proposed fix |
|---|---|---|---|---|---|---|---|
| CL-1 | **high** (silently wrong log likelihood for any weighted restricted cross-level fit with a zero-weight group) | `R/cross-level.R:317` (`individual_scale`) | Puts the rating factor's pseudo log likelihood back on the group weight scale with `sum(row_weights) / length(row_weights)`. The rating fit (`multilpa(id = NULL, weights =)`) normalizes its row weights to sum to the number of **positive-weight** rows (`R/weights.R:56`, after ML-1), so the factor that undoes it is `sum(row_weights) / sum(row_weights > 0)`. | The pseudo likelihood is `sum_j w_j log L_j` (Skinner 1989; the package's own definition). Restricted factorization: `log L_j = log sum_h omega_h f(m_j \| h) + sum_i log sum_k pi_k f(x_ij \| k)`. | `<S>/w1_crosslevel_zero.R`, `<S>/w1b_restricted_fix.R` (24 groups with weights 1 to 3, plus 5 zero-weight groups) | Unpadded LL −568.7007, which equals a direct computation of `sum_j w_j log L_j` at the reported estimates to 1e-8. Padded LL **−478.8841**, but the direct computation gives −568.7007. The error is +89.8 log-likelihood units, or `(N+/N_all − 1)` times the rating part. BIC 1194.61 becomes 1018.38. Every estimate is unchanged. | Line 317: `sum(row_weights) / sum(row_weights > 0)`. Verified on paper in w1b: the corrected formula gives −568.70065072 for both data sets. Line 413 (weights of the thin-profile diagnostic) needs the same change. |
| CL-2 | **medium** (wrong estimate; needs weights + a zero-weight group + `between_variance = "equal"` + full family) | `R/cross-level.R:48` (`.cross_level_between_update`) | The shared between-group variance divides the posterior-weighted scatter by `nrow(group_means)`, which counts every group. Under weights the posteriors are the weighted ones, whose total is the number of positive-weight groups. | M-step `sigma^2 = sum_h sum_j tau_jh (m_j − mu_h)^2 / sum_h sum_j tau_jh`, with weighted `tau` (pseudo-ML). | `<S>/w1_crosslevel_zero.R`; fix check in `<S>/w1c_crosslevel_patch.R` | Adding 5 zero-weight groups to 24 shrinks the group-mean variances by factors of 0.8247 and 0.8118 (expected 1). LL falls from −554.8704 to −555.3547, so the fit is no longer the maximum. Profile means move by 2.1e-3. Without zero weights the two denominators agree, so unweighted fits are unaffected. | Divide by `sum(counts)` instead of `nrow(group_means)`. With that patch applied in memory, the padded and unpadded fits agree exactly (LL −554.87042867 in both, max variance difference 0). |
| CL-3 | medium (criteria) | `R/cross-level.R:384-385, 551, 556` | `bic` uses `log(stats$n_groups)`, `bic_individual` uses `log(nrow(x))`, and `logLik()`/`nobs()` report `n_groups`. All of these count zero-weight groups and their rows. ML-1 fixed this for `multilpa()` but not for the cross-level families. | Count the units of positive weight. `.additive_prepare()` already returns `n_used_groups` and `n_used_obs`. | `<S>/w1_crosslevel_zero.R` | Full family, varying: LL identical (−554.796583), but BIC 1169.98 becomes 1173.57, `bic_individual` 1204.02 becomes 1207.62, and `nobs` 24 becomes 29. A zero-weight group can therefore change the class count BIC selects (ML-1's own argument). | Use `stats$n_used_groups` and `stats$n_used_obs` in `bic`, `bic_individual`, `logLik(nobs =)`, `nobs()` and the `fit` table's `n_groups`/`n_obs` (or add `n_used_*` columns). |
| ZW-3 | medium (wrong refusal) | `R/covariates.R:814-819` (`.multilpa_cov_designs`); the same call is reached from `R/covariate-inference.R:473` | The covariate design rebuild calls `.multilpa_center_indicators()` and `.multilpa_check_centered()` without `carrying`. `multilpa()` passes `carrying` (`R/multilpa.R:1734-1739`), so a zero-weight group with a degenerate person scale is excused in the base fit but refused in the covariate fit. | Zero-weight units influence nothing (the ML-1/ML-2 contract), and `R/centering.R` already implements the excuse. | `<S>/w7_person_scaling_zero.R` | `centering = "person", scaling = "person"`, one zero-weight group with a constant `y1`. The covariate-free fit succeeds (LL −349.2115). Adding `profile_covariates = "z"` raises `latents_bad_data` ("divides by a standard deviation that is zero ... in 1 group(s)"). Dropping that group by hand fits (LL −348.7927). | Pass `carrying = row_weights[first_rows] > 0` (or the group weights) to both calls in `.multilpa_cov_designs()`, and the stored weights from `.multilpa_cov_covariance()`. |
| ZW-4 | medium for spherical structures, low otherwise | `R/centering.R:145-150` (grand scale "judged on every row") and `:174` (grand offset), as used by `multilpa(weights =, centering/scaling = "grand")` | The grand standard deviation and the grand mean are taken over every row, zero-weight rows included, and unweighted. | Zero-weight units contribute nothing (`R/weights.R` header; ML-1). Grand scaling changes the solution for `EII`/`VII`, as `?multilpa` says. | `<S>/w6_scaling_zero.R`, `<S>/w6b_grand_scaling_eii.R` | Six zero-weight groups added, `scaling = "grand"`. For VVI only the units change: raw-unit means are identical, and LL moves 154.7 (Jacobian-corrected difference 2e-13), so BIC is not comparable between the padded and unpadded data. For **EII** the solution changes: profile means in raw units move by up to 5.98, group posteriors by up to 0.70, and the Jacobian-corrected LL by −25.1. With `centering = "grand"` alone, the reported (centred) means move by 0.99 while LL is unchanged. The covariate path behaves the same (LL difference 223 under grand scaling). | Compute the grand offset and scale over carrying rows. A weighted mean and SD would be more coherent still, but that needs a decision. At minimum, document that zero-weight rows set the grand scale. |
| ZW-5 | low | `R/noise.R:182` (`.multilpa_group_sample_size`), `R/multilpa-starts.R:177-189`, `R/noise.R:172` / `:209` (informative rows and groups) | Under `drop_empty = FALSE`, zero-weight groups and rows count in the group sample size (`bic`, `nobs`), in the group and individual entropy (ICL, AWE, CLC, relative entropy), and in `effective_*` row selection. `.multilpa_inference_groups()` (`R/inference-service.R:202`) excludes them "under either convention". | One convention for zero-weight units. The ML-1 rationale says "in no sample size, criterion or entropy". | `<S>/w2_zero_conventions.R` | `drop_empty = FALSE`, weighted, +5 zero-weight groups: LL is unchanged, but BIC moves 1072.37 → 1074.45, `nobs` 24 → 29, ICL(individuals) +8.42, CLC(groups) +0.44, AWE(groups) +4.60. The same holds for covariate fits. | In `.multilpa_informative_rows()`/`_groups()` and `.multilpa_group_sample_size()`, apply the positive-weight mask before the `drop_empty = FALSE` early return. |
| CV-1 | low | `R/covariates.R:1033` | Under `drop_empty = FALSE`, a covariate fit's `bic_individual` uses `nrow(data)`. `multilpa()` and the criteria table (`R/diagnostics.R:73`) use `n_informative` under both conventions (documented in `?multilpa`, `drop_empty`: "except in the individual sample size of a multilpa() fit, which excludes such rows either way"). | `n_informative`. | `<S>/w3_cov_bic_individual.R` | Three all-missing rows, FIML, `drop_empty = FALSE`: `fit$bic_individual` = 1059.0308 (log 144), while `get_results(fit, "information_criteria")` gives 1058.7782 (n = 141). The covariate-free fit agrees with its table (1058.7124). With zero weights, the field also counts zero-weight rows (1093.126 vs 1090.855, `<S>/w2_zero_conventions.R`). | Use `base$n_informative %||% nrow(data)` whatever `drop_empty` is. |
| EN-1 | low (convention; needs a decision) | `R/diagnostics.R:75-76`, `:5` | For a weighted fit, the deviance in ICL, AWE and CLC is the weighted pseudo deviance, but the entropy term `EN` is the unweighted sum over units. The criterion then has an unweighted classification term next to a weighted likelihood term. | `EN = sum_j w_j e_j`, with normalized weights, is the entropy of the pseudo classification likelihood. It equals the replicated-data entropy times the normalizing constant. | `<S>/w5_replication_entropy.R` | Integer weights 1 to 4 on 30 groups, with a covariate fit: unweighted individual EN 13.5216, weighted 11.3403, replicated data × c 11.3403 (exactly equal). Group EN 3.9509 vs 4.4123. ICL(individuals) would move by −4.36. Estimates match the replicated fit (coefficients to 3e-6, LL ratio = mean weight 2.733333). | Weight the entropy sums by `.multilpa_level_weights()` (already used by the classification table), or document the choice in `?get_results`. |
| ST-1 | low (documentation / determinism) | `R/multilpa.R:288-295` (group-class part of start 1); doc `R/multilpa.R:839-842` | The doc says "The first start is deterministic: Ward's ...". The profile part is. With `n_group_classes > 1`, though, the group-class probabilities of start 1 come from `stats::kmeans` on the group compositions with random centres, so start 1 depends on the seed. Start 1 is the protected emEM finalist and the base of the covariate model's first start. | Either a deterministic group-class start (for example Ward on the compositions, as `.cross_level_initialize()` already does for groups) or a doc that says so. | `<S>/s4_ward_start_seed.R` | Start-1 log likelihood at `max_iter = 0`, seeds 1 to 6: −810.876, −810.876, −810.876, −809.996, −812.417, −811.247 (3 group classes). With one group class it is −828.0655 for all six seeds. | Use `.multilpa_ward_assignments(group_composition, n_types)` for start 1, or correct the roxygen. |
| ST-2 | low | `R/covariates.R:308-318` (base fit inherits `search`), `:356-373` | Under `search = "emEM"` the covariate base fit (`n_starts = 1`) screens a pool of 5 and draws its random starts first. The covariate pool is therefore not "the kmeans search's own starts, drawn in the same order" (`?multilpa`, `search`). For the profile model this claim holds exactly. | Same random draws as the kmeans search for pool positions ≤ `n_starts`. | `<S>/s6_cov_search.R`, `<S>/s6b_cov_rng.R` | `n_screen = n_starts = 5`, `max_iter = 1`: start log likelihoods differ at starts 2, 4 and 5 (for example −757.6046 vs −752.6689). Start 1 is identical. Final maxima differ by ≤ 4.1e-6. | Run the base fit with `search = "kmeans"` (it is one Ward start either way), or note the exception. |
| RP-1 | low | `R/covariates.R:988-1106`, `R/cross-level.R:261-286` | Covariate and cross-level fits have no `n_best_replicated`, so enumerations and comparisons report `NA` for them. | The replication count every other family reports. | `<S>/s5_misc.R` | `fit$n_best_replicated` is `NULL` for a covariate fit and for a full cross-level fit. | Compute it from `starts$log_likelihood` (final stage, valid) with the shared `1e-6 (1 + \|LL\|)` rule. |
| DOC-1 | low (cosmetic) | `R/multilpa.R:1215-1216` (`@param weights`); `R/weights.R:118` | The doc says the weights are "scaled to sum to the number of units"; since ML-1 they sum to the number of **positive-weight** units. The Kish line prints "of 30 units" when only 26 carry weight. | Accurate text. | `<S>/s5_misc.R` | "Kish effective n 21.9 of 30 units" with 26 positive-weight units. | Reword the roxygen. Print `sum(weights > 0)`, or "30 units (26 weighted)". |

There is no defect in the core likelihood, posteriors, scores, information or sandwich
(see the next table).

## Verified correct

| Algorithm | Independent check | Max discrepancy |
|---|---|---|
| Two-level covariate likelihood: profile logits with group-class intercepts **and group-class slopes** (`profile_slopes = "group_class"`), group logits with a group covariate, mixed Gaussian + categorical indicators, FIML missing values in both kinds, sampling weights | Brute-force enumeration of all (h, z_1..z_4) configurations per group (30 groups, K = 2, H = 2), from the reported coefficients (`<S>/b1_bruteforce_cov.R`) | weighted LL 0; group LL 3.6e-15; group posteriors 4.4e-16; subject posteriors 4.4e-16 |
| Per-group covariate scores: categorical response logits, FIML, group-class slopes, group covariates, weights including zero-weight groups, full and equal-full covariance with FIML (log-Cholesky) | `numDeriv::jacobian` of `w_j log L_j(theta)` at a perturbed theta (`<S>/s1_scores_numderiv.R`) | 1.2e-8, 5.3e-9, 2.7e-9, 6.7e-9, 1.5e-8 (max \|score\| 4 to 10) |
| Observed-information vcov of the covariate fit | `solve(-numDeriv::hessian(LL))` | 2.1e-8 (categorical+FIML), 6.0e-7 (full+FIML); relative SE ≤ 3.3e-8 |
| Robust (CR1) vcov, unweighted and weighted, including zero-weight groups (G = positive-weight groups) | numDeriv H⁻¹ · G/(G−1) Σ s_j s_jᵀ · H⁻¹ | ≤ 4.4e-6 absolute; relative SE ≤ 8.8e-8 |
| Robust vcov against **`sandwich::vcovCL`** | Separated profiles (posteriors 0/1 to 1e-140), one group class, so the profile-logit block is a clustered logistic regression: `glm(binomial)` + `vcovCL(cluster = ~g, type = "HC0", cadjust = TRUE)`, unweighted and with integer weights (`<S>/s2_vcovCL.R`) | unweighted 3.7e-9 (SEs 0.1699247, 0.1339014 both); observed SE = `vcov(glm)` to 1e-7; weighted 2.2e-9 with the reference on positive-weight rows |
| — a note on that reference | `vcovCL()` with zero-weight rows is biased: `bread.glm` scales by the non-zero-weight count and the meat by all rows, which shrinks the variance by (n_ok/n)². Measured SE ratio 0.9664 against the predicted sqrt(290/300)² = 0.9667. The package was right; the reference needs zero-weight rows removed. | — |
| Zero-weight invariance (`drop_empty = TRUE`): EEI, EEE, VVV, VEI with FIML; mixed categorical + ordinal + count; covariate model with equal full covariance; covariate model with categorical and group-class slopes | +6 zero-weight groups with distorted values (`<S>/w4_zero_invariance.R`) | LL, BIC, `bic_individual`, nobs: 0. Parameters ≤ 4.4e-16. Robust SEs ≤ 5.3e-12 relative |
| Integer weights = replicated groups (covariate model) | Each group repeated w_j times (`<S>/w5_replication_entropy.R`) | coefficients 3.1e-6 (flat EM stop at tol 1e-13), means 3.2e-8; LL ratio = mean raw weight exactly (2.733333) |
| Restricted cross-level factorization `L_j = [sum_h omega_h f(m_j\|h)] prod_i sum_k pi_k f(x_ij\|k)` (no zero weights) | Direct weighted sum at the reported estimates (`<S>/w1b_restricted_fix.R`) | LL 0 (−568.70065072 both) |
| Full cross-level likelihood, group posteriors and closed-form between-level M-step (varying/equal × varying/equal) | Independent R implementation (`<S>/c1_crosslevel_full.R`) | LL 0; group posteriors ≤ 7.2e-16; between means/variances satisfy the M-step to ≤ 2.4e-7 (EM stop) |
| Cross-level parameter count | Hand count: K·d + K·d + H(K−1) + H·d + H·d + H − 1, with "equal" variants | 19/17/17/15 = hand count |
| emEM: a resumed start takes exactly the steps it would take unscreened | `n_screen = n_starts` emEM against the kmeans search: EEI, VVI, EEE × SQUAREM/plain (`<S>/s3_search.R`) | every start's LL difference 0; iterations identical; selected start and means identical |
| emEM finalists = protected Ward start + top screened (ties by pool order) | Recomputed from the starts table, n_starts = 5, n_screen = 20 | equal in all 6 settings |
| `n_best_replicated` | Recomputed: final-stage valid starts within 1e-6(1+\|LL\|) | equal in all 6 settings (5, 5, 1, 1, 5, 5) |
| Criteria under `drop_empty = TRUE` for covariate fits: BIC on informative groups, `bic_individual` on informative rows, entropy over informative rows/groups | Criteria table vs fit fields vs covariate-free fit (`<S>/w2_zero_conventions.R`, `<S>/w3_cov_bic_individual.R`) | identical (for example 1058.7782 = 1058.7782, n = 141) |

## Not checked (and why)

- **`cores > 1` and the shared enumeration pool** (`R/multilpa-pool.R`, `R/start-pool.R`). I read
  them but did not run them. Under `load_all()`, socket workers lack the development namespace,
  so parallel parity needs a fresh install (see LEARNINGS 2026-10-26). These files are plumbing
  over `.latents_select_starts()` / `.latents_select_screened()`, which I verified serially.
- **Comparisons with Mplus, Latent GOLD or multilevLCA.** I did not re-run them. The harness
  suites cover them, and round 1 reported them.
- **The covariance structures other than EEI/VVI/EEE/VVV inside the two-level model.** These are
  kernel files outside this area. I checked zero-weight invariance only for VEI.
- **`lpa.R`, `lca()`, `multilca.R` and `enumerate_lpa/lca`.** I read them; they are argument
  forwarding only, with no formula to check.
- **Severity in standard errors for CL-2 and ZW-4.** Cross-level fits refuse inference
  (`latents_unsupported_inference`), so the effect is reported in the reported quantity:
  variance ratio, log likelihood and raw-unit means.
- **The full `devtools::test()` suite.** I did not run it. This was a review only, and no code
  changed.

## Fixes (applied after the coordinator's decisions)

All twelve findings are fixed. Each has a block in
`tests/testthat/test-review-1026-multilevel-r2.R` (12 blocks, 44 expectations).

**How I checked that each test fails on the old code.** I copied the tree, restored the nine
original R files (backups in `<S>/orig/`) and ran the new file there. Every one of the 12
blocks failed: 1, 2, 4, (error), 10, 4, 1, 4, 1, 2, 4 and 1 failed expectations in file order.
I also checked ST-2 on its own: reverting only the one line of the ST-2 fix fails that block
alone (2 failures) and passes the other 11.

| Id | Change |
|---|---|
| CL-1 | `R/cross-level.R`: `individual_scale = sum(w)/sum(w > 0)`. The thin-profile diagnostic's restricted row weights are rescaled the same way. |
| CL-2 | `R/cross-level.R` `.cross_level_between_update`: the equal between variance divides by `sum(counts)`. |
| CL-3 | `R/cross-level.R`: `bic`, `bic_individual`, `logLik(nobs =)` and `nobs()` use `n_used_groups` / `n_used_obs`, which the fit now stores. `n_groups` and `n_obs` still count every unit, so the `fit` table's columns are unchanged. |
| ZW-3 | `R/covariates.R` `.multilpa_cov_designs`: passes `carrying` (groups of positive weight) to `.multilpa_center_indicators()` and `.multilpa_check_centered()`. `R/covariate-inference.R` rebuilds the designs with the fit's weights. |
| ZW-4 | `R/centering.R`: the grand offset (`.multilpa_center_only`, new `carrying` argument) and the grand scale (`.multilpa_indicator_scales`) are taken over rows of positive weight. Every row is still transformed. |
| ZW-5 | `R/noise.R`: new `.multilpa_carrying_rows()` / `.multilpa_carrying_groups()`. Under `drop_empty = FALSE`, informative rows and groups and the group sample size exclude zero-weight units. `R/multilpa-starts.R`: the fit's own group sample size does the same. The `drop_empty` doc gains the sentence. |
| CV-1 | `R/covariates.R`: `bic_individual` uses `base$n_informative` under both conventions. |
| EN-1 | `R/diagnostics.R`: `.multilpa_entropy_sum(probabilities, weights = NULL)`. The ICL, AWE and CLC entropy charge is the sampling-weighted sum (`.multilpa_level_weights()`). The `entropy` table stays unweighted. The convention is documented under "Information criteria" in `?get_results` (`R/get-data.R`, roxygen only). |
| ST-1 | `R/multilpa.R`: new `.multilpa_group_class_seed()`, Ward's clustering of the group compositions (with a deterministic tilt when compositions are fewer than classes). It is used for start 1 of every multilpa, transitions and cross-level fit with more than one group class. The `n_starts` doc says so. |
| ST-2 | `R/covariates.R`: the covariate base fit runs `search = "kmeans"`. It is one Ward start either way, so the result is the same, and the pool now draws in the documented order. Code fix; the doc is now true. |
| RP-1 | Covariate fits record `n_best_replicated` and `replication_tolerance` by `multilpa()`'s rule. The full cross-level fit records the count of final starts within `1e-6 (1 + \|L\|)`. The restricted family records the minimum of its two factors' counts. |
| DOC-1 | `R/multilpa.R` `weights` roxygen: "units of positive weight"; a zero-weight unit changes nothing, including the grand offset and scale. `R/weights.R`: the Kish line reads "of 26 units (4 more of weight zero)". |

**Re-running the reproductions** (`<S>/w1*`, `w2`, `w3`, `w5`, `w6*`, `w7`, `s4`, `s5`, `s6`):

- Padded and unpadded fits now agree exactly in LL, BIC, `bic_individual`, nobs and estimates,
  for every cross-level case, for `drop_empty = FALSE`, and for grand or person centring and
  scaling (EII included).
- Start 1 gives −809.9957634 for all 6 seeds.
- The covariate emEM and kmeans searches agree start by start, iterations included.

`man/` has not been regenerated. The roxygen changes in `R/multilpa.R` and `R/get-data.R` need
`devtools::document()`.
