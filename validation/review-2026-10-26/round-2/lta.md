# Round-2 review: latent transitions (slug `lta`)

Reviewer: round-2 mathematical review, 2026-10-26. Environment: R 4.5.2, macOS,
`NOT_CRAN=true OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`, package loaded with
`pkgload::load_all()`. Oracles: an independent brute-force path enumeration (all
K^T state paths, written from scratch with `dnorm`/`mvtnorm::dmvnorm`/`dpois`,
Rabiner 1989 eq. 17), an independent log-space forward pass, `numDeriv`
(gradients, Hessians, per-group Jacobians) and Monte Carlo calibration against
chi-square. No file under `R/`, `tests/`, `equivalence/` or `man/` was edited.
The LTA files were being edited by another agent during the review
(`R/lta-engine.R` and `R/transitions.R` last changed 22:00); every finding below
was re-run on the tree as of 22:00 and still reproduces.

Scripts: `/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-lta/`
(`helper.R` holds the brute-force likelihood; `sNN_*.R` are the checks cited below).

## Scope (files)

R/transitions.R, R/transition-inference.R, R/lta-engine.R, R/lta-starts.R,
R/lta-inference.R, R/lta-methods.R, R/lta-three-step.R, R/lta-bootstrap.R,
R/lta-compare.R, R/sequences.R, R/tna.R; plus the new regression file
tests/testthat/test-review-1026-lta-open.R (LTA-F1 rank refusals, LTA-F2
uncounted trailing logits, LTA-F3 empty rows in sequences, LTA-F5 centred and
scaled fitting chart mapped back by J V J' and A' I A).

Topics: forward-backward recursion and scaling; initial/transition logit
M-steps (orthonormal basis of `.multilpa_weighted_logits()`); occasion-specific
vs homogeneous transitions and measurement; second order; mover-stayer; the
fitting chart and its Jacobian in inference; analytic scores vs numDeriv;
parameter counts for every option; three-step (ML, step-three chain);
parametric bootstrap LRT; nonparametric bootstrap; `compare_models()` for
transition fits; tna summaries; FIML, empty occasions, `drop_empty`, sampling
weights.

## Findings

| id | severity | file:line | what the code does | what it should do | reproduction | measured effect | proposed fix |
|---|---|---|---|---|---|---|---|
| R2-LTA-1 | medium | R/lta-engine.R:1130-1136 (`.lta_finalise_general`, `isFALSE(drop_empty)` branch) | Under `drop_empty = FALSE` a weighted `lta()` sets `bic <- ... log(groups$n)`, `bic_individual <- ... log(n_observations)` and `n_informative <- n_observations`, counting groups and rows of sampling weight zero. `nobs()`/`BIC()`/`compare_models()` use `.multilpa_group_sample_size()` = `sum(sampling_weights > 0)`, so the stored criterion and the methods disagree. This is the ML-1 defect (zero-weight units must not change the likelihood or the criteria) left on the general transition path. | Sample sizes over units of positive weight only (ML-1); `fit$bic == BIC(fit)`. | s13_compare_weighted.R (also s7_samplesize.R row 4) | 50 groups, 5 of weight 0: `fit$bic` 970.4047 vs `BIC(fit)` = `compare_models()` `bic_groups` 969.2457 = the same fit with the zero-weight groups deleted (969.2457). `bic_individual` 985.6539 vs 984.4949; `n_informative` 200 vs 180. `get_results(fit, "fit")` prints the wrong 970.4047. | In the branch use `n_groups_used <- .multilpa_group_sample_size(result)` and the rows of positive weight: `result$n_informative <- sum(positive_rows)`, `bic` with `log(n_groups_used)`, `bic_individual` with `log(sum(positive_rows))`. |
| R2-LTA-2 | medium (wrong refusal) | R/lta-inference.R:485 (`.lta_inference` packs every column), 254/350 (pack/scores use all `first_columns`); tables in R/lta-methods.R:20-75, 80-125 | LTA-F2 now leaves flat coordinates out of `n_parameters` (intercepts of moves after the last observed occasion; the measurement block of an all-missing occasion), but the packed vector that inference differences and inverts still contains them. The information then has exact zero eigenvalues and every Wald path (`parameter_inference()`, `vcov()`, the SE columns of `get_results()`) is refused with `latents_singular_information` for the whole fit. The flat values (start values) are printed as estimates without a flag. | Inference over the identified coordinates only; flat coordinates reported with `NA` SE and flagged (as the homogeneous engine flags empty transition rows with `estimated = FALSE`). | s2_counts.R case (a); s12_misc.R; s16_flat_block.R | (a) `transitions = "occasion"`, occasion 5 all missing for 20 of 80 groups: `n_parameters` 15, packed 17, numDeriv Hessian rank 15 (so the count is right), inference refused. `transition_coefficients` prints `occasion_5` logits -1.927218 and -1.564924 (the start values) with SE `NA`; `"transitions"` prints occasion-5 probabilities 0.873/0.127/0.173/0.827 as model-implied. (s16) `measurement = "occasion"`, same data: `n_parameters` 35, packed 43, inference refused; the `"profiles"` table prints start-value means and variances for occasion 5 (only a `latents_thin_profiles` warning). | Store `occupied` on the fit; in `.lta_inference` build `keep <-` the identified coordinates (`.lta_identified_columns()`, `.lta_identified_blocks()`), difference/invert `information[keep, keep]`, and fill the covariance with `NA` rows/columns elsewhere; add an `estimated` column to the coefficient, transition and profile tables. |
| R2-LTA-3 | low | R/lta-engine.R:1207-1216 (`.lta_identified_columns`), 936-940 (refusal) | With `order = 2` and no indicator observed after occasion 2, every second-order logit governs only moves into all-missing occasions, so the likelihood is flat in all of them; the flat rule applies only to `transitions = "occasion"` intercepts, so they are counted (and covariate columns of the second block are counted under either `transitions`). The fit is not refused. | Not counted, or refused as the first-order analogue is ("No occasion after the first has an observed indicator"). | s2_counts.R case (b) | `order = 2`, T = 3, occasion 3 all missing: `n_parameters` 15, Hessian rank 11 (the K^2 (K-1) = 4 second-order logits flat); inference refused as singular. | Refuse with `latents_bad_transition` when `order >= 2` and `max(which(occupied)) < 3`; otherwise in `.lta_identified_columns` drop the whole second block when no occasion from the third on is occupied. |
| R2-LTA-4 | low (inconsistent refusal) | R/transitions.R:858-866 (`.lta_setup`, homogeneous engine) | The default (homogeneous) engine accepts data with no indicator observed after occasion 1. The transition matrix is then flat (start values), yet it is reported with `estimated = TRUE`, expected counts and counted in `n_parameters`; only inference refuses. The general engine refuses the same data with `latents_bad_transition` (lta-engine.R:936). | Same refusal in both engines. | s2_counts.R case (c) | `lta()` default, all occasions after the first missing: `multilpa_transitions` with probabilities 0.971/0.029/0.271/0.729, `expected_count` 160.6/4.8/20.2/54.4, `estimated` TRUE, `n_parameters` 11. | Call `.lta_occupied_occasions()` in `.lta_setup` after the layout check and raise the same `latents_bad_transition`. |
| R2-LTA-5 | low | R/lta-engine.R:60-86 (`.lta_designs`), 260-285 (`.lta_check_transition_rank`) | The new rank refusal (LTA-F1) reads rows that do not enter the likelihood: the initial design over all groups (including weight-zero groups) and the transition design over every move inside a span (including weight-zero groups and moves into trailing all-missing occasions). A covariate that varies only there passes, its slope is unidentified, and it is counted. Same class of gap as ML-2 (fixed for indicators). | Judge rank on the rows the likelihood reads: groups of positive weight; moves into occasions up to the last occupied one. | s2_counts.R cases (d), (d2), (f) | (d) transition covariate varying only in 10 weight-zero groups: `n_parameters` 13, Hessian rank 11. (d2) initial covariate, same: 12 vs 11. (f) covariate varying only on moves into an all-missing trailing occasion: 13 vs 11. Coefficients are printed (e.g. `cz` 6.199, `cv` 0.062) with no SE. | Pass row masks to `.lta_check_design_rank()`: initial rows `sampling_weights > 0`; transition rows `within & positive weight & destination <= max(which(occupied))`. Use the same rows for `.lta_design_chart()`. |
| R2-LTA-6 | low | R/lta-engine.R:941-955 (warning), 1190-1200 (doc of `.lta_identified_columns`) | An all-missing occasion followed by observed ones under `transitions = "occasion"` warns (`latents_unidentified_occasion`) but counts all 2 K (K-1) logits of the moves into and out of it. Only their product A_t A_{t+1} is identified (K (K-1) free values), so K (K-1) counted parameters are unidentified and AIC/BIC are penalised for them. The roxygen says they "stay counted". | Count the identified dimension, or document the count as an upper bound. | s2_counts.R case (e) | K = 2, occasion 2 empty: `n_parameters` 15, Hessian rank 13 (2 = K (K-1) flat directions); inference refused as singular. | Subtract K (K-1) per mover class for each maximal run of inner empty occasions (each run of r empty occasions leaves r K (K-1) unidentified), or say in `?lta` and the warning that AIC/BIC count them. |

Not defects (recorded for the contract): (i) the `"transitions"` table's
`expected_count` and `get_tna()` use posterior two-slice counts, which include
model-implied moves into all-missing occasions; at the EM fixed point the
per-class matrices equal counts / row sums (8.7e-8), so the aggregate is the
documented one (`?get_tna`). (ii) `r3step()` with transition covariates fits a
homogeneous first-order step-three chain whatever the step-one transitions
(`order = 2`, `"occasion"`); this is the documented structural model. (iii)
At covariate unit 1e-4 and origin 1e8 a fit is 5.8e-5 below the reference
log likelihood: the data themselves keep only about 4 significant digits of
the covariate (eps x 1e8 = 2.2e-8 against values of size 1e-4), and the
estimates agree within 5.8e-5 SE.

## Verified correct

| algorithm | check | maximum discrepancy |
|---|---|---|
| General-engine likelihood with transition + initial covariates (origin 5, units 3 and 0.5), mapped back from the fitting chart | fit LL vs brute force at the reported raw coefficients; numDeriv gradient of brute LL at the estimate, in SE units (s1) | LL 0; gradient 9.2e-6 SE |
| Chart Jacobian J V J' (observed) and A' I A / S A (scores), first order | `vcov()` vs `solve(-numDeriv::hessian(brute LL))` in raw coordinates; inference group scores vs numDeriv per-group Jacobian (s1) | SE rel 1.4e-8; scores 4.1e-9 |
| Robust (CR1) covariance through the chart | vs H^-1 S'S H^-1 G/(G-1) from numDeriv (s1) | SE rel 2.8e-8 |
| Second order, homogeneous, covariates (chart block `transition_map[later, later]`) | as above (s1) | LL 0; SE rel 2.6e-7; scores 1.4e-8; robust 5.6e-7 |
| Second order + occasion transitions + transition and initial covariates (first block `occasion_2` only) | as above, interior fit, n = 300 (s1c) | LL 4.5e-13; SE rel 9.0e-7; scores 1.9e-8; robust 1.8e-6 |
| Mover-stayer + covariates | as above (s1c) | LL 2.3e-13; SE rel 3.3e-7; scores 8.4e-9; robust 4.6e-7 |
| Chart invariance at extreme units/origins for order-2 occasion and mover-stayer fits (unit 1e-4 to 1e4, origin -1e6 to 1e8) | coefficients vs the mapped reference M b, SEs vs sqrt(diag(M V M')) (s3) | coefficients 5.8e-5 SE, SE rel 3.0e-4 (both at unit 1e-4, origin 1e8, see note iii); all other settings <= 2.2e-5 |
| Homogeneous engine: full covariance + categorical item + partial FIML + `occasions = "grid"` gaps + unequal lengths | LL vs brute force with `mvtnorm` marginal densities; unconstrained `vcov()` vs numDeriv Hessian; analytic group scores vs numDeriv; robust (s4) | LL 0; SE rel 3.2e-7; scores 3.0e-8; robust 6.5e-7 |
| Homogeneous engine, two group classes | analytic scores vs numDeriv (inference correctly refused: probability at bound) (s4) | scores 1.1e-8 |
| General engine full covariance (EM-only path; covariance-chart scores with finite-difference Jacobian) + categorical + partial FIML + occasion transitions + covariate; diagonal counterpart | LL, observed/robust SEs, scores vs brute force / numDeriv (s5) | full: SE 1.8e-8, scores 7.0e-9, robust 3.6e-8; diagonal: SE 2.1e-8 |
| Poisson count and ordinal indicators in the chain, invariant and occasion-specific measurement, FIML counts | LL vs brute force (`dpois`), SE and scores vs numDeriv (s10) | LL 2.3e-13; SE 1.1e-7 / 2.0e-8; scores 2.1e-8 / 5.6e-9 |
| Weighted pseudo-ML with weight-zero groups and a covariate | LL = sum w_i l_i; robust = H^-1 (sum w_i^2 s_i s_i') H^-1 G/(G-1), H = sum w_i H_i, G = groups of positive weight; default `vcov_type` robust (s11) | LL 2.3e-13; robust SE 4.1e-6 |
| Unequal lengths with occasion transitions and a covariate (centred chart turns zero placeholders past a span into -centre/scale) | LL vs brute force; SE vs numDeriv (s12) | LL 1.1e-13; SE 1.4e-5 |
| Scaled forward recursion (homogeneous) on long sequences | 30 sequences x 300 occasions vs independent log-space forward; extreme emissions (variances 1e-4, LL -4.1e7, transitions 1e-9) (s14) | 3.6e-11 absolute; 2.2e-16 relative; posteriors finite |
| Logit M-step on the orthonormal basis (`.multilpa_weighted_logits`) | algebra: design beta = (design R^-1 sqrt(n)) (R beta / sqrt(n)), monotone line search; EM + quasi-Newton end at brute-force maxima (gradient in SE units, s1/s1c/s4/s5/s10) | <= 7.7e-5 SE at interior fits (9.4e-5 at a boundary fit) |
| Parameter counts (8 specs of round 1 plus: covariates, order 2 homogeneous and occasion, stayer + covariates, full covariance, count/ordinal, occasion measurement, trailing empty occasion) | `n_parameters` vs packed length and vs numerical rank of the brute-force Hessian (s1, s1c, s2a, s4, s5, s10) | exact, except the cases in R2-LTA-3, -5, -6 |
| `r3step()` transition/initial covariates | slope and SE invariance at units 1e3, 1e-3 and origin 1e6, 50; step-three LL vs brute force with D[k, m] emissions and missing emissions for empty occasions; observed SE vs numDeriv (s9) | slope rel 6.1e-7; SE rel 3.7e-8; LL 0; SE 6.9e-7 |
| `bootstrap_lrt()` for nested transition fits | p = (1 + #{T* >= T}) / (B + 1) recomputed; calibration under the null, B = 200: homogeneous vs transition covariate (df 2) mean 2.05, var 4.49, KS p 0.86; occasion vs occasion + initial covariate (df 1) mean 0.97, var 1.53, KS p 0.16; 200/200 valid (s6, s6b) | p-value exact |
| Nonparametric bootstrap with label alignment, K = 3, covariates | bootstrap SE / Wald robust SE, B = 50 (s15) | ratio 0.91 to 1.13 (sampling error ~10%); 50/50 valid |
| tna handover | `transition_counts` vs brute-force posterior expected moves (2 classes, unequal lengths, FIML); `get_tna()` weights = class-summed counts normalised; `inits` = sum_h pi_h initial_h; `get_group_tna()` = fitted matrices (s8) | 7.8e-14; 2.2e-16; 0; 1.1e-16 |
| Sample sizes: `nobs()`, `BIC()`, `fit$bic`, `compare_models()` `bic_groups` | 50 groups, 3 wholly empty, 5 weight-zero: drop_empty TRUE gives 42 (weighted) / 47 (unweighted, both engines); drop_empty FALSE unweighted 50 (s7, s13) | consistent except R2-LTA-1 |
| Trailing all-missing occasion count (LTA-F2 fix) | `n_parameters` 15 = Hessian rank 15 (s2a) | exact (inference defect is R2-LTA-2) |
| Regression files on the current tree | test-review-1026-lta-open.R 41 passes, test-review-1026-lta.R 35, test-lta-general.R 33, test-tna.R 62; 0 failures | — |

## Not checked (and why)

- LMest, seqHMM and Mplus were not re-run: the brute-force enumeration and
  numDeriv cover the same recursions and derivatives, and round 1 compared
  depmixS4.
- OPG covariance for `multilpa_lta` (only observed and robust were compared).
- Covariance structures other than VVV/diagonal in the general engine (EEE,
  VEI, ...), and K > 3; more than one mover class with covariates (the chart
  Jacobian was exercised with one mover class plus a stayer class, whose
  blocks follow the same packing).
- `search = "emEM"` and `cores > 1` (start orchestration in lta-starts.R, no
  mathematics); `.lta_enumerate()` criteria beyond what `compare_models()` reads.
- Three-step with weights or several group classes (both refused by design).
