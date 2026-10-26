# Review: latent transition / hidden Markov models (slug `lta`)

Reviewer partition: R/lta-engine.R, R/lta-starts.R, R/lta-inference.R, R/lta-bootstrap.R,
R/lta-compare.R, R/lta-three-step.R, R/lta-methods.R, R/transitions.R,
R/transition-inference.R, R/sequences.R, R/tna.R. Base: HEAD 8ac3576 plus the lead's
uncommitted constancy fix in R/lta-three-step.R. Environment: R 4.5.2, macOS,
OMP_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1, cores = 1. depmixS4, LMest, seqHMM,
numDeriv were installed; depmixS4 and numDeriv were used.

**No production file was edited by this reviewer.** One test file was added. The
mutation check below edited R/lta-three-step.R for a moment and then restored it from a
byte copy. `git diff --stat R/lta-three-step.R` afterwards shows only the lead's
6+/1- change.

## Scope

1. Codex's diff for my files (6 of 11 changed: lta-compare.R, lta-methods.R,
   lta-three-step.R, tna.R, transition-inference.R, transitions.R).
2. The forward-backward recursions (scaled one for homogeneous chains, log-space one for
   general chains), padding to the longest span and the missing-emission convention,
   the initial/transition M-steps (closed-form and multinomial-logit), occasion-varying,
   covariate, second-order and mover-stayer chains, occasion-specific measurement,
   sampling weights, Wald SEs (observed, robust CR1), relabelling for the bootstrap,
   simulation for the bootstrap LRT, LRT nesting rule, parameter counts, nobs/BIC,
   three-step (Bolck, Croon & Hagenaars 2004; Vermunt 2010; Di Mari, Oberski & Vermunt
   2016) step-two error matrix and step-three ML with its scores and sandwich.

References for the checks: Rabiner (1989) for the forward-backward recursions and the
xi/gamma moments; Collins & Lanza (2010) and Bartolucci, Farcomeni & Pennoni (2013) for
LTA/latent-Markov parameterisation and parameter counts; Vermunt, Tran & Magidson (2008)
and Bolck et al. (2004) for the three-step correction.

## Verdict on each Codex change

| Change | Verdict | Evidence |
|---|---|---|
| lta-three-step.R `.lta_three_step_inputs`: occasions with no observed indicator get modal = NA, log emission 0 (missing emission), D computed from observed occasions only | **Correct** | D matches the hand formula D[k,j] = sum_i p_ik 1(m_i=j) / sum_i p_ik over observed rows (max abs diff 2.2e-16). Step-three LL matches a brute-force path sum with NA-modal occasions given emission 1 (max abs diff 1.3e-15). Appending 11 trailing all-missing occasions changed r3step estimates by at most 4.4e-6 and SEs by at most 1.2e-4. Those gaps match the step-one refit tolerance: D moved 1.8e-7. |
| `.lta_occasion_posteriors`: pooled shares from informative rows only; general branch computes per-block densities with the default row offset | **Correct** | Equals the independent p(k|y) proportional to w_k f_k(y), with w = colMeans of the step-one posteriors over observed rows (max abs diff 3.3e-16). The row offset is constant within a row, so the posteriors are unchanged. The mutation that drops the `observed` restriction is caught by 3 failures in the new test file. |
| `informative_groups` restricts the objective, scores and the CR1 cluster count | **Correct** | A wholly missing person has LL 0 and score 0 (checked: max abs score 2.4e-10 from numDeriv). With 2 wholly missing persons out of 80, the r3step robust SEs match an independent sandwich H^-1 S'S H^-1 G/(G-1) with G = 78 (max rel diff 7.8e-8). Using G = 80 gives a rel diff of 1.6e-4, so the restriction has a real effect and is the right one. Observed-information SEs match numDeriv to 4.3e-8. Estimates match an independent BFGS on the brute-force likelihood to 1.2e-7. |
| Start values no longer solved through the basis (`start <- natural start`) | **Correct** | With every slope at zero, beta = T theta gives beta_0 = theta_0 and beta_slope = 0, so the two charts agree. A covariate at a 1e8 offset fits with the same slope as the centred one (-0.328421 in both). |
| `.lta_three_step_basis` spread check `<= 0` (Codex) | **Regression**, already fixed by the lead | The lead's exact range test refuses `rep(0.1, n)` with `latents_bad_covariate` for both `transition_covariates` and `initial_covariates` (rechecked through `r3step()`). test-review-1026-constancy.R passes (5/5). I agree with the fix. `any(spread <= 0)` is now redundant but harmless. I found no other tolerance-based constancy test in my partition: the only other check, `length(unique(v)) == 1L` for within-person initial covariates (lta-three-step.R:56), is already exact. A related but different gap is F1 below. |
| lta-compare.R `.lta_nested` also requires the same weights column, equal normalised sampling weights, `min_variance`, `min_probability` | **Correct** (a stricter refusal, not a regression) | A likelihood-ratio statistic between fits with different bounds or weights compares different objectives. Checked: equal-variance vs varying-variance fits with the same bounds are nested (TRUE); with `min_variance = 1e-4` on the alternative they are refused, and `bootstrap_lrt()` raises `latents_bad_argument`. |
| `.lta_check_data` also checks the weights column and compares designs with tolerance 0 | **Correct** | Designs and centred indicators are recomputed from the same data the same way, so exact equality is the right test. The `tryCatch` handler re-raises a classed refusal and does not swallow the error. |
| lta-methods.R / transitions.R `logLik()` and `nobs()` use `.multilpa_group_sample_size()` | **Correct** | With 2 wholly missing persons out of 50: `BIC(fit) == fit$bic` for both engines (1027.163 homogeneous, 1039.285 general) with nobs 48 under `drop_empty = TRUE`, and nobs 50 (1027.612 / 1039.897) under `drop_empty = FALSE`. |
| tna.R `.multilpa_modal_sequences` keys groups by `match()` on ids, not character levels | **Correct** | Numeric ids that print alike no longer merge. Row order and `group_class` stay aligned (the existing tna tests pass). |
| transition-inference.R roxygen for `seed` | **Correct** | Matches `.latents_local_seed()`: a NULL seed leaves the RNG stream alone, and a supplied seed restores `.Random.seed` on exit. |

## Findings

| id | severity | file:line | problem | reproducing input; actual vs expected | status |
|---|---|---|---|---|---|
| F1 | medium | R/lta-engine.R:595-645 (`.lta_setup_general`; no design-rank check) | `lta()` fits a rank-deficient transition or initial covariate design without any message and reports meaningless, unidentified logit coefficients. `multilpa()` refuses the same situation (R/covariates.R:841-844, "Covariate design is rank deficient"). `n_parameters`, AIC and BIC count the unidentified columns. Wald inference does refuse (`latents_singular_information`), but `get_results(fit, "transition_coefficients")` still prints the estimates. | (a) `d$cst <- 0.1; lta(d, c("y1","y2"), "id", 2, "time", transition_covariates = "cst")`, with sim n = 50, T = 4, seed 11 (scratch s8.R): no warning; coefficients (Intercept) -1.3877, cst 1.3884 for origin 1 (any split of the intercept fits equally well); `parameter_inference()` gives `latents_singular_information`. Same for `initial_covariates = "cst"`: (Intercept) -0.173, cst 2.116. (b) A realistic slip: `transitions = "occasion", transition_covariates = "tt"` with `tt = time` (s24.R). Stacked design rank 3 of 4, `n_parameters` = 17 where 15 are identified, no warning. Expected: a classed refusal or warning, as `multilpa()` gives. | reported, not fixed. Choosing refusal or warning is a judgement call. A correct rank check must also follow the order-2 column subsets (`.lta_transition_columns`). Proposed: before the starts, check `qr(stacked transition design over layout$within rows)$rank` (per order block) and `qr(initial design)$rank`, and raise `latents_bad_covariate`. Well-posed fits are unaffected. |
| F2 | low | R/lta-engine.R:380-470 (M-step keeps start values when there is no information) and `.lta_count_parameters` | With `transitions = "occasion"`, an occasion whose rows all have missing indicators (allowed under `missing = "fiml"`) leaves its transition logits with no information. The likelihood is flat in them, so EM returns the start values, which are then reported as estimates and counted as free parameters. | sim n = 50, T = 4, seed 11, plus 10 rows at time 5 with y1 = y2 = NA (s6/s7/s26/s27.R). LL unchanged by the extra occasion (-510.5949 both, correct), but `n_parameters` 15 -> 17. Fitted occasion_5 logits (-1.4876, -2.0337) equal the `max_iter = 0` start exactly. Setting them to (3, -4) leaves the LL at -510.594894208 to 12 digits. Under `measurement = "occasion"` the thin-profile warning flags the empty block ("0.0 Gaussian membership mass at occasion 5"); nothing comparable flags the transition logits. Expected: a warning or an `estimated = FALSE` flag, as the homogeneous engine gives for uninformative transition rows. | reported, not fixed (choice of diagnostic). |
| F3 | low | R/tna.R:269-285 (`.multilpa_modal_sequences`) | For a cross-sectional `multilpa(time = )` fit, an occasion with no observed indicator still gets a modal profile (from the prior), which then enters `get_tna()` sequences as an observed state. | `multilpa(..., missing = "fiml")` on sim n = 40, T = 4, seed 81, with three persons all-NA at time 2 (s21.R): `subject_profiles` has 0 NA, and the TNA sequences contain profile labels at those occasions. It is arguably a design decision of the cross-sectional classifier (outside my partition); it only shows up in tna.R. | reported, not fixed |
| V1 | info | R/lta-engine.R:51-75 | Transition covariates are read at the destination occasion. depmixS4 reads them at the origin. This is documented, not a defect. | Same data (s19/s20.R): depmixS4 `transition = ~w` best LL -806.1963 vs latents -806.4047. With depmixS4 given the lead covariate, best LL -806.4047, diff 2.8e-8. | verified-correct (convention) |
| V2 | info | R/lta-three-step.R:262-268 | The step-three chain stops with `latents_no_converge` when a move logit runs to the boundary. This is a loud failure, not a defect. | sim n = 60, T = 4, seed 21 (s11.R): optim code 1, largest score 2.7e-3. In s9.R the same data after a different deletion converged with intercept -3.1 (SE 3.2). | not a defect |

## Algorithms verified correct (computed checks)

All brute-force checks sum over every K^T state path (Rabiner 1989, eq. 17), with an
all-missing occasion contributing emission 1.

| Algorithm | Check | Max abs diff |
|---|---|---|
| Homogeneous scaled forward pass, unequal lengths, length-1 sequences | total and per-person LL vs brute force (s2, s22) | 5.7e-14 total, 3.6e-15 per person, 2.8e-14 with length-1 sequences |
| Homogeneous occupancy posteriors and expected move counts (xi) | vs brute-force path posteriors (s23) | 6.7e-16 / 7.1e-15 |
| General log-space recursion: occasion intercepts, transition covariate (destination), initial covariate | LL vs brute force (s3); posteriors and pair counts (s23) | 0 (LL), 1.1e-15, 1.4e-14 |
| Second-order chain (pair-state augmentation, index (i-1)K+j) | LL vs brute force (s4) | 1.1e-13 |
| Mover-stayer (stayer share 0.745) | mixture LL vs brute force (s5) | 1.1e-13 |
| Two group classes (homogeneous engine) | LL and group posteriors vs brute force (s4) | 1.1e-13 / 2.2e-16 |
| Sampling weights (pseudo-ML) | weighted LL vs sum w_i l_i with weights scaled to mean 1 (s4) | 3.4e-13 |
| Occasion-specific measurement | LL vs brute force (s4) | 1.1e-13 |
| Missing emissions: trailing all-NA occasion; grid skip vs explicit all-NA row | refit LL (s6) | 1.3e-11 / 1.1e-13; transitions 3.8e-8 |
| depmixS4 at the same parameters (homogeneous, 2 Gaussian responses, unequal ntimes) | `setpars` LL and `forwardbackward` gamma (s18) | LL 0, gamma 7.8e-16; depmixS4's own maximum equals latents' (-606.0401) |
| Homogeneous Wald SEs (observed and robust CR1) via delta method | numDeriv Hessian / score Jacobian of the brute-force LL (s14) | match to 6 decimals on all 14 natural parameters |
| General Wald SEs (observed, robust) for logits | numDeriv on brute-force LL (s15); gradient at optimum 5.0e-6 | match to 6 decimals |
| Weighted robust SEs | H = sum w l'', meat = sum (w s)(w s)', CR1 (s25) | match to 6 decimals |
| Parameter counts | hand count for 8 specifications (H=2, occasion, covariates, order 2, order 2 + occasion, stayer, occasion measurement, equal variances) vs `n_parameters` and packed length (s17) | exact |
| Bootstrap relabelling of general logits (K = 3, covariates) | implied transition and initial probabilities permute as P[o,o] (s16) | 1.1e-16 |
| Engine log-transition softmax | vs independent softmax (s16) | 1.1e-16; row sums 1 |
| `.lta_simulate` (used by the bootstrap LRT) | 4500 simulated sequences; logistic regressions of simulated moves on w and starts on z recover the coefficients (s12) | all within about 1-1.5 sampling SE (e.g. slope -0.113 vs -0.164; SE about 0.04) |
| Three-step (D, occasion posteriors, step-three LL, scores, estimates, SEs) | see the verdict table | 2.2e-16 to 1.2e-7 |
| Duplicate (id, time) | refused with `latents_bad_time` under both layouts (s22) | n/a |

## Tests added

`tests/testthat/test-review-1026-lta.R` (9 tests, about 5 s):

- homogeneous lta() matches the brute-force path sum on unequal sequences
- covariate transitions and initial states match the brute-force path sum
- second-order and mover-stayer likelihoods match the brute-force path sum
- the brute-force likelihood is invariant to relabelling the profiles
- relabelled general-fit logits imply the permuted probabilities
- an occasion with no observed indicator leaves the likelihood unchanged
- the step-three chain matches its brute-force likelihood and numerical scores (numDeriv guarded by `skip_if_not_installed`)
- lta() agrees with depmixS4 at the same parameters (`skip_if_not_installed("depmixS4")`)
- nested transition comparisons refuse different fitting bounds (error path, class `latents_bad_argument`)

Result: `testthat::test_file("tests/testthat/test-review-1026-lta.R")` gives
`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 37 ]`.

Mutation check: replacing `colMeans(posteriors[observed, , drop = FALSE])` with
`colMeans(posteriors)` in `.lta_occasion_posteriors` gave `[ FAIL 3 | WARN 0 | SKIP 0 | PASS 34 ]`.
The file was then restored from a byte copy.

## Related existing test files run

Run with `NOT_CRAN=true`, so no `skip_on_cran()` skips, under `load_all()` with the lead's
uncommitted fix and my new file. Command: scratch runtests.R, which calls
`testthat::test_file(reporter = "silent")` per file. Counts are `test_that` blocks.

| file | tests | failed | skipped | errors | seconds |
|---|---|---|---|---|---|
| test-lta-general.R | 6 | 0 | 0 | 0 | 7.2 |
| test-lta-compare.R | 5 | 0 | 0 | 0 | 15.6 |
| test-lta-bootstrap.R | 10 | 0 | 0 | 0 | 11.0 |
| test-lta-three-step.R | 4 | 0 | 0 | 0 | 2.6 |
| test-lta-acceleration.R | 6 | 0 | 0 | 0 | 0.3 |
| test-lta-covariance-boundary.R | 3 | 0 | 0 | 0 | 0.2 |
| test-lta-thin-profiles.R | 5 | 0 | 0 | 0 | 1.6 |
| test-transitions.R | 17 | 0 | 0 | 0 | 1.7 |
| test-transitions-numerics.R | 8 | 0 | 0 | 0 | 0.7 |
| test-transition-inference.R | 5 | 0 | 0 | 0 | 2.8 |
| test-three-step.R | 16 | 0 | 0 | 0 | 1.0 |
| test-tidy-three-step.R | 15 | 0 | 0 | 0 | 0.6 |
| test-sequences.R | 9 | 0 | 0 | 0 | 0.5 |
| test-tidy-sequences.R | 6 | 0 | 0 | 0 | 0.2 |
| test-tna.R | 8 | 0 | 0 | 0 | 0.9 |
| test-tna-fits.R | 4 | 0 | 0 | 0 | 0.5 |
| test-review-lta.R | 19 | 0 | 0 | 0 | 50.3 |
| test-review-1015-transitions.R | 12 | 0 | 0 | 0 | 4.7 |
| test-review-1026-constancy.R (lead) | 5 | 0 | 0 | 0 | 0.1 |
| test-review-1026-lta.R (new) | 9 | 0 | 0 | 0 | 4.7 |

No regression. (A first run without NOT_CRAN skipped 24 blocks and two whole files, test-transition-inference.R and test-review-lta.R, and also had 0 failures.)

## What I could not check

- LMest and seqHMM were installed but not compared: depmixS4 plus brute force already
  covers the same recursions, and the repository's equivalence/lta-extensions/ directory
  holds LMest/Mplus results that I did not re-run.
- Categorical, ordinal and count emissions inside the chain were not brute-forced. Only
  Gaussian diagonal emissions were. The recursion consumes an emission matrix, so those
  families differ only in `.latents_measurement_log_density` (outside my partition).
- Full-covariance (VVV/EEE) LTA SEs and the EM-only (non-diagonal) path were not checked
  against numDeriv.
- Bootstrap SEs (`method = "bootstrap"`) and `bootstrap_lrt()` p-values were not run end
  to end (cost). Their pieces were checked: relabelling, simulation and the nesting rule.
- F1 and F2 were not fixed: they need a design decision (refuse vs warn) and, for F1, a
  rank rule that follows the order-2 column subsets.
