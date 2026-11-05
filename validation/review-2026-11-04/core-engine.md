# Accuracy review 2026-11-04: `core-engine`

Probes: `validation/review-2026-11-04/probes/core-engine-*.R`. All R ran through the two-slot wrapper, one process, `cores = 1`. Nothing in R/, man/, tests/ or DESCRIPTION was touched. Mutants were made on throwaway copies in the scratchpad.

## 1. What was reviewed

| file | depth |
|---|---|
| R/multilpa-starts.R | read in full; finalise (criteria, counts, noise, weights) checked by hand-computed fits |
| R/multilpa.R | E-step, M-step, EM wiring, setup, indicators, counting, single-level path read; roxygen of `multilpa` read as rendered |
| R/kernel-em.R, R/acceleration.R | read in full; SQUAREM and stopping exercised on 17 model set-ups |
| R/kernel-blocks.R, R/gaussian-moments.R | read in full; likelihood and FIML by hand |
| R/kernel-structures.R | nested and logit structure and `.multilpa_weigh` read and exercised. The mixture and forward-backward parts were not re-read: the 2026-11-01 audit checked them by brute force |
| R/covariance-structure.R | resolve, counts and structure list read; M-steps checked only through whole-fit stationarity (the closed forms are mclust-verified, not re-derived) |
| R/categorical.R, R/extra-indicators.R | density, M-step, counts, Poisson limit, prediction encoding read; inference coordinates and scores read only |
| R/local-dependence.R, R/direct-effects.R, R/noise.R | read; likelihoods and counts by hand |
| R/prior.R | read; not re-derived (mclust equivalence exists, including noise + prior) |
| R/centering.R, R/scaling.R (header only), R/weights.R, R/predict.R | read in full, exercised |
| R/fixed-measurement.R | read in full; held counts by hand |
| R/input-checks.R, R/search.R, R/lpa.R, R/multilca.R | read |
| R/multilpa-pool.R, R/start-pool.R | first part read only; not run, because they start workers |

Not reviewed: the additive and cross-level families, covariate models, transition engines, all inference and standard errors, and the classification statistics, which belong to other areas. No stress test of `.multilpa_structure_covariances` at the floor (the audit did this with `constrOptim`).

Note on the brief: it says kmeans is the default search. The code defaults to emEM for EII, EEI and EEE (`.latents_default_search`, R/search.R; commit c7c1cb9b, 2026-10-23; NEWS and the `search` help say so). Code and documentation agree, so this is not a finding, but the owner's note and memory entry are stale.

## 2. Findings

Counts: 0 CRITICAL, 2 MAJOR, 4 MINOR, plus EDGE items. All are CONFIRMED unless marked otherwise.

### F1. MAJOR, CONFIRMED. Person centring turns a singly observed value into a fabricated exact zero, and BIC then selects a spurious profile

- **Where:** R/centering.R:315-322 (`.multilpa_check_centered()` refuses only when *every* group is flat on an indicator); the help claim is at R/multilpa.R:1097, man/multilpa.Rd *Scale*: "Person centring needs each group observed more than once on each indicator".
- **What is wrong:** a group observed once on an indicator has deviation 0 by construction. That 0 carries no information, yet it enters the likelihood as an observed value. The code accepts this silently unless all groups are that way. `latents_atomic_groups` fires only under `scaling = "person"`, which refuses singletons anyway. Several such rows are identical all-zero vectors, so a profile can collapse onto them.
- **Evidence** (80 groups, 12 of them singletons, 482 rows, two true profiles, VVI, `centering = "person"`):
  - With 3 profiles, the fit has `profile_1` with means 0.000, 0.000, 0.000, variances 0.00000, share 0.025, and all 12 zero rows are modal in it. The only warning is the generic `latents_boundary`, which does not name the cause.
  - `boundary = "exclude"` cannot help: all five starts reach the spike (log likelihood -1766.447 either way).
  - BIC with the singletons: k = 2 gives 3992.752 and k = 3 gives 3620.534, so it picks the spurious 3, by 372. BIC with the 12 groups dropped: 3906.904 and 3914.515, so it picks 2.
  - Without a spike (k = 2), means still move by up to 0.028.
  - Through FIML, an indicator observed once in 28 groups gives 40 fabricated zeros on it (with the 12 singletons) and no warning at all.
- **Probe:** probes/core-engine-person-singletons.R
- **Affected:** `centering = "person"` on unbalanced data, without person scaling. Longitudinal and ESM data routinely have some persons with one observation.
- **Smallest fix:** a classed warning naming the indicators and the number of groups changes no number. Treating those deviations as missing under `missing = "fiml"` and refusing them under `"error"` is the principled repair, and it changes numbers for exactly those fits. Owner's decision.

### F2. MAJOR, CONFIRMED. Under weights, grand centring and grand scaling use unweighted means and SDs, so integer weights do not act as repeated units

- **Where:** R/centering.R:192 (`colMeans`), R/centering.R:116-120 (grand scale); the claim is at R/multilpa.R:1105, man/multilpa.Rd *Weights*: "so integer weights act as repeated units".
- **What is wrong:** `carrying` removes zero-weight units but ignores the size of the positive weights. A weighted fit's reported "grand-centred" means are centred on the unweighted sample. Under a structure not invariant to rescaling (EII, VII, EVE, VVE, EEV, VEV) the solution itself differs, not only its units.
- **Evidence** (n = 400, single level, one class oversampled, weights 3 and 1, against the duplicated rows, tol 1e-10):

| model | centring / scaling | reported means differ | back on raw scale | variances differ | shares differ |
|---|---|---|---|---|---|
| EII | none / none | 0 | 0 | 0 | 0 |
| EEV | none / none | 0 | 0 | 0 | 0 |
| VVI | grand / none | 0.727 | 0.00001 | 0.00001 | 0 |
| VVI | grand / grand | 0.520 | 0.00001 | 0.359 | 0 |
| EII | grand / grand | 0.557 | 0.185 | 0.101 | 0.016 |
| VII | grand / grand | 0.552 | 0.125 | 0.155 | 0.019 |
| EEV | grand / grand | 0.522 | 0.033 | 0.143 | 0.005 |
| VEV | grand / grand | 0.515 | 0.046 | 0.158 | 0.0003 |

  The population-weighted mean of the grand-centred means on the first indicator is -0.727, where a reader expects 0. Person centring and scaling are unaffected, because weights are constant within a unit.
- **Probes:** probes/core-engine-weights-grand.R, and the last table of probes/core-engine-weights-dup.R
- **Affected:** `weights` with `centering = "grand"` or `scaling = "grand"`.
- **Smallest fix:** either compute weighted offsets and scales at the two lines (this changes reported numbers for those fits, and for non-invariant structures the solution), or keep the behaviour and qualify the *Weights* and *Scale* text (no number changes). Owner's decision.

### F3. MINOR, CONFIRMED. `boundary` covers Gaussian spreads only, so it reads FALSE on fits the package itself treats as boundary fits

- **Where:** R/gaussian-moments.R:175 (`.multilpa_covariance_boundary()`), R/multilpa-starts.R:172, surfaced as "boundary status" in `get_results(x, "model")` and on the grids. Condition catalogue: R/conditions.R:337.
- **Evidence** (probes/core-engine-combos.R, probes/core-engine-floors.R):

| fit | estimate | warning | `fit$boundary`, model table | Wald |
|---|---|---|---|---|
| negative binomial | both dispersions 1e-08 | `latents_boundary` | FALSE | refused, `latents_boundary_fit` |
| LCA | a response probability at 1e-10 | none | FALSE | refused, `latents_boundary_fit` |
| Poisson, a class of zeros | mean 1e-10 | none | FALSE | refused ("a count mean is estimated at zero") |

  So an all-categorical fit can never report `boundary = TRUE`. `count_model` calls a zero dispersion "a boundary solution". The `latents_boundary` entry names only variances and eigenvalues, although the condition is also raised for dispersions, and `latents_boundary_fit` does not list count means.
- **Affected:** anyone reading the flag or the grids' boundary column on LCA, ordinal or count models.
- **Smallest fix:** say on `get_results` that the flag is about Gaussian spreads, and complete the two condition entries (no number changes). Extending the flag changes no estimate, but changes the column and any rule that reads it. I did not trace those readers; they belong to the enumeration area.

### F4. MINOR, CONFIRMED. The `drop_empty` help misstates the criteria

- **Where:** R/multilpa.R:995, man/multilpa.Rd:230 ("The information criteria of a profile or class fit count the informative rows either way") against R/multilpa-starts.R:191.
- **Evidence:** single-level fit, 200 rows, 197 informative. `bic` is 1298.392 under `TRUE` and 1298.528 under `FALSE`. Only `bic_individual` is unchanged (1298.392 both). Probe: probes/core-engine-combos.R.
- **Smallest fix:** say that `bic_individual` counts informative rows either way, while `bic` and `bic_groups` count all groups of positive weight under `FALSE`. No number changes. Open Issue 2 in HANDOFF is the transition-model relative of this.

### F5. MINOR, CONFIRMED. A separated ordinal indicator reports an arbitrary location with no flag

- **Where:** R/extra-indicators.R:308 (`.latents_ordinal_maximize()`); there is no check at finalise.
- **What is wrong:** when one class only ever uses the top category and no row in the other class's range offsets that, the location has no finite maximum. The fit returns a number set by the stopping rule, with `converged = TRUE`, `boundary = FALSE` and no warning.
- **Evidence:** the location is -21.62 at tol 1e-6 and 1e-8, -21.84 at 1e-10 and -23.13 at 1e-12, while the log likelihood stays at -832.3386. `parameter_inference()` refuses with `latents_singular_information`, so no false standard error is produced, and `get_results(fit, "ordinal")` shows the location with `NA` errors. The category probabilities are fine (1.0000 on the top category). This is the ordinal counterpart of the response-probability floor, which Wald reports as a boundary. Probe: probes/core-engine-floors2.R.
- **Note, tried to refute:** when the separating class overlaps on the continuous indicator, the location is a genuine finite maximum (4.05, SE 0.75 in probes/core-engine-floors.R). So only true separation is affected.
- **Smallest fix:** a `latents_boundary`-class warning when a fitted category probability underflows. No estimate changes.

### F6. MINOR, CONFIRMED by probe and by reading. The scoring of ordinal categories is undocumented, and gaps depend on the column type

- **Where:** R/extra-indicators.R:73 and :196-197; `ordinal` help, R/multilpa.R roxygen.
- **Evidence:** values {1, 2, 5} as integers get scores 0, 1, 2 (rank, not value; n_parameters 8, log likelihood -919.3665). The same column as a factor with levels 1:5 gets scores 0 to 4 and a warning (n_parameters 10, log likelihood -919.2848). Probe: probes/core-engine-floors.R.
- **By reading only (line 73):** an unordered `factor()` is accepted and its alphabetical level order becomes the score order, although the help says "an ordered factor or integers". "high" < "low" < "med" would be scored 0, 1, 2 without a word.
- **Smallest fix:** one help sentence on rank scoring, and optionally a warning for an unordered factor. No number changes.

### EDGE

- `predict(fit)` fails with `latents_bad_data` on a valid fit whose zero-weight group has a degenerate person scale (fit factor 1, R/predict.R:199). The fit itself tolerates this, and the help promises such rows keep their predictions. It is a refusal, not a wrong number. Probe: probes/core-engine-predict2.R.
- `missing = "error"` refuses an `NA` in a continuous value on a zero-weight row (R/multilpa.R:1915) but accepts one in a categorical value. Same probe family: core-engine-weights-grand.R.
- A new one-row unit under person centring is classified at deviation 0 whatever its values: `y1 = 99, y2 = -99` gives profile 2, posterior 0.657 (core-engine-predict2.R). This is the prediction face of F1.
- The default `min_variance` ignores weights (2.268e-06 against 2.332e-06 duplicated). It is documented, and it matters only at the bound.
- The comment at R/multilpa.R:704 says prior fits stop on the likelihood's relative change, "mclust's rule". The code uses the unit-free `tol` per observed value. Comment only.

## 3. Verified sound

Each check below was run here, against plain base R that uses none of the package's helpers unless stated.

**Likelihood and FIML**
- All 14 structures under FIML (2 profiles, 3 indicators, 12% missing cells): hand log likelihood differs from the fit by 0 in every case. BFGS from the fit, in each structure's own coordinates (volume, determinant-one shape, rotation), gains at most 1.8e-11, so the constrained M-steps are stationary under missing data, which mclust cannot check. The check can fail: fits stopped after 3 EM steps show gains of 0.009 (VVI), 0.054 (EVE) and 1.08 (VEV).
- Two-level, 2 x 2 classes, continuous + categorical + ordinal + count, 10% missing per indicator: agreement to all 12 printed digits (Poisson -1604.73373363; NB -1601.92545207; shared dispersion -1601.99536365). Finite-difference scores of the ordinal intercepts and location, count means and dispersions all round to 0 at five decimals.
- Continuous pair + FIML: hand difference 1e-13, optim gain 3e-13. Direct effects + FIML (VVI, EEI): hand difference at most 1e-13, optim gain at most 2e-12. Half-observed categorical pair: hand difference 0.
- Row order with interleaved groups: log likelihood 6e-14, row posteriors 6e-16, and group posteriors by hand 1e-15. `vars` order with a direct effect: log likelihood difference 0, means difference 0.
- Evaluate-only round trip through `starting_values()` (mixed two-level, three transforms, EVE): log likelihood difference 0, posteriors at most 2e-15.
- By derivation: ordinal density and the NB score (d log f / d log alpha) match the code, the Poisson log-ratio offset is exact, and the bounded multinomial solution is the KKT one.

**Parameter counts** (hand, and for structures also the numerical rank of the parameterization's Jacobian)
- All 14 structures: `n_parameters` = hand = rank + means + mixing.
- Mixed two-level: 21 (Poisson), 23 (NB), 22 (shared dispersion).
- Categorical pair 15; continuous pair 15 (VVI) and 11 (EEI); direct effect 14; noise 15, as mclust counts it.
- Held blocks: `fixed = "means"` 7 with measurement 11; `"measurement"` 1, with the values held exactly; VEV with held means 6.

**Weights:** 10 settings against duplicated groups (EEI, EVI, VEV and EVE under FIML; VVV with 2 group classes; mixed NB; Poisson; direct effects; both pair kinds): rescaled log likelihood within 1.6e-9, means within 3e-6, with equal counts. Person centring and scaling match too (means within 6e-5 at tol 1e-10; both runs used their own starts). The exception is F2.

**Scale:** rescaling units leaves EEI, VVI, EEE and EVI within 7e-11 of the Jacobian shift, and VEV is off by 0.047, as documented. The help's list of scale-invariant structures matches the derivation from C Sigma C.

**Noise:** hand log likelihood -1023.41726384 equals the fit; `hypervolume` equals `mclust::hypvol()` (178.4578); shares sum to 1.

**SQUAREM and stopping** (17 set-ups: 12 structures, mixed NB, weights, direct effect, pair, discrete only)
- The likelihood history never falls (smallest step +1.3e-10), and its last entry equals the reported value.
- It reaches plain EM's value within 5e-9, in 33-323 steps against 157-1395.
- `max_iter` of 1, 2, 3, 4, 5, 7 gives exactly that many iterations, with `converged = FALSE`.
- A landed point below the two plain steps is rejected, and a failed extrapolation falls back, so an accelerated cycle cannot report a lower likelihood.

**What the flags mean, as implemented**
- `converged`: the last kept step or cycle changed the log likelihood by at most `tol` times the observed values. There is no gradient test. At the default tol on one fit, the next EM step gained 1.6e-6 (threshold 6.9e-6) and the run stood 3.1e-5 below the maximum. It can be TRUE at a bound, and is FALSE at `max_iter = 0`.
- `boundary`: only the smallest Gaussian variance or eigenvalue at `min_variance * (1 + 1e-8)`. See F3.

**predict()**: 20 settings (structures, FIML, mixed, weights, direct effects, both pair kinds, seven centring and scaling combinations, noise, prior, LCA) reproduce the fit's posteriors to 8e-16 with 0 class differences. Single-level `density` sums equal the log likelihood. Zero-weight fits match to 4e-16, and permuted `newdata` to 2e-16.

**Not established:** selecting starts by likelihood while a prior fit climbs the posterior could in principle pick a lower posterior. The help states the rule. On two data sets, 55 single-start runs reached one solution each, so the probe could not separate them (probes/core-engine-floors.R, core-engine-floors2.R). A data set with competing maxima would settle it.

A probe note: the last block of probes/core-engine-squarem.R stopped with `latents_bad_start` because my script passed a one-group-class start to a two-class fit. It tested nothing about the package. Held values under SQUAREM were verified instead in probes/core-engine-mixed-ll.R.

## 4. Test gaps

I ran only 2 of 16 planned mutants, because the shared slots were saturated and I stopped the rest to leave room for the other reviewers. Everything marked "by reading" below is PLAUSIBLE until someone breaks the formula and runs the file.

- **Noise parameter count** (run): changing `+ 2L` to `+ 1L` passes tests/testthat/test-noise.R (17 expectations, 0 failed). The count is asserted only in tests/equivalence/test-noise.R:69,96 and equivalence/test-mclust-noise.R, which are not in the built package; grep found no other shipped test.
- **Pair count** (run): removing the per-profile multiplier passes tests/testthat/test-local-dependence.R (50 expectations, 0 failed). This is probably not a real gap: by reading, test-features-pair-covariance.R:63 should catch it. I did not run that file against the mutant.
- **Hypervolume, by reading:** shipped tests assert only translation invariance (test-noise.R:65), so dropping the principal-component box would pass there. Equality with `hypvol()` is in tests/equivalence only.
- **Weights with grand transforms:** test-weights.R tests duplication for structures, categorical, FIML and two levels, but never with `centering` or `scaling`. That is why F2 survived.
- **Partial singletons:** test-centering.R:131 tests only the case where everything is flat. By grep, nothing tests a mix. That is why F1 survived.
- **Discrete and count floors:** test-extra-indicators.R:253-255 asserts the warning and the Wald refusal for the NB floor, but not the `boundary` flag.
- **Separated ordinal indicators:** by grep, no test asserts what a fit reports there.
