# Round-2 review: finite mixtures of regressions (`mixreg`), 2026-10-26

Every number below was produced in this session (R 4.5.2, `load_all`,
`NOT_CRAN=true OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`). Scripts:
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-mixreg/`
(`NN_*.R`; `helpers.R` loads the package and defines the numDeriv comparison).
No file under `R/`, `tests/`, `equivalence/` or `man/` was edited.

## Scope

- `R/mixture-regression.R`, `-engine.R`, `-inference.R`, `-methods.R`,
  `-ordinal.R`, `-predict.R`, `-starts.R`, `-enumerate.R`.
- Read-only, because the mixture regression calls them: `R/kernel-structures.R`
  (E-step, multinomial membership M-step, membership scores, packing),
  `R/inference-service.R` (score information, sandwich, delta Jacobian),
  `R/extra-indicators.R` (NB dispersion derivatives).

Round-1 items MR-F1 to MR-F4 were re-checked and hold. MR-F6 (the
differencing step) is fixed. MR-F5 (origin and unit dependence of the
inversion) is fixed for units but only partly for origin; see N6.

## Findings

| id | sev | file:line | problem | reproduction | measured effect | proposed fix |
|---|---|---|---|---|---|---|
| N1 | medium | R/kernel-structures.R:247 (`.mixture_multinomial`), reached through `membership =`, `group_membership =` and the two-level class logits | The membership (concomitant) Newton step uses the raw ridge `1e-10 (1 + max diag)`. This is the defect fixed in round 1 for the regression and ordinal M-steps (MR-F3); the membership step was not changed. With a membership covariate far from its origin, the intercept–slope direction has curvature about n p(1-p) s²/(1+m̄²). The ridge is far larger than that and damps the step in that direction by a factor of about 10³. The M-step is then not a maximiser, and EM creeps. Because EM stops on the relative change in the log likelihood, it can report convergence short of the maximum. | `14_membership_origin.R`, `15_multinomial_mstep.R` | One M-step on fixed soft counts (n=400), against glm: shift 0 reaches the maximum (gap -1.7e-13). Shift 2010 stops 45.0 below the maximum of the objective (slope -0.0577 vs -1.200). Shift 1e4 stops 50.4 below. Full fits, membership `m = a + b t`: **a=2010, b=1: converged=TRUE after 225 iterations, LL 8.4e-5 low, membership slope 0.013 SE off, SE 0.13% off**. a=2010, b=6: 0.0011 SE off. a=1e4, b=1: unconverged after 1000 iterations, LL 47.3 low, slope 7.5 SE off, SE 43% off. Unit scaling alone (b up to 1e5) is fine (2.4e-6 SE). | Use the same Jacobi-standardised solve as `.mixture_update_regression` (solve `D^-1/2 I D^-1/2` when the ridge exceeds 1e-7 of a diagonal). Centring membership columns with \|mean\| >> spread before the solve and mapping the intercept back would also make it exact in origin. |
| N2 | medium | R/kernel-structures.R:275 (`separation = any(abs(coefficients) > 30)`), raised as `latents_separation` in R/mixture-regression-starts.R:203 | The membership separation flag tests raw coefficient magnitudes, which depend on units and origin. Round 1 replaced this test in the regression block with a test on fitted logits. An ordinary covariate with a mean of 30 (age), or a calendar year, raises `latents_separation` when nothing is separated. `enumerate_regressions(bootstrap > 0)` treats a replicate refit that warns `latents_separation` as invalid, and the p-value is `NA` when any replicate is invalid. So the bootstrap test is refused for such covariates. | `24_membership_separation_flag.R`, `25_blrt_membership_flag.R` | Identical fitted priors (range [0.008, 0.979]) for m = t, 0.01 t, 2010 + 6 t and 30 + t. The last three warn `latents_separation` (estimates -147.6 for the slope; 493.7 and 43.73 for the intercept). BLRT, 8 replicates, m = shift + 2t: shift 0 gives 7 valid and 1 flagged; shift 30 gives 6 and 2; **shift 2010 gives 0 valid, all 8 flagged, p = NA**, against an observed statistic of 421.5. | Flag membership separation on the fitted log-odds: \|w'gamma_k\| > 30 (or a prior within 1e-13 of 0 or 1) over supported units, as the regression block now does. |
| N3 | medium | R/mixture-regression-methods.R:648-654 (`.mixture_weighted_sample`) | For single-level fits with `id` (`class_level = "observation"`, `n_group_classes = 1`) and `weights`, `sampling_weights` holds one weight per group (R/mixture-regression.R:480-484). `.mixture_weighted_sample` therefore counts groups. Without weights the same model counts rows, and the classes and the likelihood are per row. Adding weights of exactly 1 changes `bic`, `sabic`, `icl`, `nobs()` and `logLik()`'s `nobs`. | `04_weights_units.R` | n=300 rows in 60 groups, K=2, k=7. Without weights: LL -496.747096, bic 1033.4207, nobs 300. With weights all 1: same LL, **bic 1022.1546, nobs 60** (difference k·log 5 = 11.27); sabic 1011.22 → 1000.14; icl 1257.07 → 1245.80. The per-parameter penalty differs, so the BIC-selected class count can differ between a weighted and an unweighted analysis of the same design. | For `nesting == "observation"`, count the rows with positive weight: `units = sum(row_weights > 0)` and `unit_kept = row_weights > 0`. The groups enter only the robust meat, as documented. |
| N4 | low | R/mixture-regression-engine.R:91 (`dnbinom` objective) with :417/:442, and R/mixture-regression-inference.R:215 | Near the Poisson limit, `stats::dnbinom(size = 1/alpha)` has rounding noise of about 1e-7 per 300 rows at alpha = 1e-8. The dispersion line search and the "move to the floor only if it does not lower the objective" test compare objective differences of that size. A class can therefore stop just above the floor: alpha = 1.0004e-8, log gap 4.2e-5 from the floor. The boundary test is `theta <= log floor + 1e-6`, so that class is not held. It then gets a standard error, and its `latents_boundary` warning does not name it. | `09_nb_near_floor.R`, `10_dnbinom_noise.R`; also found in 1 of 12 seeds of `08_nb_vs_poisson.R` | dnbinom minus the exact first-order expansion jitters between -2.1e-7 and +9.0e-8 across alpha = 1e-8 (1 + k·1e-5). Fit: class-1 alpha 1.000422e-8 (score at zero -29.5, so the KKT condition says floor). It reports dispersion SE 1.84e-5 instead of NA, and confint for log_dispersion.class_1 is [-3628, 3591]. Other SEs are unaffected (3.6e-9 relative). | Use the stable NB-minus-Poisson expansion (already in `.latents_negative_binomial_dispersion_derivatives`) for the objective when alpha·max(y, mu) < 1e-3. Alternatively, call a coordinate a boundary when its KKT score at zero is ≤ 0 and alpha < 1e-6. |
| N5 | low-medium | R/mixture-regression-ordinal.R:386-389 (separation check) and the inference boundary rule (R/mixture-regression-inference.R:213) | An ordinal class whose posterior assigns no row to a *middle* category has its MLE on the boundary: the gap between two thresholds is 0. The log gap runs to -∞, EM reports convergence, and no warning is raised. \|t - eta\| > 30 catches only an empty end category. The coefficient table shows two equal thresholds, each with a Wald interval. | `26_ordinal_dup.R` | K=2 ordinal, 4 categories, n=300 (weights 1-3) or 600 (duplicated rows). Class 2: log gap(2\|3) = -21.56. Posterior-weighted count of category 2 in class 2 is 1.2e-8. Moving the gap to e^-10 changes the LL by only -9.8e-5. No warning. The observed SE of that log gap is 33065; the robust SE is 2.07. Thresholds 1\|2 and 2\|3 are both -0.19 (SE 0.31). | Treat a log gap below log(1e-6) (relative to 1 + \|t\|) as a boundary, as for the NB dispersion. Hold that coordinate (NA SE), warn `latents_boundary` with the empty category named, and keep the remaining errors conditional. |
| N6 | low (residual of MR-F5) | R/mixture-regression-inference.R:236-254 (`.mixture_invert`) and :152 (information in the raw chart) | The correlation-form inversion fixed unit dependence. Origin dependence remains. The information is differenced in the raw chart, so the intercept–slope correlation is 1 - O((sd/shift)²), and the relative error of the differenced entries is amplified by about (shift/sd)². SEs lose accuracy quietly as the origin grows, then turn all-NA, including slopes that glm reports without trouble. The roxygen claims acceptance up to "a mean some 1e5 standard deviations from zero"; measured failure starts at 5e4. | `02_invariance_se.R`, `03_origin_threshold.R` | One class, n=300, slope SE against glm: shift 1e4: 1.6e-5 relative; 2e4: 1.5e-4; 3e4: 5.1e-4; ≥ 5e4: NA (gaussian and Poisson). K=2, shift 1e5: every SE NA in all four families. Shift 2010: SE error ≤ 3.8e-6. | Difference the score in the centred chart the M-step already uses, then map back exactly with V = J V_c J' (J the intercept map b0 = b0_c - c'beta). Correct the roxygen limit. |
| N7 | low | R/mixture-regression-methods.R:87-97 (`n_assigned`, `mean_posterior`) and :315-329 (`"classification"`) | ML-1 removed zero-weight units from every criterion and from entropy. The classes table's `n_assigned` and `mean_posterior`, and the average-posterior classification table, still count them. | `18_zero_weights.R` | Adding 6 zero-weight groups (30 rows): LL, BIC and SEs unchanged (≤ 1e-6), but `n_assigned` 99/101 → 125/105 (observation) and 20/20 → 26/20 (group), and classification means move by 1.7e-2. | Restrict `modal` and `posterior` to units with positive weight (`.mixture_weighted_sample()$unit_kept`, or `row_weights > 0` for the observation nesting). |
| N8 | low | R/mixture-regression.R:663 | The two-level model refuses when `n_groups <= n_classes`, but its classes sit on rows. 4 groups of 200 rows, 4 row classes and 2 group classes are refused with "4 groups cannot support 4 classes". (`n_groups <= n_group_classes` is checked separately at :670, and that check is right.) | `21_twolevel_refusal.R` | The refusal message is `latents_bad_data`; the same 4 row classes fit single-level (LL -1259.543). | Use `units <- if (nesting == "group") n_groups else n`. Keep the separate group-class check. |

Not a defect, for the record: in `08_nb_vs_poisson.R` seed 6, EM started at
the Poisson solution stays there (-640.9637). Multi-start finds -640.9622
with alpha₂ = 0.0084. The Poisson point satisfies KKT (score at zero ≤ 0),
so it is a genuine local maximum, and the start search handles it.

## Verified correct

| algorithm | check | max discrepancy |
|---|---|---|
| Observed information (today's scale-aware central differences of the analytic score; correlation-form inverse) | SEs vs inverse numDeriv Hessian of the LL, K=2, tol 1e-13 (`01`, `16`): gaussian; gaussian with common + membership; equal variance; Poisson with membership; NB varying and equal; binomial counts; ordinal; ordinal with common (± membership); NB with common + membership; gaussian group nesting with membership; Bernoulli group; Poisson group; gaussian two-level with membership + group membership; Poisson two-level K3 H2; NB two-level | 1.8e-9 to 2.6e-6 relative; 8.3e-6 for an unconverged NB two-level fit |
| Unit scores (Fisher identity, all nestings, weighted) | analytic vs numDeriv Jacobian of each unit's (weighted) log likelihood (`05`): gaussian with membership, weighted; group weighted; two-level, unweighted and weighted; NB; ordinal; observation with id | ≤ 6.6e-9 absolute (scale ~8); NB at the floor 6.5e-6 (boundary coordinate, held) |
| OPG | inverse crossprod of numerical unit scores (`05`) | 6.5e-11 relative |
| Cluster sandwich (CR1) | one-class weighted Poisson clustered on id vs `sandwich::vcovCL(type="HC0", cadjust=TRUE)`: 3.5e-7; row-weighted: 2.2e-7; K=2 group nesting vs numeric bread·meat·bread·G/(G-1): 1.6e-8 (`05`) | 3.5e-7 |
| One-class limits | gaussian (and with `common`) vs glm (SE × √((n-p)/n)); Poisson with offset; binomial counts; NB vs glm.nb (coef, alpha, LL) and numDeriv SE; ordinal (and with `common`) vs polr (`07`) | coef ≤ 3.5e-6 (polr's tolerance), LL ≤ 4.4e-9, SE ≤ 4.0e-6 |
| K=2 with concomitant model vs flexmix | Poisson + `FLXPmultinom(~m)`: ours LL higher by 3.2e-7, coefficients 1.6e-6. Grouped binomial `\| id` + concomitant: dLL 8.0e-10, coefficients 2.4e-8, df 6 = 6, bic_rows = flexmix BIC exactly (`06`) | 3.2e-7 LL |
| Weighted IRLS M-step with `common` | one M-step vs stacked weighted glm (class-specific x, shared z), Poisson, binomial, gaussian (`20`) | 3.9e-10 |
| Two-level likelihood and posteriors | brute-force enumeration of all 3^5 row-class assignments × 2 group classes, at a perturbed generic point (`12`) | LL 1.4e-14, rho 3.3e-16, tau 8.9e-16 |
| Unit and origin invariance of estimates and SEs (MR-F6 fix) | predictor ×1e4, ×1e-4, shift 2010, shift 2000 with ×1e3; gaussian, Poisson, NB, ordinal K=2 (`02`) | estimates ≤ 4e-10 SE; SE ≤ 3.8e-6 relative (shift 2010), ≤ 2e-12 (pure scale) |
| Outcome scale invariance | y × 1e-6 … 1e6 (and + 1e6 offset at ×1e3, ×1e6), gaussian K=2 with membership (`22`) | ≤ 3.9e-6 SE. The remaining cases are data-rounding limited: 1e6 + 1e-6 y₀ |
| NB dispersion KKT / GEM monotonicity (MR-F4) | 12 data sets: NB K=2 LL ≥ Poisson K=2 LL; EM started at the Poisson solution never decreases (`08`) | NB - Poisson ≥ -5.0e-7 (the floor alpha = 1e-8, not 0); EM gain ≥ 1.5e-10 |
| Delta method (scale-aware Jacobian) | sigma SE vs sigma·SE(log sigma); prior SE vs analytic softmax Jacobian (`19`) | 7.4e-11 |
| Ordinal thresholds (log gaps), threshold Jacobian | polr threshold SEs (`07`); numDeriv SE agreement on the log-gap coordinates (`01`, `16`) | 9.7e-7 |
| Criteria | AIC, BIC (units), SABIC, ICL (BIC + 2·entropy), relative entropy by hand (`13`); zero-weight units leave LL/BIC/ICL/entropy unchanged (`18`) | 0 (hand); ≤ 1e-5 (refit) |
| Sampling weights = duplication | integer weights vs duplicated rows: LL·W/n equal (Poisson 2e-12, ordinal 1.8e-11, NB 2.4e-9); coefficients equal (Poisson 1.6e-8). Ordinal and NB coefficients differ along flat boundary directions only (see N5) (`23`, `26`) | LL 2.4e-9 |
| Predictions | `response` = Σ prior·mean; `class_response` priors by hand for observation, group and two-level, ordinal and binomial counts; `posterior` on training rows reproduces fit tau; ordinal probabilities sum to 1, equal plogis differences, expected score = class_response (`11`) | 8.9e-16 |
| Enumeration table | grid rows = standalone fits (LL 5e-10, BIC 1e-6, from the seed stream); BLRT p = (1+b)/(1+B) = 0.1 for B=9 (`13`) | — |
| Separation detection (regression block) | Bernoulli group mixture with a class determined by sign(x): `latents_separation` raised, SEs 1.4e6 (`17`) | — |
| Parameter counts | `n_parameters` = length of packed theta in every fit of `01`/`16`; two-level count (H + ncol(w))(K-1) + ncol(v)(H-1); flexmix df equal | exact |

## Not checked (and why)

- K ≥ 2 NB and ordinal mixtures have no external reference installed. They
  were checked through numDeriv SEs, invariances, duplication and the
  one-class limits.
- The tau weighting of the regression separation flag (round-1 F7): my
  attempt to build a false positive (`17` b) did not reach \|eta\| > 30, so
  no new evidence.
- The definition of ICL for the two-level model (row posteriors only, no
  group-class entropy) and the BLRT design (p = NA when any replicate is
  invalid) are choices, not formula errors. Not tested further except as
  they interact with N2.
- Growth-mixture paths (`random =`, `cluster =`) belong to the growth area.
- Plots: no visual change was made, so no HTML render.
- The full test suite and R CMD check were not run (review only).

## Fixes (same day, at the coordinator's request)

All eight findings are fixed. No git was used. Regression tests are in
`tests/testthat/test-review-1026-mixreg-r2.R`: 10 tests, 49 passes.

### How the tests were checked against the old code

I copied the package to `.../scratchpad/round2-mixreg/oldpkg` with the
original `R/kernel-structures.R` and `R/mixture-regression*.R`, which are
backed up in `.../round2-mixreg/orig`. On that copy, every one of the 10 tests
fails:

| test | failures | errors | passes |
|---|---|---|---|
| N1 fit | 7 | | 1 |
| N1 M-step | 2 | | 1 |
| N2 flag | 3 | | 1 |
| N2 BLRT | 2 | | 0 |
| N3 | 6 | | 0 |
| N4 | 2 | 1 | 2 |
| N5 | 1 | 1 | 0 |
| N6 | 7 | | 0 |
| N7 | 4 | | 2 |
| N8 | | 1 (the refusal) | 0 |

### Per-file changes

| file | change |
|---|---|
| R/kernel-structures.R | N1: `.mixture_multinomial()` runs Newton in the orthonormal basis of the design's QR decomposition (coefficients R·beta, mapped back exactly) when column curvatures differ by more than the ridge allows. The gate is the MR-F3 rule, so ordinary designs keep their arithmetic. The Newton loop moved into `.mixture_multinomial_newton()`. N2: `separation` is read on fitted log probabilities (a category below exp(-30) on a unit with counts), not on \|coefficient\| > 30. |
| R/mixture-regression-engine.R | N4: new `.mixture_negative_binomial_log_density()`. Where alpha·max(y, mu) < 1e-3 it uses dpois plus the fourth-order NB/Poisson expansion; that expansion's error against Rmpfr is 4e-16, and dnbinom's noise was 2e-9 per row. The dispersion step now lands exactly on the floor, where the existing 1e-6 boundary rule holds it. |
| R/mixture-regression-inference.R | N6: new `.mixture_chart()`. It centres and scales x, z, w and v; centres are absorbed into class intercepts, ordinal first thresholds, membership intercepts or two-level group-class intercepts. It returns A and J = A^-1 in closed form (D^-1 - C D^-1). `.mixture_inference()` computes scores, information and inversion in that chart and maps back (V = J V_c J', information A' I_c A). N5: `.mixture_boundary_coordinates()` also holds log threshold gaps below `.ordinal_min_log_gap`. The `.mixture_invert` comment was updated. |
| R/mixture-regression-ordinal.R | N5: `.ordinal_min_log_gap <- -15`. The threshold table gives NA for every threshold whose Jacobian row touches a held gap. |
| R/mixture-regression-starts.R | N5: a `latents_boundary` warning naming the class and the empty category. |
| R/mixture-regression-methods.R | N3: `.mixture_weighted_sample()` counts rows of positive weight in the observation nesting, with or without `id`. N7: the classes table (`n_assigned`, `mean_posterior`) and the classification table use positive-weight units only. |
| R/mixture-regression.R | N8: classes are counted against rows except under `class_level = "group"`; group classes are still checked against groups. The roxygen Conditions section now names the ordinal gap boundary. **`man/` was not regenerated, so `devtools::document()` is needed.** |

### After the fixes (scratch scripts rerun)

| finding | before | after |
|---|---|---|
| N1 | membership at origin 2010: dLL -8.4e-5, 0.013 SE off; at 1e4: unconverged, 7.5 SE off | dLL ≤ 3.1e-12 and slope ≤ 8e-13 SE at 2010, 1e4 and 50000+1000t; 10 iterations each. One M-step reaches glm (gap 2.8e-14) at 2010 and 1e4. |
| N2 | false flags at 0.01t, 2010+6t, 30+t | no warning in any of the three; genuine separation is still flagged. BLRT (8 replicates): valid/flagged 7/1 at origins 0, 30 and 2010 alike, where 2010 gave 0/8 before. The test (seed 2, 3 replicates) gives p = 0.25 at both origins. |
| N3 | weights of 1 changed bic 1033.42 → 1022.15 | bic, sabic, icl and nobs are identical with and without weights of 1 |
| N4 | class-1 dispersion stuck at 1.0004e-8, with an SE | both dispersions at 1e-8, both held, confint NA; other SEs unchanged (2e-16) |
| N5 | gap e^-21.6, no warning, SE 33065 or 2.07 | `latents_boundary`; that log gap and thresholds 2\|3 and 3\|4 NA; threshold 1\|2 and the slope keep their SEs |
| N6 | SE error 1.6e-5 to 5e-4 at shifts 1e4 to 3e4; NA from 5e4 | slope SE equals glm (0.05945975 / 0.04404157) at every shift from 1e4 to 1e6. The K=2 invariance script gives rel dSE ≤ 3e-12 for gaussian, Poisson and NB, 1.5e-7 for ordinal at 1e5 (the fit-precision limit, F9). |
| N7 | n_assigned 20/20 → 26/20 with zero-weight groups added | unchanged (20/20 and 99/101) |
| N8 | refused | fits |

Rechecked after the fixes, all unchanged from before: numDeriv SEs (01, 16:
≤ 2.6e-6, or 8.3e-6 unconverged); scores, OPG and sandwich (05); flexmix
(06); glm, glm.nb and polr (07); predictions (11); two-level brute force
(12); delta method (19). The NB unit-score error at the floor fell from
6.5e-6 to 2.5e-10.

### Test runs

All test files matching `mixture|regression|ordinal|membership|separation|weight|review-1026`,
37 files, run after all changes:

- 0 failures, errors or warnings in 36 files.
- `test-review-1026-mixreg-r2.R`: 49 passes.
- `test-mixture-regression.R`: 78 passes.
- `test-growth-mixture.R`: 122 passes.
- `test-review-1026-mixreg.R`: 57 passes.
- The exception is `test-review-1026-inference-r2.R`: 5 failures and 3 errors. That is another agent's new person-centring file (multilpa, `R/centering.R`), outside these files. One error is "data must reproduce the original indicator data", which comes from accessors.

A first full run had transient errors in the growth, zero-weight and 1015
files while other agents were editing. Each passed on its own rerun and in
the full rerun above.

Harness `flexmix` via `run-suite.R`: 433 of 433 quantities agree (29,933
values), 0 outputs uncovered.

### Bootstrap rule (changed later, at the coordinator's decision)

Before this change, the bootstrap test treated a replicate refit that warned
`latents_separation` or did not converge as invalid, and the p-value was NA
when any replicate was invalid. Under the null, a spurious extra class with membership covariates
often drifts to a separating or slowly converging solution. On data where the
null holds, 4 of 12 K=2 refits did so after the fix, 3 with separation and 1 unconverged (`31_blrt_diag.R`).
So the bootstrap test is often NA with membership covariates. That is a policy
choice for the coordinator.

The rule now in R/mixture-regression-enumerate.R (`.mixture_blrt`):

- The p-value is (1 + b) / (1 + V) over the V valid replicates.
- The table gains the column `blrt_n_invalid`, and `print()` shows it.
  `blrt_flagged` is kept as the same count under its released name.
- `latents_failed_replicates` is warned, naming the count, whenever any
  replicate is invalid.
- The p-value is NA only when fewer than half the replicates are valid, or
  when the observed fits fail validation.
- The roxygen documents the reason: invalid replicates are expected under the
  null.

I did not use the shared `.latents_blrt_p_value()` in
R/simulation-service.R; the other families' `bootstrap_lrt()` still use it,
and I left them unchanged.

Test: in `test-review-1026-mixreg-r2.R`, with mocked refits, 1 of 3 invalid
gives p = 1/3 (V = 2), 2 of 3 invalid gives NA, and 0 of 3 gives 1/4 with no
warning.

`test-post-submission-review.R:18` used `dnbinom()` as the reference at
alpha = 1e-8, with tolerance 1e-12. That reference is 2e-9 off (Rmpfr). It
now uses the exact finite form sum log1p(j alpha) + y log mu - (y + 1/alpha)
log1p(alpha mu) - log y!. The test passes on the new code and fails on the old.

Final run of the 37 matching files: 1,606 passes, 0 errors, 0 warnings, and
4 failures, none from these changes.

- 2 in the `test-mixture-regression.R` snapshot, and 1 at
  `test-review-1015-regression.R:124`, where the expected
  `latents_separation` no longer fires: the threshold stops at 29.69 instead
  of 30.28.
  - Both come from another agent's concurrent edit to
    R/mixture-regression-engine.R: `.mixture_em_settings(stop_scale = ...)`,
    review "K2-6".
  - Shown by reverting only that edit in a scratch copy: the fit returns to
    LL -747.567194604, 113 iterations, threshold 30.281, and the warning
    fires.
  - With it: LL -747.567209039, 111 iterations, threshold 29.688, and in the
    snapshot LL -3144.0074 becomes -3144.0073 and 49 iterations become 53.
- 1 in `test-review-1026-growth-r2.R`, during the full run. It passed on
  rerun in both the working tree and a scratch copy, so it was transient
  while files were being edited.

### The same rule for `bootstrap_lrt()` (profile, additive and transition families)

The shared helper `.latents_blrt_p_value()` (R/simulation-service.R) now
applies the rule the regression bootstrap uses:

- the p-value is (1 + b) / (V + 1) over the V valid replicates;
- the Monte Carlo SE is sqrt(p(1-p)/(V+1));
- `n_invalid` is returned;
- `latents_failed_replicates` is warned, naming the count, whenever any
  replicate is invalid;
- the p-value is NA only when fewer than half the replicates are valid.

Other changes:

- **Callers.** R/enumeration.R, R/additive-compare.R and R/lta-compare.R store
  `n_invalid` in the result.
- **Output.** `get_results(x, "test")` and `as.data.frame()` gain an
  `n_invalid` column, falling back to iter - n_valid for older objects.
  `print()` reports "(k invalid)" and says why the invalid replicates are left
  out, or why the p-value is withheld.
- **Roxygen.** The `bootstrap_lrt()` description and its Conditions section
  give the reason.
- **Two existing tests** encoded the old rule and were updated:
  `test-review-1026-selection.R` (1 of 9 invalid now gives (1+4)/9; 5 of 9
  invalid gives NA) and `test-alignment-service.R:256` (1 of 3 invalid now
  gives 2/3, with n_invalid 1).

New tests in `test-review-1026-mixreg-r2.R` (now 15 tests, 84 passes):

- the shared helper on its own: 1 of 3 invalid gives a p-value, 2 of 3 gives NA;
- one test per family (profiles, additive, lta), each running the real
  verb with the replicate loop mocked to 1 of 3 invalid. Each checks n_valid
  2, n_invalid 1, the p-value over the two valid statistics, the warning, and
  "1 invalid" in `print()`.

On the old code, all 15 tests in the file fail.

Runs:

- `Rscript equivalence/run.R mplus-tech14`: 13 passes, 0 failures.
- The 37 test files matching `bootstrap|blrt|tech14|compare` plus every file
  that references `bootstrap_lrt`, `blrt` or `failed_replicates`: everything
  passes except these.
  - `test-mixture-regression.R`: the 2 snapshot failures already reported
    (the K2-6 stop rule).
  - `test-audit-reporting-repairs.R`: 1 to 2 failures, at line 17
    (variance floor) and line 45/48 (trust plot). The trust-plot failure
    also occurs in a fresh copy of the tree with my four bootstrap files
    reverted. The variance-floor failure did not reproduce in that copy, so
    it was a transient state of the tree; neither touches the bootstrap.
  - `test-informative-rows.R:42` (ICL): fails in the same reverted copy.
  - `test-search.R:175`: the parallel workers load the installed package,
    which lacks `.latents_collect_conditions`, another agent's new function.
