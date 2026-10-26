# Round-2 review: growth mixtures (slug `growth`)

## Scope

Files: `R/growth-mixture.R`, `R/growth-mixture-engine.R`,
`R/growth-mixture-inference.R`, `R/growth-mixture-methods.R`,
`R/growth-mixture-parameter-inference.R`, `R/growth-mixture-predict.R`,
`R/growth-mixture-correlated.R`, `R/growth-mixture-simulate.R`,
`R/growth-mixture-plot.R` (formulas only) and `R/trajectory-views.R`. The
growth branch of `mixture_regression()` (`R/mixture-regression.R:373-410`)
was read because it calls the new standardization. So were the shared
helpers the growth inference uses (`.inference_score_information`,
`.mixture_invert`, `.mixture_delta`, `.latents_quasi_newton`).

All runs used `NOT_CRAN=true OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1` and
`pkgload::load_all()`. The scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-growth/`
(`helper.R` holds the data generator; the checks are `c01` to `c20`). Round
1's G1, G2, G3 and ML-1 are not reported again here. The outcome-scale part
of G2 was re-verified and holds (see "Verified correct").

The growth test files on the current tree (`test-review-1026-growth.R`,
`test-growth-mixture.R`, `-correlated`, `-predict`, `-r3step`,
`-parameter-inference`, `test-multilevel-growth.R`) gave 0 failures, 0
errors and 291 passes.

## Findings

The common thread: today's standardization makes the fit free of the
**outcome's** units, but the **predictors'** units and origin still enter
three places. These are the variance floor, the Cholesky jitter and the
degeneracy rule, and the optimizer's coordinates. The time variable of a
random slope is the predictor where this matters most.

| Id | Sev. | file:line | What the code does / what it should do | Reproduction | Measured effect | Proposed fix |
|---|---|---|---|---|---|---|
| GR2-1 | **high** | `R/growth-mixture-engine.R:378-381` (M-step floor on `diag(G)`); `R/growth-mixture-inference.R:283` (L-BFGS-B bound `log L_jj >= 0.5 log(min_variance)`); `R/growth-mixture-correlated.R:217-220`; `R/growth-mixture.R:487` (`at_floor`); contract in `R/mixture-regression.R:135` | `min_variance` is documented as a floor for a **residual** variance, relative to var(y). It is also applied to every random-effect variance and to every Cholesky diagonal. Those are in units of y² per unit² of the random-design column. For a random slope on time measured in days, the true slope variance (about 1e-6 var(y)) sits below the floor (3.64e-6), so the slope variance is pinned there. A floor must be invariant to the units of the random design (Verbeke & Molenberghs 2000: the LMM likelihood is invariant to a reparametrization b = A b*, with G = A G* A'). | `c10_time_units.R`, `c11_time_units_cause.R` (two classes, `random = ~ 1 + time`, the same data with time x 365) | Days vs years, default `min_variance`: LL **11.31 lower**; slope variances x2.03 and x1.85; class-1 intercept-slope correlation 0.17 instead of 0.67; slope estimate 0.88 SE off; coefficient SEs x1.31 to x1.95. **`converged = TRUE` and no boundary warning.** With `min_variance = 1e-12` the days fit equals the years fit (LL difference < 1e-9), which isolates the cause. In hours (x8760): LL 462 lower, a boundary warning and NA SEs. | Floor each random-effect variance in the metric of its column: `min_variance_j = min_variance / mean(z_j^2)`, or `/ var(z_j)` for a non-constant column. Better: standardize the random design internally (Z* = Z A^-1, map back G = A G* A' exactly, as G2 does for y) and apply all floors in the starred coordinates. Document that the floor also governs random-effect variances. |
| GR2-2 | medium (wrong refusal / spurious boundary) | `R/growth-mixture-engine.R:220-222` (`.growth_factor`, jitter `1e-12 * max(diag(G))`), with the L_jj bound of GR2-1 | With a time origin far from 0 (calendar year), G in those coordinates is ill-conditioned: diag 1.57e6 and 0.39, smallest eigenvalue 4.9e-7. The jitter, 1e-12 x 1.57e6 = 1.57e-6, is 3x the smallest eigenvalue. The density is then evaluated at a different G, and EM is no longer monotone. The likelihood is invariant to the origin of time when the random design holds an intercept and time, so the fit should be too. | `c14_calendar_year.R`: `y ~ year, random = ~ 1 + year`, year = 2000-2004. `c13_time_origin.R` for origin +10/+100. `c15_jitter_cause.R` | 1 and 2 classes: **error `latents_em_decrease` "EM likelihood decreased by 0.384 ... This is a defect"**. lme4 fits the same model (LL -1233.46, with convergence warnings). Origin +100: both classes flagged degenerate (correlation -1.00 printed as "boundary: no interval"), `latents_unconverged`, slopes 0.077 SE off, LL 3.4e-3 lower. Origin +10: fine (4.7e-5 SE). With a per-column jitter patched in, the error disappears, but the K=1 fit is then 87 LL units below the time-origin-0 fit while reported converged: the L_jj floor (GR2-1) binds the conditional variance. | Standardize/centre the random design (see GR2-1), or at least use a jitter relative to each diagonal entry (`1e-12 * diag(G)`) and an L_jj bound expressed in standardized coordinates. Detect degeneracy in the centred coordinates (GR2-3). |
| GR2-3 | medium | `R/growth-mixture.R:491-492` (`min(values) <= 1e-6 * max(values)`) | The eigenvalue-ratio rule compares eigenvalues of variances in different units (intercept in y², slope in y²/time²). It flags a well-identified covariance as degenerate when time is in small units. The flag turns off every random-effect SE and interval and raises `latents_random_boundary`. The 0.995 correlation rule is unit-free; this rule is not. | `c11_time_units_cause.R` (last block): the years fit's G re-expressed in other time units (`D G D`, D = diag(1, 1/c)), with a negligible floor | Correlation 0.67 and -0.37, plainly interior: flagged "class 1, class 2" degenerate at c = 1000 and c = 8760 (hours); not at c <= 365. | Apply the eigenvalue rule to `cov2cor(G)` (or to G in standardized-Z coordinates), like the MR-F5 fix of `.mixture_invert`. |
| GR2-4 | medium | `R/growth-mixture-inference.R:278-302` (L-BFGS-B in raw coefficients and Cholesky elements); the EM relative stop | The finish is judged converged by L-BFGS-B's relative reduction of the objective. Its coordinates are in the predictors' raw units and origin, which today's standardization of y does not touch. With a poorly scaled time variable it stops short while reporting `converged = TRUE`. The correlated path instead finishes with Newton steps and a criterion in SE units (displacement < 1e-4), and reaches about 1e-9 SE. | `c12_time_conditioning.R` (`min_variance = 1e-14`, which removes GR2-1); `c18_year_fixed_only.R` (`random = ~ 1`, calendar year in the fixed part only) | Time x 8760: LL 0.0146 lower, coefficients 0.12 SE off, remaining Newton step 0.157 SE, `converged = TRUE`. Time x 1e5: LL 0.106 lower, 0.076 SE (Newton step 0.35 SE), `converged = TRUE`. Time x 365: 1.1e-3 SE off, honestly `converged = FALSE`. Calendar year in the fixed part: slopes 0.019 SE off, LL 3.5e-4 lower, SE 0.15-0.18% off, `converged = TRUE`. | Extend `.growth_standardize()` to centre and scale the design columns (x, z, random design) with an exact back-map. Alternatively, give L-BFGS-B `parscale` from an SE scale and end with the correlated path's Newton steps and SE-unit criterion. Report `converged` from that criterion. |
| GR2-5 | low | `R/growth-mixture-methods.R:189-209` (`.growth_class_table`) | `assigned_persons` (`tabulate(modal)`) and `average_posterior` count persons of sampling weight 0. NEWS (ML-1) states that such units are in no "classification table". `expected_persons`, `share` and the criteria correctly leave them out. | `c09_zero_weights.R`: 8 of 120 persons at weight 0, vs the fit without them | All criteria, shares and SEs agree (to 3.5e-7). Assigned persons 69/51 vs 64/48; average posterior 0.9508/0.9208 vs 0.9493/0.9189. | Restrict `modal`/`posterior` to `.growth_weighted_persons(x)$kept` in those two columns, or document the convention. |
| GR2-6 | low | `R/growth-mixture-methods.R:341` (two-level `icl`) | Uses `bic + 2 (sum of marginal person entropies + H(rho))`. That is neither the entropy of the joint posterior of (group class, person classes), H(rho) + sum_i sum_h rho_jh H(p_i\|h), which is the ICL-BIC definition of Biernacki, Celeux & Govaert (2000) applied to the complete data, nor the two-level regression family's convention (`R/mixture-regression-methods.R:244`: person entropy only, no group term). | `c19_icl_twolevel.R` | Entropy used 26.87, joint entropy 25.16. ICL 2461.17 (reported) vs 2457.74 (joint) vs 2447.11 (regression convention). | Use the joint-posterior entropy in both two-level families, or document one convention for both. |
| GR2-7 | low (read, not reproduced) | `R/growth-mixture-parameter-inference.R:93-94` | `residual_held` is `block == "residual_sd" & any(startsWith(on_bound, "log_sigma."))`. If one class's log sigma is held, every class's residual-SD row loses its SE and interval. It should match per class (`sprintf("log_sigma.%s", class)`, or `"log_sigma.shared"`). | read; it needs a correlated fit with one class's residual variance on its floor | — | Match on the class's own coordinate name. |
| GR2-8 | low (read) | `R/growth-mixture-inference.R:281-283` with `R/kernel-em.R:658` (`theta <- pmax(theta, lower)`) | EM floors variances (`diag(G)`), but the finish bounds the Cholesky diagonals L_jj, which are conditional SDs. An EM point with a strongly correlated G (L_22^2 = G_22(1 - r^2) below the floor) is moved onto the bound before L-BFGS-B starts, and "improved" is then judged against the moved point, not EM's. The two stages maximize over different parameter sets. This contributes to GR2-2 (an 87-unit deficit under the patched jitter). | read; see `c15_jitter_cause.R` | — | Use one constraint set in both stages (preferably in standardized coordinates), and compare the finish with EM's own likelihood. |
| GR2-9 | info | NEWS "Growth mixtures in any units" ("to 8e-7 standard errors") | Outcome-scale equivariance holds. In a 2-class fit with membership and common covariates, though, every reported column agrees across y scales from 1e-5 to 1e6 (with shifts up to 1e8) only to 1.7e-5 to 3.7e-5 SE. That is the precision of the uncorrelated quasi-Newton finish at any scale: the Newton step from the fit is 1.6e-5 to 2.5e-5 SE at s = 1 as well (`c02_precision.R`). Not a scale defect; the stated precision depends on the data. | `c01_equivariance.R`, `c02_precision.R` | see "Verified correct" | Same remedy as GR2-4 (a Newton finish with an SE-unit criterion). |

Round 1 already noted that `.growth_permute_params` does not renormalize the
proportional scales. It is still harmless: the natural values read
`.growth_covariances()`, and permuted values equal the relabelled originals
exactly (0 difference, `c17_bootstrap_align.R`). It is a latent hazard,
though: a `pack`/`unpack` round trip of the permuted parameters changes the
likelihood by -86.1.

## Verified correct

| Algorithm | Check | Max discrepancy |
|---|---|---|
| One-class marginal ML (V = Z G Z' + sigma^2 I) | vs `lme4::lmer(REML = FALSE)` and `nlme::lme(ML)`, unbalanced, intercept + slope (`c03`) | LL 1.96e-9 (lme4), 2.4e-8 (nlme); beta 5.7e-8; G 4.5e-6; sigma^2 1.8e-6 |
| E-step BLUPs and conditional covariances (G^-1 + Z'Z/sigma^2)^-1 | vs `lme4::ranef(condVar = TRUE)` (`c03`) | means 1.2e-5, covariances 5.2e-6 (lme4's own tolerance) |
| Several classes, equal G and proportional G (`nwg`) | `lcmm::hlme` refitted from the latents optimum (`c04`) | LL 1.21e-5; coefficient SEs relative 4e-5 |
| Correlated (AR) density, BLUPs G Z'V^-1 r and conditional covariances G - G Z'V^-1 Z G | dense `mvtnorm::dmvnorm` and `solve()` (`c16`) | 7.1e-15, 4.0e-15, 1.0e-15 |
| Analytic scores: varying/equal/proportional/diagonal G, equal sigma, common + membership covariates, sampling weights, two-level clusters, AR, Brownian, AR + slope | `numDeriv::grad` of the (weighted) log likelihood at perturbed points (`c05`) | relative 2.5e-9 to 1.6e-7 |
| Per-person scores | `numDeriv::grad` of each person's log likelihood (`c03`) | 7.0e-9 absolute |
| Observed information (central differences of the score, SE-scaled step) | `numDeriv::hessian` / `numDeriv::jacobian` of the analytic score (`c03`, `c05`, `c06`) | SE relative 1.3e-8 to 7.5e-8. The diagonal-G case is 9.3e-5 against `hessian` but 6.7e-8 against the score Jacobian, so `hessian` is the less accurate there. For held correlated coordinates the free block is 2.7e-7 to 5.7e-7. |
| Robust sandwich (n/(n-1)), weighted; OPG | hand B M B with `numDeriv` bread (`c03`, `c05`, `c08`) | 3.0e-8 (unweighted), 1.8e-6 (weighted), 1.9e-5 (clusters), OPG 2e-16 |
| Delta method (variances, SDs, correlations, residual SD, shares) with the SE-scaled step | `numDeriv::jacobian` with numDeriv covariance (`c03`); `parameter_inference()` Wald vs `get_results()` (`c17`) | 1.2e-8 and 1.1e-9 relative; 1.4e-11 |
| Two-level growth: person densities, cluster likelihood sum_h pi_h prod_i sum_k p(k\|h) f_ik, rho, person posteriors | brute force with `dmvnorm` (`c08`) | 1.9e-10, 1.1e-12, 8.8e-13, 2.4e-12 |
| Two-level criteria | `bic` with log(n_clusters), `bic_persons`, `nobs` = clusters, `n_parameters` = length(coef) (`c08`) | exact |
| Outcome standardization (`.growth_standardize`/`.growth_unstandardize_fit`): every column of every `get_results(fit, "all")` table, the starts table and the history | y -> a + s y, s in {1e-5, 1e-2, 10, 1e3, 1e6}, a up to 1e8, 2 classes, membership and common covariates (`c01`) | LL 4.4e-9; BIC 8.8e-9; estimates <= 3.7e-5 SE; SEs <= 3.1e-5 relative; posteriors 2.5e-5 (convergence precision, see GR2-9) |
| The same, correlated residuals (AR, Brownian), including predictions on new persons | y -> a + s y, s in {1e-5, 1e3, 1e6} (`c07`) | LL 7.3e-12; estimates <= 1.9e-9 SE; SEs <= 1.6e-6 relative; `predict(type = "conditional")` 5.5e-10 x s; newdata posteriors 2.2e-11 |
| The same, two-level, weighted proportional, equal diagonal (`c20`) | s in {1e-5, 1e6} | <= 2.2e-5 SE (two-level); <= 1.2e-11 SE (weighted); <= 2.3e-7 SE |
| Criteria with zero weights | fit with 8 zero-weight persons vs the same fit without them: LL, BIC, bic_rows, SABIC, ICL, entropy, smallest share, robust SEs, nobs (`c09`) | identical to printed precision; SEs 3.5e-7 |
| Trajectory three-step inputs | `three_step(method = "bch")` and `"proportional"` vs a hand error matrix D, D^-1 weights and weighted means (`c16`) | 1.6e-15, 6.7e-16 |
| predict() on new data | `type = "posterior"` and `"random_effects"` with `newdata` = fitted data vs the fit's E-step (`c16`) | 0 |
| Bootstrap class alignment | `.growth_permute_params` + natural values = relabelled originals; membership logits rebased (`c17`) | 0 |
| ECM M-steps (beta given E[b]; sigma^2 with tr(Z'Z S); G = E[bb']/n_k; proportional closed forms G given s_k, s_k^2 = tr(G^-1 M_k)/(q n_k); diagonal restriction) | read against Laird & Ware (1982); EM monotonicity is asserted by `.latents_em`, and no fit above other than GR2-2 tripped it; optima agree with lcmm | — |
| Penalized-residual fallback, high variance-ratio score branch (projection correction of d_beta) | read: r'V^-1 r = \|\|r - Z m\|\|^2/sigma^2 + \|\|a\|\|^2 with a = A^-1 L'Z'r/sigma^2; X'V^-1 r = B'(Z'V^-1 r) + X_perp' r/sigma^2 | — |

## Not checked (and why)

- Plots: not rendered (formulas were read only; project rule on visual output).
- The high variance-ratio score branch numerically: round 1 checked it against central differences. No new data was built for it.
- `simulate()` moments and `.growth_draw`: round 1 checked these. They were only read again here (the factor carries the same jitter as GR2-2, which is negligible in well-conditioned coordinates).
- Bootstrap coverage: only a 6-replicate run, to check that it executes and aligns. Coverage needs a simulation study.
- VLMR/LMR/BLRT for growth (`R/mixture-regression-enumerate.R`) belong to another partition.
- GR2-7 and GR2-8 were not reproduced. Each needs a specially constructed boundary fit; both are reported from reading the code.
- No proposed fix was implemented or tested inside the package (review only). The causes of GR2-1 and GR2-2 were shown by changing `min_variance` and by a monkey-patched `.growth_factor` in a scratch session only.

## Fix status (same day)

All of GR2-1 to GR2-8 are fixed; GR2-9 is resolved by the fix and its NEWS
sentence corrected.

| Id | Fix | File |
|---|---|---|
| GR2-1, GR2-3, GR2-4 | `.growth_standardize()` also centres and scales every class-specific, common and random-design column except the intercept (centring only where an intercept absorbs it; random columns uncentred for a diagonal covariance). The fit, floors, jitter and degeneracy rule run in those coordinates; `.growth_unstandardize_fit()` maps back exactly (beta by the affine coefficient map, G = N G* N' with N = scale M^-1) and carries the E-step over (log densities shifted, BLUPs and conditional covariances by N). Inference is computed in standard coordinates and mapped by the analytic Jacobian (V = J V* J', I = J^-T I* J^-1; the Cholesky chart by the Cholesky differential). Every fit now ends with Newton steps and an SE-unit criterion (`.growth_newton_finish`, shared with the correlated path). Degeneracy: conditional variance at the floor, |r| > 0.995, or eigenvalue ratio <= 1e-6 of G* or of cov2cor(G*), on the free matrices. | `R/growth-mixture.R`, `R/growth-mixture-methods.R`, `R/growth-mixture-correlated.R` |
| GR2-2 | `.growth_factor` is exact when cov2cor(G) has smallest eigenvalue > 1e-10, else jitters each diagonal by 1e-12 of itself | `R/growth-mixture-engine.R` |
| GR2-5 | class table counts positive-weight persons only | `R/growth-mixture-methods.R` |
| GR2-6 | multilevel ICL = BIC + 2 x joint (cluster-class and person-class) entropy | `R/growth-mixture-methods.R` |
| GR2-7 | residual-SD row held per class | `R/growth-mixture-parameter-inference.R` |
| GR2-8 | EM floors each conditional variance (squared Cholesky diagonal) at `min_variance`, the same bound the quasi-Newton and Newton stages use | `R/growth-mixture-engine.R` |
| prediction / bootstrap matching | E-step of new persons computed in the fit's standard coordinates; standard parameters re-derived by the exact inverse map when `fit$params` was changed after fitting | `R/growth-mixture.R`, `-predict.R`, `-parameter-inference.R` |

After the fix (scripts `c01`, `c02`, `f01`, `f03`): outcome scale 1e-5 to 1e6
with shifts to 1e8, every table <= 5.5e-12 SE, SEs <= 1.4e-9 relative,
distance from the root <= 1e-11 SE; time x 365, 8760, 1e5 and +2000:
coefficients <= 9e-14 SE, LL <= 2.3e-13, same flags; multilevel, proportional,
diagonal (scale only) <= 3e-13 SE. Tests: `tests/testthat/test-review-1026-growth-r2.R`
(9 tests, 71 expectations); each of the 9 fails on the pre-fix code (old
definitions restored over the new ones in a scratch copy, `round2-growth/oldcopy/R/zzz-old-growth.R`).

Follow-up: `test-engine-fixes.R` ("a growth mixture is not trapped where the
level-only starts lead") failed after the fix (LL -1096.4 against -1085.0
from the truth). Cause: `.growth_trajectory_seed()` grouped the persons'
coefficients in the standard coordinates, where each person's level is read
at the mean of time; there the two slope-defined classes have the same
level (54.5 against 55), so the level column was pure noise and Ward's split
matched the true classes no better than chance (34/36 against 26/24). The
seed now reads levels at each predictor's minimum (baseline). This equals the
old seed whenever time starts at 0, and it stays free of the predictors'
units and origin. The seed recovers 57/1 against 3/59, and the fit reaches
-1085.03 with slopes 0.017 and 3.006.
