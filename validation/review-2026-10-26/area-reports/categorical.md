# Review 2026-10-26: categorical partition

## Scope

Files: `R/categorical.R`, `R/extra-indicators.R`, `R/local-dependence.R`,
`R/direct-effects.R`, `R/bivariate-residuals.R`, `R/multilca.R`, `R/weights.R`.
In commit 8ac3576 Codex changed only `R/extra-indicators.R` (native numeric
ordinal values, Poisson log-ratio offset) and `R/weights.R` (relative
within-unit weight check). I read every algorithm in the seven files, not
only the diff. All checks were run under `devtools::load_all()` with
`OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1` (R 4.5.2, macOS arm64). poLCA,
numDeriv, Rmpfr and MASS are installed.

## Verdict on each Codex change

| Change | Verdict | Evidence |
|---|---|---|
| `.latents_extra_log_density(offset = TRUE)`: Poisson log ratio `y log(mu_c/mu_ref) - (mu_c - mu_ref)` plus a per-row offset equal to the reference class's `dpois` | **Correct** | Algebra checked: the `lgamma(y+1)` term cancels, so the expression equals the log-pmf difference. Numerically, on 50 rows (3 classes, ordinal plus Poisson, NAs): offset + centred density vs plain: max abs diff 3.6e-15; plain vs hand-computed: 2.2e-15; posteriors from offset vs plain: 3.3e-16. A full mixed fit (Gaussian + binary + ordinal + Poisson) reproduces hand-computed posteriors to 3.9e-16 and the log likelihood to 1.1e-13. |
| Same, accuracy at huge counts | Correct, with a trade-off (low, informational) | Against Rmpfr (200 bits): at y = 1e17 with means 1e17 and 1e17(1+1e-9) the offset gives -0.0500000119 against exact -0.0500000160 (abs err 7.9e-9); plain `dpois` differences are better there (err 1.5e-15) because `dpois` uses the saddle-point `bd0` form. Where every log pmf is about -1e16 (y = 1e17+160, means 3e16 and 3e16+24) the offset is exact (56) and plain `dpois` is wrong (80). Both regimes need counts near 1e17, so ordinary fits are unaffected. |
| Native numeric ordinal values (`ordinal_values`, 17-digit labels when `as.character` collides) in fitting and `newdata` | **Correct but incomplete** (fixed) | Values 1e17 + {0,16,32}, which all print as "1e+17": three categories, log likelihood identical to a fit on codes 1:3 (diff 0), and `predict(newdata)` identical. **Incomplete:** `.latents_as_column_type()` (used by the LTA parametric draw in `R/lta-compare.R`) mapped labels back with `as.numeric(as.character())`. On this platform R's string-to-number conversion is not exact at 18 digits: `as.numeric("100000000000000016") - (1e17+16)` = -16, and `"100000000000000032"` gives 1e17+48. So a redrawn column merged two categories and created a value that did not exist (finding F3). |
| `R/weights.R`: relative within-unit spread `1 - min/max > sqrt(eps)` | **Correct** | Units with weights 1e-3 and 1.0001e-3 are now refused (`latents_bad_weights`); the old absolute floor accepted them. Rounding-only differences (0.1 and 0.3-0.2) are accepted and normalized to 0.375, 0.75, 1.875 (sum 3). A unit of {0, 1e-300} is refused, as it should be. Weights 1e300/1e308/1 normalize without overflow. Scaling the weights by 1e-6 changes the log likelihood by 3.4e-13 and the probabilities by 5.9e-14. With integer weights, weighted LL x sum(w)/n equals the LL on replicated rows (diff 1.3e-10). |
| Constancy-test pattern (sd <= 0) | Not present in partition | The only `sd() <= 0` test (`.latents_extra_clustering`) uses `stats::sd`, which returns exactly 0 for exact constants (checked on 3,000 values x 39 lengths: max sd 0). The real instance of the pattern in my partition is a different one (F1): a covariate that is exactly constant on a row subset, after a weighted-mean centring, leaves 4.2e-17 noise. |

## Findings

| id | severity | file:line | problem | reproduction (actual vs expected) | status |
|---|---|---|---|---|---|
| F1 | **high** | R/direct-effects.R:117-137 (`.latents_center_direct`) | A direct effect is identified only by covariate variation among rows that observe its target with a positive weight. `.latents_prepare_direct` checks constancy over all rows only. If the covariate is constant where the target is observed, the centred column is rounding noise (4.163e-17), and `qr.solve` returns an arbitrary effect. The uncentred means are then garbage, and the fit still reports `converged = TRUE` with no warning. | n=200, `y1` NA in rows 1:30, `z = 0.1` except rnorm in rows 1:30; `lpa(..., direct_effects = "y1 ~ z", missing = "fiml")`. Actual: `y1:z = -3.635928e+13`, means of `y1` 3.635928e+12 in both profiles (the same LL, -618.8746, as the fit without the effect). Expected: a refusal. | **fixed**: raises `latents_bad_data` when the covariate has fewer than two distinct values in those rows |
| F2 | medium | R/local-dependence.R:408-412 (`.latents_local_dependence_table`) | A joint category that nobody gave (it is common with sparse items) is in the fit's levels but not in `local_dependence_categories`, which is read from the data. Its rows get NA `category_1`/`category_2`/`independent`/`ratio`, and its floor mass is dropped from the marginals. | a,c pair with (a=1,c=3) never observed: `get_results(fit, "local_dependence")` rows 6 and 12 had `category_1 = NA`, `independent = NA`. Expected "1","3" and the product of the marginals. After the fix, each profile's `independent` sums to 1 (tolerance 1e-12). | **fixed**: unseen labels are split into their two halves (colon-free categories are guaranteed at join time) |
| F3 | low | R/extra-indicators.R:983 (`.latents_as_column_type`) | Incomplete native-value support (see the verdict table): the label-to-number round trip is inexact for 17-digit labels. | `.latents_as_column_type(drawn, 1e17+c(0,16,32)[code])` returned values that were not identical to the input. Expected: identical. | **fixed**: labels map back through the column's own sorted values via a shared `.latents_ordinal_labels()` helper. If the labels do not match, the old behaviour is kept. |
| F4 | low (informational) | R/extra-indicators.R:203-219 | Poisson offset is less accurate than plain `dpois` differences when y is about 1e17 and the means are near y (7.9e-9 vs 1.5e-15 abs). It is more accurate when all log pmfs are about 1e16. | See the verdict table. | reported-not-fixed (judgement call; negligible statistically) |
| F5 | low (informational) | R/bivariate-residuals.R:165-193 | The categorical BVR is the Pearson X^2 of the observed table against the model-implied one, which is Latent GOLD's numerator (LG prints X^2/df). The p-value uses chi^2(df), which ignores parameter estimation (Oberski, van Kollenburg & Vermunt 2013). The doc says "approximate". | Hand-computed Pearson tables equal the function's statistic exactly (diff 0) for overall and by-profile. | verified-correct; caveat reported |
| F6 | low | R/categorical.R:94-96 | Character indicators are ordered by `sort()`, which depends on the locale's collation. Category codes and labels can therefore differ between locales. The likelihood does not. | not reproduced across locales | reported-not-fixed |

## Algorithms verified correct

| Algorithm | Check | Max abs difference |
|---|---|---|
| Categorical LCA E/M (unrestricted multinomial) | poLCA, 400 rows, 4 items, 2 classes, tight tol | LL 3.1e-11; probabilities 8.3e-7; npar 13 = 13 |
| LCA E-step | hand Bayes posteriors from fitted parameters | < 1e-12 |
| LCA M-step fixed point | one closed-form M-step from fitted posteriors | about 1e-8 (fit stopped at tol 1e-12) |
| Bounded multinomial (`min_probability`) | KKT: common multiplier on free categories, bound-category gradient <= multiplier | multiplier spread 0; closed form matches 1e-15 |
| Categorical parameter count C·sum(K-1); joint item C·(K1K2-1) | poLCA npar; LD fit 15 = 1 + 2·(5+1+1) | exact |
| Adjacent-category ordinal (LG default): log P(k)/P(k-1) = a_k - a_{k-1} + eta_c | direct identity; rows sum to 1 | 2.2e-16 / 4.4e-16 |
| Ordinal Newton M-step | numDeriv gradient at solution; BFGS optimum | grad 1.2e-9; objective +5.4e-13 over optim |
| NB2 log-dispersion score and curvature (regular, small-alpha series, large-size asymptotic branches) | Rmpfr 300-bit exact derivatives on a 140-point grid (y 0..200, mu 0.5..300, alpha 1e-9..5) | rel 2.0e-10 (score), 9.8e-11 (curvature) |
| NB Hessian (varying and shared dispersion) | numDeriv | 5.9e-9 (scale 74), 2.8e-9 |
| NB modified-Newton M-step | BFGS on the same objective | objective +3.4e-13 over optim |
| NB Poisson-limit KKT (score at alpha=0 is sum w((y-mu)^2 - y)/2) | under-dispersed data | dispersion = 1e-8 floor, mean = 4 = mean(y) |
| Extra-indicator scores (Fisher identity) for Poisson, NB varying, NB equal | numDeriv gradient of observed mixture LL | 2.5e-7, 1.1e-6, 1.7e-7 (scales 288, 98, 199) |
| Extra parameter counts | widths 8/12/10 = (K-1)+(C-1) + C means + C or 1 dispersions | exact |
| Direct-effect ECM M-step | weighted normal equations solved directly; brute-force optim | 1e-12 (test); scores at solution 2.3e-14 |
| Categorical bivariate residuals | hand Pearson tables (overall and by-profile) | 0 |
| Sampling weights | scale invariance; replicated-rows equivalence | 3.4e-13; 1.3e-10 |
| Categorical log density | label-permutation equivariance | < 1e-15 |

## Tests added

`tests/testthat/test-review-1026-categorical.R` (17 `test_that` blocks):
Poisson offset vs hand log pmf; Poisson offset vs Rmpfr at 1e17; ordinal
identity and Newton M-step; NB derivatives vs Rmpfr; NB Hessian/M-step vs
numDeriv/optim; extra scores vs numDeriv; mixed fit vs hand posteriors;
LCA vs poLCA plus EM fixed point; label-permutation equivariance; bounded
shares KKT; BVR vs hand; weights relative check and scale/replication
invariance; direct-effect normal equations; **regressions** F1 (direct-effect
refusal), F2 (unseen joint category) and F3 (huge ordinal round trip).

Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 57 ]` (about 3 s).

I temporarily reverted the three fixes by hand. The result was
`[ FAIL 6 | WARN 0 | SKIP 0 | PASS 51 ]`: the failures were at lines 360
(F1), 381/382/384/389 (F2) and 411 (F3). After I restored the fixes the
file passes again.

poLCA and Rmpfr are not in Suggests. They are reached through
`getExportedValue()` behind `skip_if_not_installed()`, so `R CMD check`
does not flag an unstated test dependency. If the lead prefers, add them
to Suggests and switch to `::`.

## Related existing test files

I ran these with `NOT_CRAN=true`, so the `skip_on_cran` blocks also ran. Each count is test blocks, then failed/errors/skipped:

- test-bivariate-residuals 9/0/0/0, test-categorical-inference 6/0/0/0, test-categorical-suspect 5/0/0/0, test-categorical 11/0/0/0
- test-direct-effects 10/0/0/0, test-features-direct-effects 12/0/0/0, test-extra-indicators 9/0/0/0, test-local-dependence 11/0/0/0 (re-run after the final edit)
- test-multilca 4/0/0/0, test-negative-binomial-regression 4/0/0/0, test-ordinal-regression 6/0/0/0, test-review-extra 5/0/0/0
- test-review-1020-zero-weights 2/0/0/0, test-weighted-logits 1/0/0/0, test-weights 14/0/0/0, test-zero-weight-initialization 6/0/0/0
- Other files that touch direct effects or local dependence: test-features-joint-items 11/0/0/0, test-features-pair-covariance 9/0/0/0, test-lmr 16/0/0/0, test-parallel-profile-pool 12/0/0/4 (own skips), test-review-0970-estimation 4/0/0/0, test-review-1020-inference 6/0/0/0, test-review-1020-integration 1/0/0/0, test-review-pairs-trees 7/0/0/0, test-thin-profiles 8/0/0/0

No regressions.

## Files edited (uncommitted)

- `R/direct-effects.R`: F1 refusal in `.latents_center_direct`
- `R/local-dependence.R`: F2 labels for unseen joint categories (`sub()` on the colon-free halves, so an empty category is safe)
- `R/extra-indicators.R`: F3. Codex's label rule is factored into `.latents_ordinal_labels()` (fitting output unchanged), and `.latents_as_column_type()` maps labels through native values
- `tests/testthat/test-review-1026-categorical.R`: new file

## Not checked, and why

- No Latent GOLD run for ordinal, NB or BVR values: LG under Wine is outside a fast single-core review. The formulas were checked against their definitions and exact arithmetic instead.
- The LTA parametric-draw path (`R/lta-compare.R`) was tested only at the helper level for F3. That file is outside my partition.
- Gaussian-pair residuals (Fisher-z statistic): the code was read but not calibrated.
- Weighted design-based BVR inference is not implemented (NA by design).
- The cross-locale ordering of character categories (F6) was not run.
