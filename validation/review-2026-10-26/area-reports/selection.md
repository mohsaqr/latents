# Review 2026-10-26 — partition `selection`

## Scope

Files: R/three-step.R, R/three-step-ml.R, R/diagnostics.R, R/lmr.R,
R/compare-models.R, R/enumeration.R (BLRT, grid criteria, AvePP), R/candidates.R,
R/search.R, R/selection.R, R/descriptives.R. Codex commit 8ac3576 touched six of
them (compare-models, descriptives, search, selection, three-step-ml, three-step).
All checks run under `devtools::load_all()` with OMP/VECLIB threads = 1 and
`NOT_CRAN=true`. Other reviewers were editing other R files concurrently, so
the test runs below load their in-progress edits too.

## Verdict on each Codex change

| File | Change | Verdict | Evidence |
|---|---|---|---|
| three-step-ml.R | Equal-variance BFGS candidates in the chart variance = floor + exp(2a) | **Correct math; error handling was a defect (fixed)** | Chain rule: d log sd / da = exp(2a)/variance = `exp(2*point[k+1] - 2*decoded[k+1])` — matches line 176. Start `log(1 - smallest)/2` gives variance exactly 1, and `values` *is* standardized (lines 38-40: `(values - mean)/sd`), so variance 1 is the outcome's own variance — appropriate. Independent ML3 (optim + numDeriv) agrees with `three_step(method = "ml", variance = "equal")`: max abs mean diff 8.8e-8, SE relative diff 6.6e-8. The blanket `tryCatch(..., error = function(condition) NULL)` around the whole candidate was a silent-failure violation; narrowed (see F1). |
| three-step.R | Roxygen text on BFGS, Box-Muller cache | Correct | Doc only; matches the code. |
| compare-models.R | `.latents_same_comparison_data()`: exact (tolerance 0) data identity, 1e-12 on normalized weights, forward-map proof for centered LTA storage | **Correct but incomplete** | `all.equal(tolerance = 0)` is exact (mean relative difference must be <= 0). Forward map verified: for homogeneous vs occasion LTA fits, `sweep(raw, centers, "-")` equals the general engine's `x` with max diff 0. But the identity still failed on matrix row names (F2, pre-existing, also failed with the old `all.equal` default), so the forward-map path was never reached in that real pairing. Fixed. |
| descriptives.R | `.multilpa_icc()` uses integer indices | Correct | ICC(1) = (MSB − MSW)/(MSB + (n0 − 1)MSW), n0 = (N − Σn_j²/N)/(k − 1) vs `anova(lm())` on 8 unbalanced groups (sizes 3–20): diff 5.6e-17. Native ids 0.3 vs 0.1+0.2 are now two groups (ICC 0.17989 = hand ANOVA on the two groups). Row-order and relabel invariance diff 0. `split(values, index)` returns groups in index order 1..k, aligned with `tabulate(index)`. |
| search.R | Softened claim that shared-covariance likelihood maxima are never degenerate | Correct | Doc only; the new wording (shared covariance removes the isolated-component collapse, floor/boundary checks still needed) is accurate (McLachlan & Peel 2000 §3.8). |
| selection.R | Trust view: `NA` replication → "replication unavailable", blank label | Correct logic; low style issue | Status ordering is right. The new fill `#BDBDBD` is not an Okabe-Ito colour (project rule); the state is also distinguished by its blank label, so not colour-only (F4, reported). |

Same-pattern search (lead's `sd <= 0` constancy regression): the only variance
test in this partition is `stats::var(values) <= 0` (three-step-ml.R:30).
`var()` uses R's two-pass mean, and I computed `var(rep(x, n))` for
x ∈ {0.1, 0.3, 0.7, 1/3, 2/3, 0.2, 1e8+0.1, 123.456} × n ∈ {3,5,7,10,50,100,333,1000}:
all 64 are exactly 0. `three_step(..., method = "ml")` on y = 0.1 for every row
raises `latents_bad_outcome`; BCH returns 0.1 exactly with SE 0. Other
constancy tests here are exact (`diff(range()) <= 0` in selection.R:113,
`length(unique(v)) == 1L` in three-step.R). No instance of the pattern.

## Findings

| id | severity | file:line | problem | reproduction (actual vs expected) | status |
|---|---|---|---|---|---|
| F1 | medium | R/three-step-ml.R:179-213 (was 179-199) | Blanket `tryCatch(..., error = function(condition) NULL)` around every BFGS candidate: any error (a bug in scores, Hessian, `slim()`) silently dropped the candidate. | Probe: the handler recorded every condition across all three-step test files (0 caught) and 60 ill-posed random cases (4 caught, all "the leading minor of order 2 is not positive" from `chol()`). Regression test injects a classed error into the first `.multilpa_step_ml_log_likelihood()` call (a BFGS start): before the fix the fit returns silently (test fails: 1 failure); after, `latents_test_injected` reaches the caller. | **fixed**: start checked with `is.finite(objective(start))` (the only `optim` BFGS error is a non-finite initial value), non-finite information rejected explicitly, and only `chol()` on a finite symmetric matrix is caught (its sole failure is non-positive-definiteness). 61 cases (60 random + rare-class fixture) give `identical()` fit lists before and after. |
| F2 | medium | R/compare-models.R:521-545 (`.latents_comparison_signature`) | Data identity compares `indicators` including row names. A homogeneous `lta()` fit (`multilpa_transitions`) keeps the frame's row names on its indicators, the general engine's `x` has none, so `compare_models()` refuses two fits of the *same* data as "fitted to different data". Pre-existing (old `all.equal()` default also FALSE) — not a Codex regression, but it defeats Codex's forward-map path. | `course_engagement`, students <= 40 (and <= 12), `lta(..., n_profiles = 2)` vs same with `transitions = "occasion"`: `compare_models(h, o)` → `latents_incomparable_models`; field-by-field diff shows only `dimnames[[1]]` ("1".."547" vs NULL); forward-map diff 0. | **fixed**: row names dropped from `raw` and `centered` in the signature. Regression test fails before (error) and passes after; a 1e-9 change to one value is still refused. |
| F3 | low | R/diagnostics.R:89 | ICL is the *soft* ICL-BIC (BIC + 2·EN, McLachlan & Peel 2000), while mclust's `icl()` uses the hard MAP entropy (Biernacki et al. 2000). Documented as "2 EN", so not a bug, but numbers differ from mclust. | VVI, 3 profiles, n = 240: latents ICL (mclust sign) −1891.17 (soft) vs hard version −1854.13 vs mclust −1853.73 (mclust's own loglik −869.994 vs latents −869.989). | reported-not-fixed (convention) |
| F4 | low | R/selection.R:410 | New trust-view fill `#BDBDBD` is not Okabe-Ito (project palette rule). Not colour-only (blank label). | n/a (visual) | reported-not-fixed (plot change needs rendered review) |
| F5 | low | R/three-step.R `.multilpa_r3step_frame` | `sqrt(pmax(diag(covariance), 0))` turns a negative variance into SE 0 → statistic/p `NA` without a warning. | Not triggered on test data; observed by reading. | reported-not-fixed |

## Algorithms verified correct (computed checks)

| Algorithm | Check | Max abs difference |
|---|---|---|
| AIC, KIC, BIC, SABIC, CAIC, AWE, ICL, CLC (individual convention) | Hand formulas from fitted L, q, n and EN = −Σ p log p of the posteriors (3-profile VVI, n = 240) | 0 |
| BIC sign/parameter count | `mclust::bic("VVI", latents L, n, d = 2, G = 3)` = −latents BIC; `nMclustParams` = q = 14 | 0 (test tol 1e-12) |
| Relative entropy 1 − EN/(n log K) (Ramaswamy et al. 1993) | hand | 0 |
| AvePP (Nagin 2005) | mean posterior within modal class | 6.7e-16 |
| OCC = [AvePP/(1−AvePP)]/[π/(1−π)], π = mean posterior | hand | 3.8e-12 |
| Label/row permutation invariance | entropy, relative entropy, AvePP matrix (permutes as P[order, order]), OCC | 7.1e-15 |
| BCH (Bolck, Croon & Hagenaars 2004; Vermunt 2010) | D[t,s] = Σ p_it 1(W_i=s)/Σ p_it; weights = D⁻¹[W_i, ]; class means; influence wᵢₜ(yᵢ − μₜ)/Σw; group-summed CR1 SE; pairwise SE; Kish effective n; stored `bch_weights` | est 5.3e-15, SE 5.6e-17 (cluster), 1.1e-16 (independent) |
| Proportional / modal | posterior- and modal-weighted means | 1.8e-15 / 1.1e-16 |
| ML step three, varying variances (Vermunt 2010; Bakk, Tekle & Vermunt 2013) | Independent likelihood Σ log Σ_t π_t D[t,Wᵢ] φ(y; μ_t, σ_t), BFGS→Nelder-Mead→BFGS, numDeriv Hessian/Jacobian sandwich, clustered CR1 | means 4.8e-8, SE rel 6.2e-7, variances 5.0e-8 |
| ML step three, equal variance | same | means 8.8e-8, SE rel 6.6e-8, variances 2.3e-7 |
| ICC(1) with n0 (unbalanced) | `anova(lm(y ~ factor(g)))` mean squares | 5.6e-17 |
| Imhof tail (Imhof 1961) | CompQuadForm::farebrother (positive weights, 3 cases); conditioning integrals for (1, −0.3) at q = 1 and q = −0.5, (0.7, −0.2, 0.05) at q = 0.4; `pchisq(3, 3)` | 1e-15 (farebrother), 3.6e-14, 0 (q = −0.5), 2.0e-14, 5.6e-17 |
| VLMR weights (Vuong 1989, Thm 3.3) | Independent W = [−B₁A₁⁻¹, −B₁₀A₀⁻¹; B₁₀′A₁⁻¹, B₀A₀⁻¹] from numDeriv scores/Hessians in my own chart (means, log sd, logits), univariate 2 vs 3 profiles VVI, n = 300 | reference mean 4.659535 vs 4.659539, sd 5.313557 vs 5.313562, VLMR p 0.4342634 vs 0.4342636, adjusted p 0.4540150 vs 0.4540152 (numDeriv-limited) |
| LMR adjusted LR | LR/(1 + 1/(d log n)), d = parameter difference; tidyLPA uses d = 3Δk (= parameter difference for univariate varying-variance mixtures, LMR 2001 eq. 15) | exact (test tol 1e-12) |
| BLRT p-value | (1 + #{T* ≥ T})/(B + 1) (Davison & Hinkley 1997), ties counted, NA + `latents_failed_replicates` when any replicate invalid | exact |
| Kneedle chord distance | hand perpendicular distance in rescaled units; invariant to affine rescaling of the criterion | 1.4e-14 (test tol) |
| r3step score | read: ∂/∂β_s log Σ_t π_t(z) D[t,W] = (posterior_s − π_s) z — correct | by derivation only |

## Tests added

`tests/testthat/test-review-1026-selection.R`, 13 tests / 60 expectations:
information criteria vs definitions and mclust; entropy/AvePP/OCC vs hand;
relabel invariance; BCH means/SE/effective n/stored weights vs hand; ML step
three vs independent likelihood + sandwich (skip_on_cran, numDeriv);
unexpected BFGS error not swallowed (F1 regression); constant outcome;
ICC vs ANOVA; Imhof vs independent integrals; VLMR weights vs independent W
(skip_on_cran, numDeriv); BLRT p-value; knee; comparison identity ignores
row names (F2 regression).

Result (`testthat::test_file()`, NOT_CRAN=true): **13 tests, 60 expectations,
0 failed, 0 errors, 0 skipped, ~4.6 s.** With both fixes reverted by hand:
test 6 failed (1 failure) and test 13 errored, all others passed.

## Related existing tests run (after the fixes)

All with `NOT_CRAN=true`, `testthat::test_file()` under `load_all()`:

| file | tests | expectations | failed | errors |
|---|---|---|---|---|
| test-three-step-ml.R | 10 | 69 | 0 | 0 |
| test-audit-ml3-optimizer.R | 3 | 29 | 0 | 0 |
| test-features-ml-three-step.R | 17 | 85 | 0 | 0 |
| test-three-step.R | 16 | 74 | 0 | 0 |
| test-tidy-three-step.R | 15 | 77 | 0 | 0 |
| test-growth-r3step.R | 4 | 16 | 0 | 0 |
| test-three-step-clusters.R | 8 | 41 | 0 | 0 |
| test-lta-three-step.R | 4 | 22 | 0 | 0 |
| test-diagnostics.R | 13 | 133 | 0 | 0 |
| test-lmr.R | 16 | 134 | 0 | 0 |
| test-features-vlmr.R | 11 | 176 | 0 | 0 |
| test-descriptives-safety.R | 9 | 37 | 0 | 0 |
| test-tidy-diagnostics.R | 11 | 56 | 0 | 0 |
| test-gaussian-diagnostics.R | 7 | 41 | 0 | 0 |
| test-search.R | 17 | 232 | 0 | **1** |
| test-default-search.R | 7 | 28 | 0 | 0 |
| test-enumeration.R | 27 | 321 | 0 | 0 |
| test-enumerate-verbs.R | 7 | 30 | 0 | 0 |
| test-enumeration-type.R | 9 | 122 | 0 | 0 |
| test-lta-compare.R | 5 | 31 | 0 | 0 |
| test-additive-compare.R | 4 | 29 | 0 | 0 |
| test-parallel-candidates.R | 1 | 16 | 0 | 0 |
| test-audit-inference-repairs.R | 17 | 85 | 0 | 0 |
| test-audit-recent-api.R | 8 | 42 | 0 | 0 |
| test-informative-rows.R | 13 | 57 | 0 | 0 |
| test-post-submission-review.R | 17 | 60 | 0 | 0 |

The one error (test-search.R:175, "emEM gives the same fit on one core and on
two") is environmental: its `cores = 2` PSOCK workers load the *installed*
latents 0.95.1, which lacks `.latents_collect_conditions` (defined in source
R/parallel.R, version 0.97.1): "could not find function
\".latents_collect_conditions\"". Unrelated to this partition's edits;
reinstalling the package would be needed to run it.

Before the edits, the same three-step files (ml, audit-ml3, features-ml,
three-step, tidy, growth-r3step, clusters) passed with identical counts under
the probe that instrumented the old handler.

## Could not check / not done

- Latent GOLD / Mplus were not run; the VLMR Latent GOLD calibration in
  tests/testthat/test-lmr.R (printed values to 5e-5) passed.
- BLRT end-to-end reference distribution (simulation size) not re-validated;
  only the p-value arithmetic and replicate handling.
- Sampling-weighted entropy/AvePP conventions not compared with an external
  program (entropy uses unweighted posteriors; class sizes are weighted).
- No `R CMD check`; no plot rendering (F4 left unchanged for that reason).
