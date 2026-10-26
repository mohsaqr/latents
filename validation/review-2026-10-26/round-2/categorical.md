# Round 2 review: categorical and extra indicators

All scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-categorical/`
(written `round2/` below). They were run with `pkgload::load_all()`,
`NOT_CRAN=true OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`, R 4.5.2 on macOS arm64.
Nothing under `R/`, `tests/`, `equivalence/` or `man/` was edited. The only
patch (`nb_shared_patched.R`) was applied to the in-session namespace.

## Scope (files)

`R/categorical.R`, `R/extra-indicators.R`, `R/local-dependence.R`,
`R/direct-effects.R`, `R/bivariate-residuals.R` and
`R/structure-inference.R`. I also read the call sites that decide what these
files compute: `R/gaussian-moments.R:75-165` (the covariance-pattern M-step),
`R/inference.R:315-360` (free index of a pattern), the docs in
`R/get-data.R`, `R/report.R` and `R/multilpa.R`, and the BVR passages of the
vignettes. Round 1 (`area-reports/categorical.md`) is not repeated here,
except where these checks add a measurement to one of its findings.

A note on the brief: the ordinal block is **not** a cumulative logit. It is
the adjacent-category logit, documented as such (`R/multilpa.R:1227-1236`)
and matching Latent GOLD's default ordinal model, so `polr` is not a valid
reference. I checked it against the equivalent Poisson log-linear
row-effects model (`glm`) instead.

## Findings

| id | sev. | file:line | what the code does / what it should do | reproduction | measured effect | proposed fix |
|---|---|---|---|---|---|---|
| C2-1 | **medium** | `R/extra-indicators.R:1265-1274` (`.latents_poisson_limit`), reached from `.latents_negative_binomial_maximize` (l.1248) | With `count_dispersion = "equal"` (shared NB dispersion), the shared alpha is set to the 1e-8 floor whenever the **summed** score at alpha = 0 is <= 0. The objective is never compared. For a single profile, a non-positive score at 0 does imply a global maximum at 0 (the Anscombe/Aragon condition; no counterexample in 4,000 random weighted cases, `nb_varying_weighted_scan.R`). The summed rule across profiles is only a **local** condition. One overdispersed profile and one under-dispersed profile give a bimodal profile likelihood: the floor is a local maximum and the interior maximum is higher. This is the defect fixed in round 1 as MR-F4 for `mixture-regression-engine.R:425-441`. The extra-indicator copy was not fixed. | `round2/nb_shared_repro.R` (M-step alone), `round2/nb_shared_fit2.R` (user API), `round2/nb_shared_scan.R` | M-step started **at** the interior maximum (alpha = 1.3228, objective -570.78) returns alpha = 1e-8 with objective -631.29, a loss of **60.5 log-likelihood units**. Through `lpa(..., count = "y", count_model = "negative_binomial", count_dispersion = "equal")` on well-separated data, every start fails: `latents_all_starts_failed: EM likelihood decreased beyond numerical roundoff`. The same data with `"varying"` fit fine (LL -937.42). With the fix below, the `"equal"` fit converges to alpha = 1.32278, which equals the hand profile-likelihood maximum (1.322783) on the same (separated) partition. The likelihood is at least 60 units better than at the floor. The same M-step serves `lta()` and the two-level fits. | Do what MR-F4 does: compute each class's objective (the summed objective when shared) at the Newton point and at the KKT boundary candidate, and move to the floor only if `gain >= -1e-12 (1 + |objective|)`. The in-session patch in `round2/nb_shared_patched.R` is this rule, about 12 lines. Add a regression test with the two-profile data from `nb_shared_repro.R`. |
| C2-2 | **medium** | `R/local-dependence.R:208-233` (`.latents_join_pairs`): the joint label is `NA` unless **both** items are observed | Under `missing = "fiml"` a row that observes one item of a categorical pair loses that item. The joint indicator is missing, so the row's likelihood leaves out `P(item1 = a \| c) = sum_b P(a, b \| c)`. `?multilpa` says `"fiml"` "maximizes the observed-data likelihood under ... MAR". This is not that likelihood. When the missingness of item 2 depends on the observed item 1 (MAR), deleting item 1 makes the pair's missingness depend on a value the likelihood no longer sees, so the estimator is inconsistent. A `latents_partial_pair` warning is raised, but it says only that "the pair is treated as missing", not that the estimates are biased. | `round2/partial_pair_bias.R` | n = 40,000, 2 classes; u1:u2 dependent; u2 missing with probability 0.7 if u1 = 1, else 0.05 (MAR). Estimates: class shares 0.640/0.360 against truth 0.6/0.4 (the complete-data fit gives 0.599/0.401). P(u1 = 1 \| class) **0.068 / 0.521** against truth 0.2 / 0.8. At this n the sampling SE is about 0.005, so the bias is more than 50 SE. | Integrate out the missing half rather than dropping the row. In the E-step a half-observed row contributes `log sum_b P(a, b \| c)`. In the M-step the expected joint counts spread that row's posterior over `(a, b)` in proportion to `P(a, b \| c) / P(a \| c)`. This is closed-form EM, the same as the missing-cell treatment elsewhere. The score (Fisher identity) changes in the same way. Until then, reword the warning to say that under MAR this biases the estimates. |
| C2-3 | **medium** | `R/bivariate-residuals.R:22-23, 41-46` | `get_results(fit, "residuals")` computes residuals only for Gaussian-Gaussian and categorical-categorical pairs. Every pair involving an ordinal or count indicator, and every Gaussian-categorical pair, is dropped silently. The contract says "Bivariate residuals for **every indicator pair**" (`R/get-data.R:103`), and `diagnostics()` presents the table as the check of "the within-profile independence the model assumes" (`R/report.R:8`). | `round2/bvr_coverage.R` | 6 indicators (2 Gaussian, 2 binary, 1 ordinal, 1 count) give 15 pairs. **2** are reported. No message or condition says that 13 were skipped. | Either implement the missing kinds (ordinal and count pairs as categorical tables of the model-implied joint, the same Pearson table, as Latent GOLD does) or state the coverage in `?get_results` and add a `message()` listing the pairs not assessed. |
| C2-4 | medium (measured update of round-1 F5, listed as "approximate") | `R/bivariate-residuals.R:182-199` | The categorical BVR p-value refers Pearson X² of the observed table against the model-implied table to chi²((R-1)(C-1)). It ignores that the expected table is fitted from the same data (Oberski, van Kollenburg & Vermunt 2013), so the test is very conservative. The LCA vignette filters pairs on it (`subset(get_results(fit, "residuals"), p_adjusted < 0.01)`, `vignettes/articles/workflow-lca.Rmd:142`). | `round2/bvr_calibration.R` (200 replications, n = 500, data from a true 2-class local-independence model) | Categorical, `by = "profile"`: rejection rate **0.25%** at nominal 5% (0.04% at 1%), mean p 0.685. `by = "overall"`: **0.08%** at 5%, mean p 0.746. The Gaussian Fisher-z test on the same design is calibrated: 4.8% (profile) and 5.8% (overall) at 5%, mean p 0.52/0.50. | Replace the p-value with a parametric bootstrap of X², or with the score test of Oberski et al. (2013) using the fit's information. At a minimum, document that the chi² reference is conservative, by about a factor of 20 at the 5% level in this design. |
| C2-5 | low | `R/bivariate-residuals.R:3-4` vs `R/get-data.R:103-106, 496-506` | The roxygen says that "the statistic per pair kind, what the test assumes and the references are documented on `?get_results`". They are not. `?get_results` does not say that `observed`/`expected` are correlations for Gaussian pairs and Cramer's V for categorical pairs, that `statistic` is a Fisher-z difference times sqrt(n_eff - 3) or a Pearson X², or what the reference distributions are, and it cites no reference. | `grep -i "fisher\|cram\|pearson" man/get_results.Rd` finds nothing about residuals | documentation only | Add the two statistics, their nulls and the references (Fisher 1921; Cramer 1946; Oberski et al. 2013) to the *Bivariate residuals* section. |
| C2-6 | low | `R/multilpa.R:1081-1083` vs `R/local-dependence.R:236-239` | The docs say that the joint item's categories are "the **observed** combinations of the two items' categories". The code uses the **full grid** on purpose (comment at `R/categorical.R:34-35`), so an unobserved combination gets a probability on the floor and counts as a parameter. Separately, l.1091-1092 says that a pair "not both categorical" is refused, but continuous pairs are accepted; the refused case is a *mixed* pair. | `round2/joint_grid.R` | 3 of 4 combinations observed: levels `0:0 0:1 1:0 1:1`, npar 11 (full grid) rather than the 9 the documentation implies. BIC therefore carries 2 parameters the reader does not expect. | Reword to "every combination of the two items' categories (unobserved ones are estimated on the probability floor)" and "a pair that mixes a categorical and a continuous indicator". |

No finding is high severity under the brief's definition for `lpa()`. There,
C2-1 makes EM fail loudly rather than return a wrong answer. Through `lta()`,
however, C2-1 is **silently wrong**: on a separated two-state design the old
rule returned the 1e-8 floor with a log likelihood 23.5 below the maximum
(-1856.93 against -1833.45; `fix/lta_probe_old.R` against `fix/lta_probe.R`),
so it is high there. C2-2 is wrong with a warning that
understates it. C2-3 and C2-4 give silently incomplete or miscalibrated
diagnostics, not estimates.

## Verified correct

| Algorithm | Check (script) | Max discrepancy |
|---|---|---|
| Bounded multinomial M-step (`.multilpa_bound_probabilities`, active-set KKT) | Water-filling solution `p = max(m, c / lambda)` with lambda from a log-scale root; 3,000 random count vectors, K 2-12, bound 1e-10 to 0.9/K, zero counts included (`bound_probabilities.R`) | 7.8e-16 in p; min(p) - bound = 0 |
| Ordinal adjacent-category M-step (`.latents_ordinal_maximize`) | Equivalent Poisson log-linear row-effects model, `glm(n ~ row + col + row:score, quasipoisson)`, 50 random fractional C x K tables (`ordinal_checks.R`) | 6.6e-9 (the glm convergence tolerance) |
| Ordinal likelihood, Fisher-identity scores and observed-information SEs (mixed Gaussian + 2 ordinal, FIML) | Hand mixture LL; numDeriv gradient and Hessian (`ordinal_checks.R`) | LL identical to 13 digits; gradient 6.3e-6 at the fit; SE relative error 5.6e-7; npar 12 = 12 |
| NB count block: mean, dispersion, SEs and LL (one profile) | `MASS::glm.nb` and `theta.ml(eps = 1e-14)`, and numDeriv on (log mu, log alpha) (`nb_glmnb.R`) | mean exact; alpha 4e-11 relative; SE(mean) 1e-9 relative; SE(alpha) 0.05353652 = 0.05353652 against tight `theta.ml`. The 1e-4 gap to `glm.nb`'s default comes from its loose theta tolerance (`eps^0.25`) |
| NB "varying" Poisson-limit rule (score at 0 <= 0 implies the floor) | 4,000 random weighted single-profile cases, grid search of the alpha profile (`nb_varying_weighted_scan.R`) | no case where an interior alpha beat the floor (largest gain 0) |
| Poisson mean on its floor (a profile with all-zero counts) | `poisson_zero_profile.R` | correctly refused: `latents_boundary_fit`, with a message in the count table |
| Direct effects: likelihood, npar, effect and mean SEs (FIML, two targets, shared covariate) | Hand mixture LL with shifted indicators; numDeriv Hessian (`direct_checks.R`) | LL identical; npar 16 = 13 + 3; SE relative error 2.1e-9 (effects), 6.8e-8 (means) |
| Direct-effect centring (`.latents_center_direct` / `.latents_uncenter_direct`) | Refit with z + 2000, z + 1e6, z x 1e-4 and z x 1e4 (`direct_checks.R`) | LL changes by at most 1.5e-10; effect changes by at most 1.6e-7 SE; SE relative change at most 5.3e-8 |
| Direct effects in `predict(newdata)` | Posteriors of a refit on the fitting data (`bvr_gaussian_hand.R`) | 4.4e-16; assignments identical |
| Continuous local-dependence pair (block-diagonal covariance M-step; free index on both scales) | Hand LL with `mvtnorm`, FIML, `variance_model` varying and equal; numDeriv Hessian (`local_dependence_checks.R`) | LL identical to 13 digits; npar 19/14 as expected; off-pattern entries exactly 0; pair-covariance SE identical to 7 digits; means SE relative error 4.9e-8 |
| Categorical joint item (likelihood, npar) | Fit on a hand-joined factor column (`local_dependence_checks.R`) | LL identical to 13 digits; npar 13 = 13 |
| Gaussian BVR (weighted residual correlation, Fisher z, n_eff, direct effects removed) | Hand computation per profile and pair with FIML rows (`bvr_gaussian_hand.R`) | 0 (observed, statistic and effective_n) |
| Gaussian BVR size under H0 | 200 replications (`bvr_calibration.R`) | 4.8% (profile) and 5.8% (overall) at 5%; 1.5%/1.2% at 1% |
| BVR model-implied correlation for structured covariances | VVV, VEE, EEV, VEI, EEE: `expected` against the fitted Sigma_k correlation (`bvr_coverage.R`) | identical to 6 digits |
| Structure-inference charts (all ten: EII, VII, VEI, EVI, VEE, EVE, VVE, EEV, VEV, EVV): volume, det-one shape, Cayley rotation, det-one log-Cholesky, their Jacobians, the chain-rule score and the parameter counts | Chart-invariant SEs (means, every covariance entry, shares) recomputed in an **independent** chart: matrix-exponential rotation, Helmert shape contrasts, `C0^1/2 expm(M) C0^1/2` for the Cholesky forms. Plateau central-difference Hessian, h = 1e-3 checked against 3e-4. Complete data (`structure_se_check_fd.R`) and FIML with 25% missing rows (`structure_se_check_fiml.R`) | relative SE error at most 4.0e-6 (complete) and 2.8e-6 (FIML), the same as the FD noise of the reference (2e-6 to 4e-6); npar equal to mclust's count for all ten; package SEs unchanged by `step` from 1e-3 to 1e-6 |

A note on method. numDeriv's default Richardson Hessian was **not**
accurate enough as a reference for the rotation charts: it moved the SEs by
15% under other settings. The first comparison (`structure_se_check.R`, up to
9e-4 relative) was its error, not the package's. The plateau
central-difference Hessian removes it.

## Not checked (and why)

- Ordinal or NB scores and SEs inside `lta()` and two-level fits: not in this
  partition. C2-1 affects them through the shared M-step, but I did not
  reproduce it there.
- A latent-variable-level (Latent GOLD) run for joint items, NB or BVR: LG
  under Wine is outside a single-core review. The formulas were checked
  against closed forms, glm/glm.nb and numerical derivatives instead.
- Power of the BVR tests under real local dependence: only the size under
  H0 was measured.
- Locale dependence of character category order (round-1 F6): not rerun.
- Structure charts at exactly repeated eigenvalues (EEV/VEV with a spherical
  shape, or VVE/EVE where the anchor combination has ties): the Cayley chart
  is then not locally injective, and I did not test what the inference does
  there. This is a degenerate case.
- The size of the categorical BVR in sparse tables or with many categories
  (the reported calibration uses binary and 3-category items with n = 500).

## Fixes applied (same day, at the coordinator's request)

Scripts for this section are in `round2/fix/`.

| id | change | files |
|---|---|---|
| C2-1 | The MR-F4 rule now also guards `.latents_poisson_limit`: a dispersion moves to the floor only if its own objective does not fall (for a shared dispersion, the objective summed over profiles). | `R/extra-indicators.R` |
| C2-2 | A half-observed pair row is labelled `"a::"` or `":::b"` by `.latents_join_pairs`. It keeps an `NA` code, and its set of joint categories travels as the `"latents_partial"` attribute of the code matrix (`.multilpa_attach_partial`). The density sums the set's probabilities. The M-step (`.multilpa_partial_maximize`) is an inner EM for a concave objective over the bounded simplex: a minorize-maximize step, tolerance 1e-13, at most 10,000 steps, and `latents_no_converge` if that bound is hit. The scores replace the indicator by its expectation within the set. predict(), the rebuilt model frame (bootstrap, sensitivity and other refits) and `.latents_categorical_labels` carry the rows. The `latents_partial_pair` warning is gone; its entry in `?latents-conditions` now says it is no longer raised and why. | `R/categorical.R`, `R/local-dependence.R`, `R/inference.R` (`.multilpa_response_scores`), `R/predict.R` (categorical block), `R/accessors.R` (`.multilpa_model_frame` categorical block), `R/conditions.R` |
| C2-3 | Every pair now has a row. Categorical and ordinal pairs both go through the Pearson table, with an ordinal's levels as its categories. Gaussian–discrete pairs (`kind = "mixed"`) and pairs with a count (`kind = "count"`) get the posterior-weighted within-profile correlation of the values, with `expected = 0` and `statistic`, `df` and `p_value` set to `NA`. | `R/bivariate-residuals.R` |
| C2-4 | The categorical chi-square is documented as conservative, with the measured sizes, in `?get_results` and in `vignettes/articles/workflow-lca.Rmd` (the prose before the `p_adjusted < 0.01` filter). The optional parametric-bootstrap p-value is **left for later**. Each bootstrap draw needs a full refit, so a 100-replication calibration test would take about 10^4 fits, and the argument plumbing sits in `get_results()` code outside the roxygen entries I was allowed to edit. | `R/get-data.R`, the vignette |
| C2-5 | The *Bivariate residuals* section of `?get_results` now gives each kind's statistic and null distribution, with Fisher (1921), Cramer (1946) and Oberski, van Kollenburg & Vermunt (2013) in the references. | `R/get-data.R` |
| C2-6 | The joint-item paragraph now says "the full grid; unobserved combinations on the floor and counted" and "a *mixed* pair is refused", and the `"local_dependence"` table entry says the same. | `R/multilpa.R` (that paragraph only), `R/get-data.R` |

**man/ was not regenerated.** `devtools::document()` writes `man/`, which
was outside my permissions, so `man/get_results.Rd`, `man/multilpa.Rd` and
`man/latents-conditions.Rd` still show the old text.

Tests that asserted the old C2-2 behaviour were updated, with a comment
that the old expectation was the bias itself:
`tests/testthat/test-features-joint-items.R` (a half-observed row now
predicts like the hand-computed marginal) and
`tests/testthat/test-local-dependence.R`.

### Verification

| Check | Result |
|---|---|
| Fits without half-observed rows are unchanged bit for bit. The snapshot covers an LD fit (log-likelihood, posteriors, probabilities, observed and robust inference tables, local-dependence table), a plain FIML LCA (with its residual table) and a mixed Gaussian/categorical/ordinal/NB fit; `fix/baseline.R`, saved before any edit and again after all of them | `identical()` TRUE for all 13 objects |
| C2-1 through `lpa()` and `lta()` | alpha equals the hand profile maximum (1.32278 and 0.747408). The old rule failed every start in `lpa()` and **silently** returned the floor in `lta()`, with a log likelihood 23.5 lower |
| C2-2 bias simulation (n = 40,000, MAR) | P(u1 = 1 \| class) 0.2003 / 0.7962 (truth 0.2 / 0.8), shares 0.6001 / 0.3999. Every joint cell is within 3 SE of the truth (tested) |
| C2-2 likelihood, observed and robust SEs against a hand likelihood with numDeriv (`fix/partial_se_check.R`) | log likelihood equal to 1e-10; observed SEs 4.5e-8 relative; robust SEs 8.7e-8 relative, using the package's CR1 n/(n-1) factor |
| C2-2 `predict(newdata = fitting data)` | posteriors identical to the fit's (difference 0) |
| C2-2 enumeration, refit from the rebuilt data, bootstrap | `enumerate_lca` cell log likelihood equal to `lca`'s; the refit from the model frame reproduces the log likelihood at the fit's seed (1e-9); a bootstrap with `iter = 3` runs |
| New file `tests/testthat/test-review-1026-categorical-r2.R` | 8 tests, 28 passes, 0 failures or warnings |
| Does each new test fail on the old code? | Yes. I copied the package tree, reverted C2-1, C2-2 (encoder, M-step, join) and C2-3 in the copy (`fix/revert.py`) and ran the file there: **all 8 blocks fail** (2+1E, 1, 3, 3, 1, 2, 6 failures). The C2-4/5/6 changes are documentation only; I checked them by reading the roxygen source, and no test covers them |
| Every file matching `categorical\|ordinal\|count\|local\|joint\|direct\|residual\|bvr\|extra\|review-1026` (35 files) | 313 tests, 1,628 passes, 0 warnings. Failures appear in only one file: `test-review-1026-inference-r2.R` (5 failures, 3 errors). Another agent wrote that file at 21:46 while these ran, and it fails the same way (5/3) on the tree with my fixes reverted, so these failures are not mine |
| Harness suites `polca`, `hard-polca` and `glca` (`run-suite.R`) | 2,198 quantities (294,593 values) compared, all agree; 0 outputs not covered |
