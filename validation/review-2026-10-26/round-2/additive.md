# Round-2 review: additive group-class families (slug: additive)

Date: 2026-10-26. Every check was run in this session against the working tree, with
`pkgload::load_all()` and `NOT_CRAN=true OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`.
The reproduction scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-additive/`.
In the tables below, `S/` stands for that folder. `S/setup.R` holds the shared data generators.
No external implementation of this model exists. The independent checks are therefore closed forms,
Gaussian conditioning, dense multivariate-normal densities, nlme for one class, and numerical derivatives
(numDeriv, plus Richardson extrapolation of the analytic score).

## Scope (files)

- R/additive.R covers the sufficient statistics, exact density, E-step, M-step, KKT boundary pass, starts,
  fit object and criteria.
- R/additive-inference.R covers packing, scores, information, observed/robust/OPG covariance, the
  natural-scale Jacobian, Newton polish and effective groups.
- R/additive-methods.R covers the tables, `nobs`, `parameter_inference(data =)` and the identity check.
- R/additive-compare.R covers nesting, simulation and `bootstrap_lrt()`.
- The intended model is described in validation/ADDITIVE_MODEL_DESIGN.md, additive-design-check.R and
  ADDITIVE_SIMULATION.md.
- The round-1 findings A1 (zero scatter flag) and A2 (identity tolerance) are fixed in the tree. I did
  not re-report them; A2's fix was re-verified (see V17).

## Findings

| id | severity | file:line | what the code does | what it should do | reproduction | measured effect | proposed fix |
|---|---|---|---|---|---|---|---|
| ADD2-1 | medium | R/additive-inference.R:491-492 (`information_ratio`), R/additive-inference.R:280 (`latents_weak_class` text) | The information ratio is `w(1-w) / (x$n_groups * Var(w))`. The warning's "expected groups" is `x$n_groups * w`. `x$n_groups` counts every group, including groups of sampling weight zero. | Since today's change the weights sum to the groups of positive weight (`weighted_n_groups` = `n_used_groups`). The complete-data variance of `w` uses that count. Zero-weight groups must change nothing (ML-1 contract, `?multilpa` `drop_empty` text: "Units of sampling weight zero are excluded from every count and entropy term"). | `S/t3_zero.R` | Adding 40 zero-weight groups to 80 weighted groups leaves LL, estimates, SEs and `effective_groups` unchanged (V13). However, `information_ratio` drops from 0.6786 to 0.4524 (additive), 0.7057 to 0.4704 (additive equal), 0.6506 to 0.4337 (dispersion) and 0.7631 to 0.5088 (additive-dispersion). The warning then reads "28.8 effective of 63.7 expected groups, 45% of the information kept" instead of "of 42.5 expected groups, 68%". That blames poor separation for what is only padding. | Use `x$sufficient_statistics$weighted_n_groups %||% x$n_groups` (equivalently `x$n_used_groups`) in both places. |
| ADD2-2 | medium (entropy); low (rest) | R/additive-methods.R:184 (`entropy`), :121-124 and :167-178 (`n_assigned`, `mean_posterior`, `classification`), :183-184 (`n_groups`, `n_obs` columns); R/additive.R:812 (`small_classes`, which also gates inference at R/additive-inference.R:201); R/additive.R:886 (print) | These statistics use every group, unweighted, including zero-weight groups. | The package contract says zero-weight units are in no count and no entropy term. The profile family's relative entropy is unweighted over units of positive weight (`?get_results`, criteria section). The group-classes table already weights `count`, so it mixes conventions inside one table. | `S/t3_zero.R` | Relative entropy goes from 0.9057 to 0.9334 (additive), 0.9116 to 0.9281, 0.7943 to 0.8449 (dispersion) and 0.9578 to 0.9719. `n_assigned` goes from 46/34 to 46/74, `mean_posterior` from 0.9721 to 0.9862, and the fit table reports `n_groups` 120 while `nobs()` is 80. `entropy` is a column of the `enumerate_classes(family = )` grid, so model comparison sees it. `small_classes` can pass, or refuse inference, because of zero-weight groups. | Restrict these statistics to `sampling_weights > 0` rows (`.additive_entropy(posterior[used, ])`, `tabulate(assigned[used])`, and so on). Set `small_classes <- any(colSums(posterior * (w %||% 1)) < 1)`. Report `n_used_groups` and `n_used_obs` in the fit table and in print, or add them next to `n_groups` and `n_obs`. |
| ADD2-3 | medium (b); low (a, c, d) | R/additive.R:67-79 (constant-indicator and all-singleton refusals run before the weights are read at :101); :88-90 (`zero_within_scatter`); :728 (`n_group_classes > stats$n_groups`) | The identification refusals look at all groups, including zero-weight groups. | Identification is a property of the units that carry weight. Without the zero-weight groups the same fit is refused. | `S/t8_refusals.R` | **(a)** An indicator constant in every positive-weight group but varying in zero-weight groups is accepted. It only gets `latents_boundary`, while the reference refuses with `latents_bad_data`. **(b)** Every positive-weight group is a singleton and only zero-weight groups are replicated. The fit returns `converged = TRUE` with W = (0.825, 0.938) and T = (1.599, 0.514), with no warning. That split of W+T is arbitrary, because the likelihood has a ridge. Only inference refuses (eigen ratio 5.7e-15). The reference refuses with `latents_unidentified`. **(c)** 4 group classes from 3 positive-weight groups are fitted, with weights 0.666/0.333/4e-4/3e-4 and BIC 104.94, instead of `latents_unidentified`. **(d)** y2 is constant within every positive-weight group: `zero_within_scatter` = FALSE, FALSE (the `within_floor` boundary is still flagged). | In `.additive_prepare()`, compute the weights first. Then run the constant, singleton and zero-scatter tests on rows whose group weight is positive. In `.additive_fit()`, compare `n_group_classes` with `n_used_groups`. |
| ADD2-4 | low | R/additive.R:344-409 (`.additive_initialize`) | Starting values ignore sampling weights and use every group. This covers the Ward features of start 1, random centres drawn from all groups, initial weights `counts / n_groups`, pooled within `pooled_scatter / (N - J)` and the between floor. | Zero-weight groups should not move the starts either. Otherwise the reported local maximum can depend on them. | `S/t13_starts.R`, `S/t13b_starts.R` | 3-class additive, 60 groups plus 30 far-away zero-weight groups, `n_starts = 1`: the LL differs in 8 of 12 datasets, by up to 2.49 (−959.656 vs −962.148). With the default starts all four affected datasets reach the same maximum (|ΔLL| ≤ 1.2e-5). | Run the initializer on the positive-weight groups only, with weighted moments, and give zero-weight groups any assignment. |
| ADD2-5 | low | R/additive.R:812-874 (no fit-level check) | A start can converge to two identical classes: dispersion fit, class within variances equal to 3e-5, logit direction flat (LL change 1e-9 over ±0.2). It is reported as `converged` with no warning. | The design (ADDITIVE_MODEL_DESIGN.md, "Identification") asks for classed diagnostics for duplicate components. Inference does refuse correctly (`latents_singular_information`, eigen ratio −0.0523). | `S/t7_saddle.R` | Weighted fit: within (1.12781, 0.98757) vs (1.12784, 0.98768), Hessian eigenvalue 2.5e-8 against a largest of 66. The unweighted fit is exactly collapsed (eigenvalue 2.1e-9). | Warn (for example `latents_duplicate_classes`) when two classes' parameters agree to a tolerance on the SE scale, or when the information has a null direction at the fit. |
| ADD2-6 | info | R/additive.R:764 | Newton polish is applied only to interior fits. A boundary fit (one T held at 0) is left where EM stopped. | Not a defect: EM stops on the declared `tol`, as every reference program does. | `S/t12_oneclass_boundary.R` | One class, balanced, y2 T at the boundary: the fit is 1.1e-6 below the closed-form maximum. The differences are T1 +0.0015 SE and W1 −0.0006 SE. | Optional: also run Newton steps on the free coordinates of a KKT boundary fit. |
| ADD2-7 | low | R/additive-compare.R:115 | The additive `bootstrap_lrt()` accepts `iter >= 1`. | `?bootstrap_lrt` says `iter` must be at least 2, and the profile path enforces `iter >= 2L`. | `S/t9_blrt.R` | `bootstrap_lrt(n1, f, iter = 1)` returns a result instead of refusing. | Use `iter >= 2`. |
| ADD2-8 | low (docs) | R/additive.R:25-26 and :34-43; R/multilpa.R:1205-1207 | (1) The `.additive_prepare()` roxygen says the sampling weights are "scaled to sum to J". They now sum to the positive-weight groups, and the `n_used_*` and `weighted_n_groups` fields are undocumented. (2) The `.additive_identity_tolerance` block is spliced between `.additive_prepare`'s roxygen and its function, so the prepare docs attach to the constant. (3) `?multilpa` says these families "take only data, vars, id, ..., seed and weights", but `search` and `n_screen` are accepted (`.additive_arguments()`). | Make the docs match the code. | (reading) | None numerical. | Edit the text, and move the constant above the prepare block. |
| ADD2-9 | info | R/additive-methods.R:315-318, R/additive-compare.R:78-81 | `all.equal(tolerance = 64 eps)` is a mean relative difference over the elements that differ, not an elementwise bound. A genuine change in one mean is diluted by ulp-level differences in the others. | Acceptable as it stands. The diluted change has to be below about `64 eps * J * d` relative, far below any printed digit. | `S/t6_ulps.R` | See V17: permutation reaches at most 18 eps of the 64 allowed, and a rating changed by 1e-6 is refused. | None needed. An elementwise `abs(a-b) <= 64 eps * (|a| + scale)` would be stricter. |

## Verified correct

| id | algorithm / formula | independent check | max discrepancy | script |
|---|---|---|---|---|
| V1 | Component log density `ell_jh` (raw-rating, intercepts integrated out), class-specific W and T, sizes 2-8 | `mvtnorm::dmvnorm` of each group's ratings under `W I + T 11'` | 7.1e-15 | S/t1.R |
| V2 | Weighted pseudo log likelihood `sum_j w_j logsumexp_h` | Hand computation from V1 | 1.1e-13 | S/t1.R |
| V3 | E-step conditional intercept mean `mu + nT/(W+nT)(a-mu)` and variance `TW/(W+nT)` | Gaussian conditioning of (B, y) with the dense covariance | 4.4e-16 (mean), 5.6e-17 (variance) | S/t1.R |
| V4 | Marginal intercept moments (mixture over the class posterior) | Hand formula from V3 | 4.4e-16, 8.3e-17 | S/t1.R |
| V5 | M-step, weighted, every structure (additive varying/equal, dispersion, additive-dispersion varying/equal), with and without held-at-zero coordinates (one cell, a whole column): means, between (per class, or pooled over `weighted_n_groups`), within (per class, or pooled over `weighted_n_obs`), weights `counts / weighted_n_groups` | The gradient of the weighted expected complete-data log likelihood Q with respect to each block is zero at the update, in ECM order (means at the old W where held coordinates are pooled) | Gradient ≤ 1.1e-7 (numDeriv noise on Q of about 1e3); weights exact (0); every step raised LL (+7.2 to +30.2) | S/t2_mstep.R |
| V6 | Analytic per-group scores (Fisher identity), weighted, away from the optimum | `numDeriv::jacobian` of the per-group log likelihood times w | ≤ 2.1e-8 relative | S/t4_scores.R |
| V7 | Observed information (5-point difference of the analytic score in a centred, unit-scaled chart) | Richardson Jacobian of the analytic score; `numDeriv::hessian` of the LL | 9.7e-10 against Richardson. The numDeriv Hessian itself is off by 7.7e-6, its own error. | S/t4b.R, S/t4_scores.R |
| V8 | Robust covariance: `bread (G/(G-1)) sum(w_j s_j)(w_j s_j)' bread`, with G the positive-weight groups (CR1) | Hand sandwich with the Richardson bread | SE ratio 1 to 7 digits (≤ 1.6e-5 with the noisier numDeriv bread) | S/t4b.R, S/t4_scores.R |
| V9 | Observed and OPG covariances (unweighted) | `solve(-numDeriv::hessian)`, `solve(crossprod(numDeriv scores))` | ≤ 2.7e-7, ≤ 1.7e-9 relative | S/t4_scores.R |
| V10 | Natural-parameter Jacobian (identity, `diag(variances)`, softmax with class 1 dropped) | `numDeriv::jacobian` of `.additive_natural()$estimate` | ≤ 2.4e-9 | S/t4_scores.R |
| V11 | Parameter-table SEs, `sqrt(diag(J V J'))` | Hand delta method with the hand sandwich | ≤ 8.0e-6 relative (numDeriv bread) | S/t4_scores.R |
| V12 | One class, balanced, interior closed forms `mu = Σwa/Σw`, `W = ΣwS/(Σw(n-1))`, `T = Σw(a-mu)²/Σw - W/n`, unweighted and weighted with two zero weights, all 3 families. Also BIC = −2LL + p log(n_used_groups) and `nobs` | Closed form | mean 4.4e-16, W 4.4e-16, T 3.3e-16, LL 1.1e-13, BIC 0, `nobs` 38 = positive groups | S/t5_oneclass.R |
| V13 | One class, balanced, boundary closed form (T<0 gives T=0, `W = [ΣwS + nΣw(a-mu)²]/(nΣw)`) | Closed form | Weighted: LL 1.4e-11, W 1.9e-7. Unweighted: see ADD2-6 | S/t12_oneclass_boundary.R |
| V14 | One class, unbalanced | `nlme::lme(method = "ML")` per indicator | The package LL is 1.9e-9 higher than nlme's; T 7.0e-6, W 2.6e-6, mean 5.8e-8; intercept SE ratio 1.00003 and 1.0018 (observed vs nlme's expected-information GLS SE) | S/t5_oneclass.R |
| V15 | Zero-weight invariance of the estimator and criteria: 40 zero-weight groups (sizes 1 and 12) added to 80, 4 family configurations | Refit and compare | LL ≤ 4.5e-13, estimates ≤ 1.0e-10, robust SEs ≤ 3.5e-11 relative, AIC/BIC/`bic_individual` ≤ 9.1e-13, `nobs` equal, `effective_groups` ≤ 2.2e-10 (other outputs: ADD2-1, ADD2-2) | S/t3_zero.R |
| V16 | KKT boundary pass, weighted | (i) Analytic `dLL/dT` at 0 against a one-sided numerical derivative. (ii) L-BFGS-B with T ≥ 0 started from the fit | (i) −22.9 vs −22.9 and −10.8 vs −10.8; (ii) gain ≤ 6.9e-9 | S/t10_kkt.R |
| V17 | Identity check at 64 ulps (`parameter_inference(data =)`, `bootstrap_lrt`) | Permute rows within groups, for J = 10-60 and n = 5-20000, with centred and offset (1e3) data, 3 replications | Permutation difference ≤ 18 eps (accepted). End to end at J = 30, n = 3000: permuted rows accepted, one rating changed by 1e-6 refused (`latents_bad_inference_data`) | S/t6_ulps.R |
| V18 | Newton polish reaches the maximum (weighted, 3 families) | BFGS with the analytic gradient, `reltol = 1e-16`, from the fit | Gain ≤ 4.6e-13 | S/t11_newton_scale.R |
| V19 | Weighted scale and origin equivariance, y → c + 1000y | Refit | LL shift −5207.457555, which equals −2 Σ w_j n_j log 1000 exactly; estimates ≤ 2.7e-10 SE; SE ratio ≤ 1.0e-10 | S/t11_newton_scale.R |
| V20 | `bootstrap_lrt()` for this family: simulation, statistic, p-value | Moments of 400 simulated data sets; hand p-value | Mean of S/(n-1) 1.3347/1.2546 vs model 1.3319/1.2451; mean of a −0.3895/0.3407 vs −0.3816/0.3372 (Monte Carlo SE 0.0115/0.0082); statistic = 2ΔLL exactly; p = (1 + #{T* ≥ T})/(iter + 1) = 0.1, as computed by hand | S/t9_blrt.R |
| V21 | Parameter counts (means + between + within blocks, plus H−1) equal `length(theta)` | Counted | 11, 9, 9, 13, 11 for the 5 structures | S/t4_scores.R |
| V22 | Nesting rule (`min_variance`, class count, block freedom) and weighted refusal in `bootstrap_lrt()` | Reading; weighted fits are refused upstream (R/enumeration.R:1325) | — | — |

## Not checked (and why)

- I did not rerun Wald coverage or the `latents_weak_class` threshold of 50. That needs the simulation study
  (ADDITIVE_SIMULATION.md), which is too long for this pass.
- I made no external (Mplus or Latent GOLD) comparison. Mplus is absent, and no Latent GOLD syntax for this
  model exists in the repository.
- `search = "emEM"` and `cores > 1` were not exercised for this family. The only start paths run were the
  default kmeans starts and `n_starts = 1`.
- The `recovery` and `conditional_intercepts` tables and the plots were not checked. They have no formula
  beyond V3 and V4.
- Cross-level and structure-inference code is outside this round's file list (round 1 verified both).
- The BLRT's size or power was not checked: 9 replicates only check the mechanics, not the calibration.

## Fixes applied (same day, at the coordinator's request)

Coordinator decision: groups of sampling weight zero are absent from every count, entropy, classification
table, starting value, refusal and warning. They still receive a posterior and a per-group row.

| id | status | change |
|---|---|---|
| ADD2-1 | fixed | R/additive-inference.R: new `.additive_sample_groups(x)` (`n_used_groups %||% n_groups`). It is used in `information_ratio` and in the "expected groups" of `latents_weak_class`. |
| ADD2-2 | fixed | R/additive-methods.R: `entropy`, `n_assigned`, `mean_posterior`, the `classification` table and `recovery` now use positive-weight groups only. The fit table reports `n_groups` = `n_used_groups` and `n_obs` = `n_used_obs`. R/additive.R: `small_classes` counts positive-weight groups, and `print()` shows the used counts. |
| ADD2-3 | fixed | R/additive.R `.additive_prepare()`: weights are read before the checks. The constant-indicator, all-singleton and `zero_within_scatter` tests run on positive-weight groups. `.additive_fit()` compares `n_group_classes` with `n_used_groups`. Unweighted fits are unchanged. |
| ADD2-4 | fixed | R/additive.R `.additive_initialize()`: starts are built from the positive-weight groups alone (unweighted moments, as before). Fits without zero weights get identical starts and identical random draws. |
| ADD2-5 | fixed (narrowly) | R/additive.R: new `.additive_collapsed_classes()`. Two classes collapse when every class-specific parameter agrees to 1e-8 relative (variances relative to the pair, means relative to `sqrt(W + T)`). The result is a new flag `collapsed_classes`, a `latents_weak_class` warning at fit time, and `latents_boundary_fit` from `.additive_inference()`. Limitation: the weighted example in `S/t7_saddle.R` collapsed only to 1.1e-4 relative (Newton cannot run on its singular information), so the 1e-8 rule does not flag it. Its inference is still refused (`latents_singular_information`). The unweighted example (9.5e-11) is flagged. |
| ADD2-6 | documented | Boundary fits are not Newton-polished. The measured shortfall (1.1e-6 log likelihood, at most 0.0015 SE) is stated at the Newton call in `.additive_fit()` and pinned by a test (bound 0.002 SE). |
| ADD2-7 | fixed | R/additive-compare.R: `iter >= 2`. |
| ADD2-8 | fixed (two of three) | The `.additive_prepare()` roxygen now documents the weight scaling and `weighted_n_groups`/`n_used_*`, and the tolerance constant has its own block. The `?multilpa` wording (R/multilpa.R, not edited here) is passed to the coordinator. |

Tests: `tests/testthat/test-review-1026-additive-r2.R` has 7 tests and 31 expectations, all passing.
Against the old code (a copy of the package with the four original additive files restored), the tests
for ADD2-1, -2, -3, -4 and -7 fail (2, 7, 4, 1 and 1 failures) and the ADD2-5 test errors. The ADD2-6
test passes on both, because that item is documentation.

The 32 test files matching `additive|weight|review-1026` give 268 tests, 1,300 expectations, 0 skipped.
The only failures (5) and errors (3) are in `test-review-1026-inference-r2.R`, which contains no additive
code. It fails identically with the old additive files restored, so it belongs to another agent's
in-progress work. `Rscript validation/additive-design-check.R` passes all 9 checks.
