# Astra core-model review contribution

Date: 2026-10-27. Scope: committed estimator changes from release commit `24773e6` (0.97.1) to `5596a15` (0.97.4), inspected against the current working tree. No package source or package tests were changed. Independent scripts and artifacts are under `/tmp/astra-*`. The two new regressions below occur in committed code; the overlapping uncommitted edits do not repair them.

## Confirmed findings

### CORE-1 — P1: Excluded persons change LTA covariate conditioning and the fitted answer

**Locations:** `R/lta-engine.R:125-139`, especially lines 128-129 and 135-138; caller at 981-990. Also apply the eventual correction to the inference chart caller in `R/lta-inference.R`.

**Introduction:** `.lta_design_chart()` was introduced in `d5d4ae8` (0.97.2). The preceding `.lta_designs()` rank check correctly receives `.lta_used_rows(...)`; the chart does not receive this mask. It computes means/SDs over every person's initial design and every `layout$within` transition row. Zero-weight persons and moves beyond a person's last informative observation consequently influence the optimization coordinates.

**Reproduction:** 35 persons × three observations, two separated Gaussian profiles, one initial and one transition covariate, one random start, seed 1. Original active data give LL **−226.9712124143479**, `converged=TRUE`, `boundary=FALSE`. Append five persons with **weight zero**, unchanged ordinary indicator values, and `z=w=1e8`. LL falls by **1.4621** to approximately **−228.4333**, yet `converged=TRUE`. The transition coefficient with the largest change moves **1.577333 standard errors**, measured using the original fit's robust covariance on the same 14 estimation coordinates. Thus this is materially different estimation, not harmless stopping noise or a units-only parameter re-expression.

An in-memory diagnostic replacement that computes the chart from the original 35 carrying persons, while retaining all observations throughout the actual fit, restores the original LL **exactly** with `converged=TRUE`. This isolates the cause from random initialization and ordinary mixture-mode variability. No on-disk implementation was changed.

**Repair direction:** use the already available active initial/move masks for chart estimation, both fitting and inference; preserve the complete designs for prediction of excluded persons. Add deletion-versus-zero-weight-padding regressions with extreme covariates and trailing noninformative occasions. Compare likelihood, coefficients in SE units, convergence and inference.

**Evidence:** `/tmp/astra-core-more.R`, `/tmp/astra-core-more.log`; `/tmp/astra-lta-chart-proof.R`, `.log`; `/tmp/astra-lta-magnitude.R`, `.log`; saved fits `/tmp/astra-lta-1e+08.rds` and `/tmp/astra-lta-chart-corrected.rds`.

### CORE-2 — P2: New constant-response refusal rejects valid regression models

**Locations:** `R/mixture-regression.R:628-638` and weighted counterpart `:507-515`.

**Introduction:** unweighted guard added in `d5d4ae8`; weighted counterpart in `285eca6`.

The statement that a constant Gaussian outcome necessarily produces zero residual variance is false for a model without an intercept and for a varying fixed offset. Both specifications are supported by the existing API and had valid solutions in the release.

Use `set.seed(1027)` and `d <- data.frame(y=rep(2,60), x=seq(-2,3,length.out=60), o=rnorm(60))`.

| Model | Independent ML residual variance | Independent log likelihood | Released 24773e6 | Current tree |
|---|---:|---:|---|---|
| `y ~ 0 + x` | 3.5840188 | −123.4308542 | Same estimates/LL, converged | `latents_bad_data` |
| `y ~ 1 + offset(o)` | .9624063 | −83.9867552 | Same estimates/LL, converged | `latents_bad_data` |

Each fit uses one class, `n_starts=0`, `vcov_type="none"`. References are `lm()` coefficients, `mean(residuals(lm_fit)^2)` and `logLik(lm_fit)`. The released source was extracted with read-only `git archive 24773e6 R DESCRIPTION NAMESPACE man data` to `/tmp/astra-core-release` and loaded in a separate R process. This is a verified regression, not an inference from git alone.

**Repair direction:** detect a truly degenerate residual problem conditional on the actual design and offset; do not categorically reject constant raw responses. Apply any fix consistently to weighted and unweighted paths and ensure a positive residual floor remains well defined.

**Evidence:** `/tmp/astra-core-probes.R`, `.log`; `/tmp/astra-constant-release.R`, `.log`.

### CORE-3 — P2: Zero-weight persons independently disable otherwise valid LTA inference

**Locations:** `R/lta-inference.R:662-685`; call in `R/lta-engine.R:1407-1408`.

**Status:** existing behavior retained through the review commits, and a remaining gap in the newly documented rule that zero-weight units are absent from refusals and warnings. The boundary checker receives neither sampling weights nor the active-observation masks.

Use the CORE-1 fixture but set the five excluded persons' `z=w=10000`. This moderate version preserves LL to **1.99e-13**, and both fits converge. Nevertheless the original fit has `boundary=FALSE` and returns a **14 × 14** covariance, while the padded fit has `boundary=TRUE` and `vcov()` raises `latents_boundary_fit`. Probabilities close to zero only for the excluded covariate values trigger the blanket refusal. The CORE-1 in-memory chart repair likewise leaves the padded fit's boundary flag true, confirming these are independent code paths.

**Repair direction:** evaluate initial and transition boundary checks only where the weighted observed likelihood uses those probabilities; keep structural stayer handling. Assert covariance availability is invariant to zero-weight padding.

**Evidence:** `/tmp/astra-lta-10000.rds`; `/tmp/astra-core-final-checks.R`, `.log`.

### CORE-4 — P2: Categories observed only on zero-weight rows still increase model dimension and BIC

**Locations:** `R/categorical.R:33-42` and `:47-63`. The new `contributing` mask checks the minimum number of active categories but does not determine the fitted category set.

**Status:** category discovery from all rows predates 0.97.1; it remains inconsistent with the new positive-support convention. Distinguish this remaining semantic gap from a newly introduced regression.

One-class LCA on 60 rows, binary item `a` and ternary item `b`, all weights one:

```r
d <- data.frame(a=rep(c("a","b"),30), b=rep(c("x","y","z"),20), w=1)
p <- rbind(d, data.frame(a="absent-category", b="x", w=0))
f <- function(z) lca(z, c("a","b"), 1, n_starts=1, seed=1, weights="w")
```

| Quantity | `f(d)` | `f(p)` |
|---|---:|---:|
| Parameter count | 3 | 4 |
| LL | −107.505568153683285 | −107.505568159683364 |
| BIC | 227.294169994032870 | 231.388514568255118 |

A truly absent row changes BIC by **4.094344574222248** (approximately `log(60)`), despite unchanged active responses. The negligible LL change is the additional probability-floor charge; the meaningful issue is the artificial extra category and parameter. This can affect model comparisons. The ordinary-item discovery policy already drops unused declared factor levels, so this is not a consistent policy of retaining a prespecified full category space. Declared joined-item grids need separate treatment.

**Repair direction:** resolve the ordinary fitted categories on active support and separately specify how excluded rows with unseen categories receive predictions. Cover ordinary factors/characters/numeric categories and paired items explicitly.

**Evidence:** `/tmp/astra-core-more.R`, `.log`; `/tmp/astra-core-final-checks.R`, `.log`.

## Lower-priority confirmed limitation, not a practical-case release blocker by itself

Growth starting covariances and trajectory partitions still read absent persons: `R/growth-mixture.R:132-149`, `:217-236`. On 30 carrying persons × four observations, one random-intercept class, appending four zero-weight persons with outcomes ±1e100 makes fitting raise `latents_em_decrease` (reported decrease 1.32e4), while the carrying data alone converge at LL approximately −115.2865. This closes the earlier review's unclaimed absent-start concern only at an extreme numeric scale. It does **not** establish a typical-data or multi-class local-maximum defect.

Controls matter: the same one-class padding at ±1e5, ±1e8, ±1e10 and ±1e20 preserved LL within 1.9e-13. Eight additional two-class fixtures (40 carrying persons, four observations, 40 appended zero-weight persons shifted ±20, seeds 1–8) all converged to their deletion-reference LL within 1.2e-13. Do not overstate the extreme failure as a general search failure.

Evidence: `/tmp/astra-core-probes.R`, `.log`; `/tmp/astra-growth-zeros.R`, `.log`, and `/tmp/astra-growth-seed-1.rds` through `-8.rds`.

## Coverage and assessment of additions

This pass inspected the committed diffs and relevant implementations for the shared EM stopping/counting changes; stable posterior normalization; weight normalization/support; joined categorical partial observations; additive and cross-level likelihood/count corrections; Gaussian regression refusals; growth standardization, coefficient/covariance maps, starts and Brownian bounds; and LTA covariate charts, active designs and inference boundaries. The previous consolidated review and round-three report were read, including their limitations. There are 80 changed R files in the overall committed range: this contribution is a core-focused audit, not a claim that all 80 were line-audited here.

The main design changes are coherent and valuable: stopping by observed values avoids outcome-unit dependence; raw-positive weight counts preserve the meaning of extreme positive weights; ordinary labels require explicit pair provenance; half-observed joined items contribute a marginal likelihood; growth fitting uses standardized coordinates and maps back; and LTA aliases are checked where observations actually inform them. The confirmed defects are failures to propagate these rules to every consumer, plus an overbroad input guard, rather than evidence that all new estimators are unsound.

No new full Monte Carlo calibration or external-program equivalence campaign was run by this subreview. No independent maximal-fit conclusion was reached for the earlier AR multilevel-growth probe that remained nonconverged in both unit systems. The parent is running the full package suite. This contribution's evidence consists of small deterministic synthetic probes and direct baseline/current comparisons, not an independent full-suite run.

All current probes used R 4.5.2 with `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`, generally `NOT_CRAN=true`, and either `pkgload::load_all(".", quiet=TRUE)` or the parent's fresh installed snapshot at `/private/tmp/astra-review/lib`. The magnitude and final categorical/boundary reproductions were rerun against that installed snapshot. No network access or extra package installation was needed.
