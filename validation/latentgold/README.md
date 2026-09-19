# Latent GOLD comparison

Latent GOLD 6.1 is **Windows only** and offers a free academic single-user
licence. It cannot be run on the machine this package is developed on, so this
directory follows the same pattern as `validation/mplus/`: the datasets and the
targets are prepared here, someone runs Latent GOLD on a Windows machine, the
output files are committed, and `compare.R` checks them offline.

## Why Latent GOLD rather than more Mplus

Three features in this package have no external anchor at all, and Latent GOLD
is the reference implementation for all three.

- **Bivariate residuals for continuous indicators.** Mplus's TECH10 is
  categorical only. Latent GOLD's BVR is the standard.
- **Step-3 with the BCH correction.** Jeroen Vermunt wrote the three-step
  papers and co-authored Latent GOLD, so its Step-3 is as close to the seminal
  source as software gets. The current check is against `tidySEM`, a
  third-party reimplementation.
- **Step-3 with covariates**, the R3STEP equivalent.

## What to run

Two datasets are written here by `make-fixtures.R`, tab separated with a header.

### 1. `bvr.dat` — bivariate residuals

400 rows, 40 groups, columns `id`, `group`, `y1`, `y2`, `y3`.
`y1` and `y2` share a term the profiles do not explain; `y3` does not.

Fit: **two-level cluster model**, 2 latent classes at the case level and 2 at
the group level, three continuous dependents, local independence, and request
bivariate residuals in the output.

This package's values, for the same data:

| quantity | value |
|---|---|
| log likelihood | -1990.324253 |
| free parameters | 15 |
| BVR, `y1` with `y2`, profile 1 | 0.5841 |
| BVR, `y1` with `y2`, profile 2 | 0.5664 |
| BVR, `y2` with `y3`, profile 1 | 0.0543 |

Latent GOLD reports BVR as a chi-square-like quantity per pair rather than as a
residual correlation, so the two will not be numerically identical. What must
agree is the **ranking and the separation**: `y1` with `y2` far above every
other pair, and the other pairs near zero. If Latent GOLD flags a pair this
package does not, or the reverse, that is the finding.

### 2. `threestep.dat` — BCH and covariates

600 rows, 50 groups, columns `id`, `group`, `x`, `y1`, `y2`, `distal`.
Class membership was generated as `plogis(-0.3 + 1.2 * x)`; the distal outcome
is 10 in class 2 and 0 in class 1.

Fit: **single-level, 2 latent classes**, dependents `y1` and `y2`, then Step-3
twice — once with `distal` as the dependent under the BCH correction, once with
`x` as a covariate predicting class membership.

This package's values:

| quantity | value |
|---|---|
| log likelihood | -2024.808558 |
| classification error, true 1 assigned 2 | 0.046048 |
| classification error, true 2 assigned 1 | 0.026016 |
| distal mean, class 1, BCH | 10.0755 |
| distal mean, class 2, BCH | 0.2927 |
| distal mean, class 1, modal | 9.70111 |
| covariate slope on `x` | 1.21877 (SE 0.13546) |

Class labels are arbitrary and will differ; align them by the means of `y1`
before comparing. The generating slope is 1.2 and the generating distal
difference is 10, so both packages can also be scored against the truth.

## Files to commit back

Whatever Latent GOLD writes: the `.lgs` syntax, the `.lst` or `.txt` output,
and any saved parameter files. Then run:

```
Rscript validation/latentgold/compare.R
```

## Syntax

A sketch is in `syntax-sketch.lgs`. It has **not been run** and the exact
keywords should be checked against the Latent GOLD 6.1 manual before use; it is
a statement of the models to fit, not tested code.
