# Class enumeration and the bootstrap LRT against mclust and poLCA

`enumerate_lpa()` is compared with `mclust::mclustBIC()`, `enumerate_lca()`
with poLCA, and `bootstrap_lrt()` with `mclust::mclustBootstrapLRT()` and the
statistics printed in mclust's vignette.

```
Rscript equivalence/enumeration-packages/compare.R   # comparison.csv, search.csv, bootstrap.csv, RUN.txt
```

`equivalence/test-enumeration-packages.R` runs the same comparison under
testthat.

**Result (2026-10-29): 168 quantities compared, 164 agree, 4 known
differences (below), 0 unexplained.**

| function | reference | grid |
|---|---|---|
| `enumerate_lpa()` | `mclustBIC()`, EM to 1e-12 | EEI, VVI, EEE, VVV at 1-4 profiles; `faithful`, `iris`, mclust's `diabetes` |
| `enumerate_lca()` | poLCA, 30 starts, tol 1e-12 | 1-4 classes; poLCA's `carcinoma`, `values`, `cheating` (4 items) |
| `bootstrap_lrt()` | `mclustBootstrapLRT()`, 199 replicates, EM to 1e-12; vignette | diabetes VVV, 1 vs 2, 2 vs 3, 3 vs 4 |

latents fits every cell with 50 starts, seed 1, tol 1e-12. Every cell's
parameter count is compared. The log likelihood and BIC are compared where
both programs reach the same maximum (`search.csv` lists each cell's maxima).

Tolerances, fixed before the comparison: the log likelihood to 1e-6 plus
1e-10 of its size, the BIC (−2 LL plus a constant) to twice that; parameter
counts exactly; a bootstrap p-value within four Monte Carlo standard errors
of the two estimates combined, plus 1 / B; printed statistics to half a unit
of their last digit.

## mclust's tolerance

mclust's EM stops by default at a relative change of 1e-5 in the log
likelihood. Here that stopped up to 0.04 short of the maximum (faithful EEI,
four profiles), so 31 of the grid's 60 cells looked like different maxima;
converged, 12 are.
Every mclust fit is therefore run with `control = emControl(tol = c(1e-12,
1e-12))`, which `mclustBootstrapLRT()` also passes to its replicate fits.

## Search (`search.csv`, not pass/fail)

36 of the 48 Gaussian cells and all 12 latent class cells reach the same
maximum. In the other 12 latents' maximum is higher, by 0.32 to 108 (iris VVV
with four profiles). mclust is never the higher. For diabetes VVV, latents
reaches −2343.012 with two profiles and −2281.257 with four, against mclust's
−2355.638 and −2285.018. These are proper solutions: mclust's own EM started
from latents' posteriors stays there, the smallest class has 12 rows, and
the smallest covariance eigenvalue is 57.

## Bootstrap LRT (`bootstrap.csv`)

mclust's hierarchical start reaches lower maxima at two and four profiles,
so the bootstrap is compared at mclust's own solutions: `lpa()` is started
there (`start =`, one start) and each replicate is refitted with its 50-start
search. The observed statistics then agree with mclust's to 6e-8 (370.05113,
121.09004, 20.15106).

The data are the vignette's. mclust 6.1.2's vignette renames the columns
(glucose, insulin, sspg become fpg, glucose, insulin) and sets
`glucose[104] <- 455`, which is the installed data's `insulin[104]`. The grid
uses the same corrected data.

## Known differences (`known-differences.csv`)

- **The vignette's three statistics** (369.98537, 121.10652, 19.90071) are
  mclust's at its default tolerance. mclust at that tolerance reproduces them
  exactly; converged, it gives the values latents agrees with.
- **The 3 vs 4 p-value**: latents 0.965, mclust 0.425. `mclustBootstrapLRT()`
  refits each replicate's four-profile alternative by one EM from the
  observed data's four-profile estimates. `bootstrap_lrt()` refits with the
  50-start search, which reaches higher alternative maxima. The mean bootstrap
  statistic is 38.5 for latents against 19.8 for mclust. Over the 104 of 199
  replicates not on the variance floor it is 31.8 (median 28.4), with p 0.933.
  At 1 vs 2 and 2 vs 3 both p-values are 1/200.

## Warnings

Under testthat the kit passes with 480 warnings (2026-10-29). All are latents'
own diagnostics on deliberately over-fitted models, mostly the bootstrap
replicates and the four-class cells: low Gaussian membership mass (292), a
covariance on the variance floor (159), effective class size below one (24),
failed starts (3), and an unidentified model (2: poLCA's `values` at four
classes has 15 cell degrees of freedom for 19 parameters; poLCA fits it too, and
only its maximized likelihood, which is still determined, is compared). They are
left visible, not muffled.
