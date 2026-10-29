# mixture_regression() against independent implementations and published results

Every case is fitted by `mixture_regression()` (50 starts, `tol = 1e-14`) and by
independent maximum-likelihood implementations, each from 20 starts, and
compared on the maximized likelihood, class shares, every coefficient, the
residual standard deviations and every posterior. Three published analyses are
reproduced to their printed precision.

```
Rscript equivalence/regression-packages/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-regression-packages.R` runs the same comparison under testthat.

**Result (2026-10-29): 76 quantities (3,231 values) compared, 71 agree,
5 known differences (MoEClust's search on one dataset), 0 unexplained.**

| reference | family | cases |
|---|---|---|
| mixtools `regmixEM` (20 random starts, `epsilon = 1e-14`) | Gaussian | `CO2data`, `tonedata`, `NOdata`, flexmix's `NPreg` (`yn ~ x + I(x^2)`) |
| MoEClust `MoE_clust` (model "V", 20 random starts) | Gaussian, with a membership covariate | the same four, and a simulated set with `membership = ~z` (its `gating`) |
| flexmix `FLXMRglm` / `FLXMRglmfix` (20 starts) | Poisson; binomial with classes by group | `NPreg` (`yp ~ x`); the salmonella assay with slopes shared by the classes (`common =`); `betablocker` with classes by Center and a shared Treatment effect |
| SAS/STAT 14.2 PROC FMM Example 40.3 | Poisson | Margolin et al. (1981) salmonella assay (18 rows printed in the documentation); the two-component model of Wang et al. (1996) with shared `dose` and `logd` slopes: -2LL, both intercepts, both slopes, every standard error, the mixing probability and its logit and error; and the one-class model |
| Grun and Leisch (2008), JSS 28(4) | binomial, classes by group | `betablocker`, three classes: log likelihood -159.3605, Treatment -0.2582 |
| Benaglia et al. (2009), JSS 32(6) | Gaussian | `CO2data`: log likelihood -66.93977, residual SDs 2.0493 and 0.8094, share 0.7549 |

Tolerances, fixed before any comparison was run: software references, the log
likelihood to 1e-6 plus 1e-10 of its size and the rest to 1e-4 of
max(1, |value|); published values, half a unit of the last printed digit.

Gaussian flexmix is not used: its M-step divides the residual sum of squares by
n - rank, so its fixed point is not the maximum likelihood estimate (the
harness's flexmix suite replaces that step; see `harness/suites/flexmix.R`).

## Known differences (`known-differences.csv`)

- **tonedata, MoEClust.** The data have several maxima (De Veaux 1989).
  MoEClust stops below the maximum: 141.187 from its default
  hierarchical-clustering start, 132.070 from 20 random starts. mixtools' 20
  random starts and `mixture_regression()` reach 141.198 with the same
  estimates.

Nothing in latents disagreed with any package or published value.
