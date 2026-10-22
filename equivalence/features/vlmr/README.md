# VLMR feature kit

`get_results(grid, "lmr")` (Vuong 1989; Lo, Mendell and Rubin 2001; Vermunt
2024) against Latent GOLD 6.1, the Mplus 9 Demo's TECH11 and an independent R
reference. Test: `equivalence/test-features-vlmr.R`; unit and property tests:
`tests/testthat/test-features-vlmr.R`.

## Rerun

From the project root:

```sh
Rscript equivalence/features/vlmr/make-kit.R           # kit/ and kit/cases.rds (deterministic)
JOBS=2 sh equivalence/features/run.sh vlmr all          # returned/: .lst (LG) and _k<K>.out (Mplus)
Rscript equivalence/run.R features-vlmr                 # the equivalence test (fits live)
Rscript equivalence/features/vlmr/summarise.R           # summary.csv, imhof-stress.csv, comparisons.csv
NOT_CRAN=true Rscript equivalence/features/vlmr/sabotage.R
Rscript equivalence/features/vlmr/repro-missed-maximum.R
```

`VLMR_FIT_CACHE=<dir>` keeps fitted grids between runs; unset, every grid is
fitted live (about 50 minutes on 2 workers). Each case is fitted with seed 1
(60 starts, 30 for n = 5000, tol 1e-12); if a candidate stays below Latent
GOLD's maximum, seeds 2 and 3 with 150 starts are tried and the grid
reaching the most maxima is kept.

## Design (simulab 0.4.10; 107 datasets, 271 neighbouring pairs)

62 datasets run through both programs (Latent GOLD `Cluster nominal 1:4` or
`2:4`, 200 start sets; Mplus one input per K with `ESTIMATOR = MLR`,
`STARTS = 200 40`, `TECH11`), 45 reference-only:

| Configuration | Model | n | True k | Cases |
|---|---|---|---|---|
| eei_k2 | EEI, 4 items, well separated | 500 | 2 (null true for 2 vs 3) | 5 |
| eei_k3poor | EEI, poorly separated | 500 | 3 | 5 |
| eei_small | EEI, 3 items | 80 | 2 | 5 |
| eei_large | EEI, 3 items, range 2:4 | 5000 | 3 | 5 |
| eei_k1 | one class (simulate_correlated) | 300 | 1 (1 vs 2 under the null) | 5 |
| vvi_k2 / vvi_k3 | VVI | 500 | 2 / 3 | 5 / 5 |
| eee_k2 | EEE, shared correlation | 500 | 2 | 5 |
| lca_k2 / lca_k3 | LCA, 5 / 6 binary items | 500 / 1000 | 2 / 3 | 5 / 5 |
| weights | EEI, weights 0.5 / 1.5 summing to n | 400 | 2 | 3 |
| missing | EEI, MCAR 10% (FIML); case 3 with 20 empty rows | 500 | 2 | 3 |
| vvi_small | VVI, over-extracted near the boundary | 60 | 2 | 3 |
| eei_tiny | a 3% class far from the rest | 500 | 3 | 3 |
| ref_* (reference only, range 1:3) | EEI n 50-4000, VVI n 100-2000, LCA n 150/3000, EEE n 300, one class n 2000, VVI with random sampling weights | 50-4000 | 1-3 | 45 |

## References

- **Latent GOLD**: LL, LR, adjusted LR, both p-values, mean and SD of the
  reference distribution (4 decimals). The reference evaluated at Latent
  GOLD's own printed estimates, and at latents' estimates jittered by
  +-5e-5, measures how far estimates agreeing to 4 decimals move the weights
  (up to 0.02 for over-extracted 4-class models); that movement is the only
  allowance beyond printing.
- **Mplus TECH11**: H0 LL, LR, parameter difference, adjusted LR (3 decimals).
- **Independent R reference** (`reference.R`, no latents code): own
  likelihoods (last-class logits, log variances, log-Cholesky, item logits),
  numDeriv scores and Hessian in coordinates scaled to each parameter's
  natural size, Vuong's W from the published formula, CompQuadForm Imhof /
  Davies / Farebrother and 1e6 Monte Carlo draws.
- **Imhof stress** (`vlmr_imhof_stress()`): 460 weight vectors against exact
  chi-square, F and convolution probabilities and Davies.

## Results (measured 2026-10-21)

Equivalence test (`Rscript equivalence/run.R features-vlmr`, live fits): the
first 10 failures are the missed maxima below, then testthat's failure limit
stops the run. With the failure limit lifted (cached grids from the same
code): 5318 expectations, 13 failures, all of them finding 1. Unit file:
176 passed, 0 failed (10 s). Sabotage: the chi-square tail, unsigned
weights and an `n` for `log n` adjustment fail 51, 8 and 2 unit expectations.

| Reference | Pairs | Quantity | Max difference |
|---|---|---|---|
| Latent GOLD | 122 computed at shared maxima | LL | 4.9e-5 |
| | | LR / adjusted LR (Latent GOLD's n) | 4.9e-5 / 5.0e-5 |
| | | VLMR p / aLMR p | 4.4e-4 / 4.5e-4 (95% of pairs within 5e-5, median 2.4e-6) |
| | | mean / SD | 5.5e-3 / 5.1e-3 (81% within 5e-5; all within the measured estimate movement) |
| | 12 LCA computed | VLMR p / mean | 8.3e-7 / 5.9e-4 |
| | 9 weighted | p / mean of the frequency-weight reference | 2.7e-5 / 7.0e-4 (latents' sampling-weight mean differs by 1.99, by design) |
| Mplus | 130 of 131 pairs where Mplus reaches both maxima | H0 LL / H1 LL / LR / adjusted LR | 4.9e-4 / 4.9e-4 / 5.0e-4 / 5.0e-4 |
| | | VLMR p (not tested) | max 0.41, median 5.5e-3 |
| | | mean (not tested) | max 41.5, median 1.99 |
| Independent reference | 241 pairs, 107 datasets | LR (relative) | 2.8e-12 |
| | | weights (eigenvalue gap) | 4.4e-6 |
| | | mean / SD | 6.8e-6 / 5.3e-6 |
| | | VLMR p vs Imhof / Davies; aLMR p vs Davies | 7.4e-7 / 7.4e-7; 7.6e-7 |
| | | VLMR p vs 1e6 Monte Carlo | 2.4 binomial SEs |
| Imhof stress | 460 vectors | exact (chi-square, F, convolution, 1e-8..1e4 range, tiny/huge q) | 3.2e-13 |
| | | 50 random weights vs Davies | 1.1e-12 |
| | | joint scaling by 1e-6 / 1e6 | 1.1e-16 |

Farebrother never applied: every fitted pair had at least one negative
weight.

## Conventions found

- Latent GOLD's `caseweight` is a frequency weight (B = sum w s s'); latents'
  weights are sampling weights (B = sum w^2 s s'). A frequency-weight
  reference reproduces Latent GOLD (p within 2.7e-5).
- Latent GOLD counts rows with nothing observed in the adjustment's n
  (`missing_03`: 500), latents counts informative units (480).
- Latent GOLD's degenerate VVI solutions (20 pairs: a class of one or two
  rows with variances printed 0.0000) exceed latents' likelihood, which its
  variance floor keeps regular; these pairs are not compared.
- latents refuses 17 of 29 LCA pairs (bound-active response probabilities;
  note in the table) where Latent GOLD prints a VLMR.
- Mplus TECH11: same LR, H0 and adjusted LR; its mean, SD and p-values come
  from another construction of the weights, not reproduced (W with observed,
  OPG or mixed information, single blocks or no cross block were tried on
  `eee_k2_01`). Its H0, re-estimated from the K model without class 1, is
  far below the K-1 maximum in several VVI pairs (by 162 at most).

## Findings

1. **Random starts miss regular maxima (latents).** 9 candidates of 8
   EEI/EEE datasets stay 0.15 to 4.85 below a non-degenerate maximum that
   Latent GOLD and Mplus both reach (`eei_k2_02` k3, `eei_k3poor_02` k4,
   `eei_k3poor_04` k3, `eei_large_02` k4, `eei_k1_01` k4, `eei_k1_02` k3,
   `eei_k1_05` k4, `eee_k2_04` k3 and k4). Repro:
   `repro-missed-maximum.R` (`eee_k2_04`, 3 classes: latents -2222.2043 from
   60, 150 and 600 starts; Latent GOLD -2220.9820; latents started at Latent
   GOLD's estimates converges to -2220.9822). The VLMR computation at a given
   maximum agrees with the reference to 1e-6; the LRs of these pairs inherit
   the wrong maximum. Not fixed, on purpose: the start strategy is shared by
   every released model, whose results were checked against Mplus, Latent
   GOLD and mclust (the same limitation as heavy_02 in
   `../direct-effects/FINDINGS.md`). Re-run 2026-10-21 after the session's
   other fixes: 5318 expectations, the same 13 failures.
2. The weights of over-extracted models are sensitive to the fourth decimal
   of the estimates (0.02 in the mean), so Latent GOLD's printed moments
   cannot be checked to 5e-5 there; the test allows the measured movement.
3. Scaling the indicators by 1e-4 changes the weights by up to 3e-5
   (relative) through EM's endpoint (estimates 5e-6 apart); the reference
   agrees with latents at both endpoints. A 1e6 shift is exact to 1e-9.

Not done: multilevel (group-class) comparisons have no external or
independent reference here.
