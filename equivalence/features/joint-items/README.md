# Joint-items feature kit

Local dependence between categorical items
(`lca(..., local_dependence = "a:b")`, the pair fitted as one joint item)
against Latent GOLD 6.1 (`a <-> b | Cluster`), the Mplus 9 Demo (the joined
item as one NOMINAL variable, the same model), poLCA on the joined item, and
an independent EM (`reference.R`, no latents code). Test:
`equivalence/test-features-joint-items.R`; unit and property tests:
`tests/testthat/test-features-joint-items.R`.

## Rerun

From the project root:

```sh
Rscript equivalence/features/joint-items/make-kit.R           # kit/ (deterministic)
JOBS=2 sh equivalence/features/run.sh joint-items all          # returned/
Rscript equivalence/run.R features-joint-items                 # writes report.csv, study.csv
Rscript equivalence/features/joint-items/findings/squarem-early-stop.R
Rscript equivalence/features/joint-items/findings/bootstrap-rare-combination.R
```

Tolerances and the measurements behind them are in the header of
`compare.R`.

## Design (simulab 0.4.10; 54 cases, 62 Latent GOLD runs, 51 Mplus runs)

The conditionally independent items are drawn by simulab; the second item of
each pair is then redrawn from a stored class-specific table
P(second | first, class), and the case keeps that table and the joint table
it implies as its truth.

| Configuration | Cases | k | Items | Pair | n |
|---|---|---|---|---|---|
| bin_small | 5 | 2 | 5 | 2 x 2 | 80 |
| bin_mod | 5 | 2 | 5 | 2 x 2 | 800 |
| cat3_mod | 5 | 3 | 6 | 3 x 3 | 800 |
| bin4_mod | 5 | 2 | 6 | 2 x 4 | 800 |
| twopair_mod | 5 | 3 | 6 | two 2 x 2 pairs | 800 |
| k4_large | 5 | 4 | 7 | 2 x 2 | 5000 |
| huge | 1 | 3 | 8 | 3 x 3 | 20000 (Latent GOLD only) |
| sparse | 4 | 2 | 5 | 3 x 3 with two combinations never observed | 600 |
| perfect | 3 | 2 | 5 | second item equal to the first | 600 |
| boundary | 2 | 2 | 5 | 2 x 3, one category absent in class 1 | 700 |
| partial | 3 | 2 | 5 | 2 x 3, rows missing one item of the pair | 800 |
| indep | 3 | 2 | 5 | 2 x 2, truly independent | 1000 |
| strings / numeric | 2 / 2 | 2 | 5 | string labels; codes {1, 11} x {0, 10, 1} | 600 |
| weights | 2 | 2 | 5 | case weights 1 to 3 (and duplicated rows) | 400 |
| multilevel | 2 | 2 (+ 2 group classes) | 5 | 2 x 2 | 60 groups x 15 |

Every single-level unweighted case except `huge` is also fitted by poLCA.
The reference study adds 140 simulab datasets (bin_small 40, indep 40,
bin4_mod, cat3_mod and twopair_mod 20 each) against the independent EM.

## Results (re-run 2026-10-21)

`Rscript equivalence/run.R features-joint-items`: 83 expectations, 9
failures, 0 errors. Every failure is finding 1: eight are case
`twopair_mod_05`, the ninth is the reference study's seven fits that stop
1.8e-6 to 1.0e-5 short. Unit file: 59 passed.

Largest differences outside `twopair_mod_05` (`report.csv`):

| Quantity | Latent GOLD | Mplus |
|---|---|---|
| log likelihood | 5.0e-5 | 4.2e-4 (8 printed digits) |
| BIC | 5.1e-5 | 4.0e-4 |
| parameter count | 0 (under Latent GOLD's full-grid count) | 0 |
| class sizes | 6.7e-5 | 5.8e-5 |
| item probabilities | 1.3e-4 | |
| joint probabilities | 1.0e-4 | 8.2e-5 |
| posteriors | 4.4e-4 | 4.4e-4 |

Independent reference (`study.csv`, 140 datasets): its likelihood at
latents' estimates equals latents' (relative difference 0); the dependence
table recomputed independently differs by 0; the joined model's likelihood
is never below the plain model's; every parameter count matches.

## Conventions found

- **Unobserved combinations:** Latent GOLD counts the full C1 x C2 grid, so
  it reports k x (unobserved cells) more parameters (perfect dependence: 23
  against 11). The likelihoods agree; counts are compared under Latent
  GOLD's convention.
- **Partial pairs:** Latent GOLD's `includeall` keeps the observed item's
  marginal for a row missing the other item of the pair; latents drops the
  pair for that row. The reference reproduces Latent GOLD's value
  (-2712.4457), and latents matches Latent GOLD run on the data with the
  pair masked.
- **Case weights:** Latent GOLD sums the raw weights; latents rescales them
  to n.
- **Two-level BIC:** Latent GOLD's BIC is latents' `bic_individual`.
- **`predict()` on an unseen combination** raises `latents_bad_data`; a pair
  item constant in the data raises `latents_bad_categorical`.

## Findings

See FINDINGS.md.
