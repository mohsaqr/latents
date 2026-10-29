# Pair-covariance feature kit

Covariances for chosen pairs of continuous indicators
(`lpa(..., local_dependence = "y1:y2")`) against Latent GOLD 6.1, the Mplus
9 Demo, mclust and an independent implementation written from the formula
(`reference.R`, no latents code). Test:
`equivalence/test-features-pair-covariance.R`; unit and property tests:
`tests/testthat/test-features-pair-covariance.R`.

## Rerun

From the project root:

```sh
Rscript equivalence/features/pair-covariance/make-kit.R      # kit/ (deterministic)
JOBS=2 sh equivalence/features/run.sh pair-covariance all     # returned/: .lst (LG), .out (Mplus)
Rscript equivalence/run.R features-pair-covariance            # writes results/*.csv
Rscript equivalence/features/pair-covariance/findings/asymmetric-covariance.R
```

The comparisons read the retained outputs offline; latents and the 120
reference-study datasets are refitted on two R workers.

## Design (simulab 0.4.10; 62 cases, 122 reference runs)

Each profile's correlation matrix is nonzero on the fitted pairs only (the
misspecified cases add one correlation outside them). Latent GOLD fits
`y1 <-> y2 | Cluster` (without `| Cluster` for one shared covariance); Mplus
`y1 WITH y2;` in each class section (in `%OVERALL%` for a shared one).

| Configuration | Model | n | k | Pairs | Cases |
|---|---|---|---|---|---|
| vvi_one_k2_m | varying, 4 items | 500 | 2 | y1:y2 | 5 |
| vvi_two_k3_m | varying, 5 items | 600 | 3 | y1:y2, y3:y4 | 5 |
| eei_one_k2_m | equal, 4 items | 500 | 2 | y1:y2 | 5 |
| eei_two_k3_m | equal, 5 items | 600 | 3 | y1:y2, y4:y5 | 5 |
| vvi_two_k4_m | varying, 5 items | 1000 | 4 | y1:y3, y2:y5 | 5 |
| vvi_one_k2_s | varying, 3 items | 60 | 2 | y2:y3 | 5 |
| eei_two_k2_s | equal, 4 items | 80 | 2 | y1:y2, y3:y4 | 5 |
| vvi_one_k3_l | varying, 4 items | 5000 | 3 | y1:y2 | 5 |
| near_one | correlation +-0.98 | 500 | 2 | y1:y2 | 2 |
| zero_cov | true covariance 0 | 500 | 2 | y1:y2 | 2 |
| floor_unpaired | floor 0.01 binds on an unpaired variance | 400 | 2 | y1:y2 | 2 |
| floor_paired | floor 1e-3 binds on the pair's eigenvalue | 400 | 2 | y1:y2 | 2 |
| small_class | class shares 0.88 / 0.09 / 0.03 | 400 | 3 | y1:y2 | 2 |
| scales | y1 in units of 1e-4, y2 of 1e4 (`min_variance = 1e-14`) | 500 | 2 | y1:y2 | 2 |
| two_vvi / two_eei | two indicators, also against VVV / EEE | 500 | 2 | y1:y2 | 2 + 2 |
| missing | MAR missing, FIML | 500 | 2 | y1:y2 | 3 |
| misspecified | y3:y4 correlated but not fitted | 500 | 2 | y1:y2 | 2 |
| huge | | 20000 | 2 | y1:y2 | 1 |

Reference runs: Latent GOLD 58 (100 or 200 start sets), Mplus 60
(`STARTS = 200 40`, 400 80 for the 4-class and small-class cases), and 4 more
Mplus runs that hold the binding floor as an equality. The reference study
adds 120 simulab datasets (n 40-2000, k 2-3, 2-5 items, 1-2 pairs, 60 equal
and 60 varying, 12 with missing values) checked against `reference.R` alone.

## Results (re-run 2026-10-29, latents fitting at tol 1e-12)

`equivalence/test-features-pair-covariance.R` on 0.98.1: 1291 expectations,
0 failures, 0 errors. latents now fits at Latent GOLD's tolerance, 1e-12 (it
was 1e-10): at 1e-10 the flat case vvi_one_k3_l_02 stopped with posteriors
2.0e-5 from Mplus, just over the 2e-5 bound; at 1e-12 within 1.0e-6.

| Against | Quantity | Largest difference (1e-10, 2026-10-21) |
|---|---|---|
| Latent GOLD (58 runs) | log likelihood | 1.7e-9 (1.9e-8) |
| | posteriors | 2.4e-6 (1.6e-5) |
| | BIC | 4.9e-5 (4.9e-5) |
| | estimates (standard-error units) | 3.6e-5 (1.7e-4) |
| | standard errors (relative), complete data | 1.0e-6 (6.2e-6) |
| | standard errors (relative), missing data | 6.3e-2 (a convention, below) |
| | parameter count | equal in all 58 |
| Mplus (56 runs reaching latents' maximum) | log likelihood | 3.7e-3 at n = 20000, 5.0e-4 at n = 5000 (its 8-digit print) |
| | posteriors | 1.1e-5 (1.7e-5) |
| | estimates (SE units) / SEs (relative) | 6.0e-5 / 1.0e-5 (1.7e-4 / 1.1e-5) |
| Mplus, floor held as an equality (4) | log likelihood | 4.4e-5 |
| Independent reference (120 datasets) | likelihood at latents' estimates (relative) | 2.0e-15 (7.0e-15) |
| | its own maximum above latents' | 2.9e-8 (1.3e-8) |
| | Newton decrement at latents' estimates | 5.8e-8 (2.7e-8) |
| | standard errors (relative) | 8.0e-6 (4.9e-6) |
| | entries outside the pattern | exactly 0 in all 120 |
| | exact symmetry | all 120 (largest asymmetry 0) |
| | parameter count against the formula | equal in all 120 |
| Floor solver (200 random problems, binding in 179) | KKT stationarity | 2.6e-11 |

mclust (`Mclust` and `me` from latents' posteriors) agrees on the
two-indicator cases, where a pair covariance is a full covariance.

## Conventions found

- **Latent GOLD's standard errors with missing values** differ from the
  observed-data Hessian of the FIML likelihood by up to 6.3% (`missing_01`
  to `_03`). latents matches the reference's Hessian of that likelihood to
  1e-4, so there it is compared with the reference instead.
- **Mplus refuses the 1e-4-scale indicator** ("One or more variables have a
  variance of zero"); the two `scales` cases are compared with Latent GOLD
  only. The default `min_variance` (1e-6) is an absolute bound and binds at
  that scale, so the kit fits them with 1e-14.
- **Mplus's inequality-constrained floor** runs stop 1.6e-4 to 1.57 below
  latents (three of the four not replicated). Held as an equality (one
  parameter fewer), Mplus agrees to 4.4e-5.
- **numDeriv's Richardson Hessian** is wrong by up to 13% on a perfectly
  separated class; the closed-form value agrees with latents, so the
  reference uses plain central differences.
- The floor solver of a40d5c7 is not used by pair fits: their floor is an
  eigenvalue clip, block by block.

## Findings

See FINDINGS.md.
