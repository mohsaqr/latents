# ML three-step feature kit

The ML step-three model for a distal outcome
(`three_step(fit, data, "y", method = "ml")`) against Latent GOLD 6.1
(`step3 modal ml` on latents' own posteriors, and end to end from Latent
GOLD's own step one), the Mplus 9 Demo (DU3STEP and DE3STEP) and an
independent implementation (`reference.R`: BFGS and Newton on its own
likelihood, numDeriv information, its own sandwich; no latents code). Test:
`equivalence/test-features-ml-three-step.R`; unit and property tests:
`tests/testthat/test-features-ml-three-step.R`.

## Rerun

From the project root:

```sh
Rscript equivalence/features/ml-three-step/make-kit.R               # kit/ and kit/cases.rds
JOBS=2 sh equivalence/features/run.sh ml-three-step all             # Latent GOLD and Mplus
Rscript equivalence/features/ml-three-step/make-pipeline-step3.R    # LG end to end, step three
JOBS=2 sh equivalence/features/run.sh ml-three-step latentgold _s3
Rscript equivalence/features/ml-three-step/make-random-starts.R     # LG with 100 random start sets
Rscript equivalence/run.R features-ml-three-step                    # writes maxima.csv, comparisons.csv
NOT_CRAN=true Rscript equivalence/features/ml-three-step/sabotage.R
```

The full equivalence run takes about 35 minutes.

## Design (simulab 0.4.10; 73 kit cases, 120 sweep datasets)

| Configuration | Family | k | n | Separation | Outcome | Cases |
|---|---|---|---|---|---|---|
| lpa_k2_hi / lpa_k2_lo | LPA | 2 | 500 | 2.2 / 1.2 | varying / equal | 5 + 5 |
| lpa_k3_hi / lpa_k3_lo | LPA | 3 | 500 / 600 | 2.0 / 1.15 | varying | 5 + 5 |
| lpa_k4_hi | LPA | 4 | 800 | 2.2 | equal | 5 |
| lpa_k5_hi | LPA | 5 | 1500 | 2.3 | varying | 5 |
| lca_k2_hi / lca_k3_lo | LCA | 2 / 3 | 500 / 800 | 0.85 / 0.70 | varying / equal | 5 + 5 |
| small_k2 / small_k3 | LPA | 2 / 3 | 60 / 100 | 2.2 / 2.4 | varying / equal | 5 + 5 |
| large_lca_k3 | LCA | 3 | 5000 | 0.85 | varying | 5 |
| huge_lpa_k2 | LPA | 2 | 20000 | 1.3 | varying | 1 |
| edge cases | | | | | offset 1e6, units 1e-4, skewed, t(3), null, perfect classification, a class no unit is modal for (2), constant outcome, NA outcome, weighted fit | 11 |
| two-level | VVI, 2 x 2 group classes | 2 | 100 groups | | individual level with `psuid`, group level | 6 |

Relative entropy runs from 0.38 to 1 (perfect classification). Mplus runs on 10 regular cases; Latent GOLD
end to end on 6. The sweep crosses LPA/LCA, k = 2-4, n = 50, 200, 1000 and
3000, high and low separation, and five outcome kinds (varying, equal,
skewed, null, heavy-tailed); it is compared with the reference only.

## Results (re-run 2026-10-21)

`Rscript equivalence/run.R features-ml-three-step`: 2011 expectations, 2
failures, 0 errors. Both failures are sweep_108 under equal variances
(FINDINGS.md, finding 1). Unit file: 82 passed; the two existing ML files
(`test-three-step-ml.R` 69, `test-three-step.R` 74, pinned to Latent GOLD)
pass unchanged.

Largest differences (`maxima.csv`):

| Against | Quantity | Largest difference |
|---|---|---|
| Latent GOLD `step3 modal ml`, latents' posteriors (139 fits) | log likelihood | 1.1e-9 |
| | means (SE units) | 5.5e-6 |
| | mean SE (relative) | 1.1e-7 |
| | pair Wald (relative) | 8.4e-7 |
| | variances (relative) | 2.2e-7 |
| Latent GOLD with `psuid` (6) | clustered mean SE (relative) | 7.9e-9 |
| Latent GOLD end to end (10) | means (SE units) / SE (relative) | 1.9e-5 / 6.6e-6 |
| Mplus DU3STEP / DE3STEP (20) | means (SE units) | 1.0e-2 (3 printed decimals) |
| | SE without CR1 (relative) | 3.3e-4 |
| Independent reference, kit (140 fits) | log likelihood | 6.1e-10 |
| | means / 90% CI (SE units) | 3.2e-7 / 3.2e-7 |
| | mean SE (relative) | 3.6e-8 |
| | effective n (relative) | 5.7e-10 |
| Independent reference, sweep (238 fits) | log likelihood (outside sweep_108) | 4.4e-11 |
| | means / 90% CI (SE units) | 4.4e-7 / 4.4e-7 |
| | mean SE / pair SE (relative) | 3.8e-8 / 3.7e-8 |
| | single-unit class (sweep_041): reference sandwich / model SE | 3.5e-11 (latents reports NA) |
| Invariance (lpa_k2_hi_01) | outcome + 1e6, outcome x 1e4 (SE units) | 1.7e-11, 0 |

## Conventions found

- Latent GOLD's robust standard errors include N/(N-1), the same CR1 factor
  as latents; Mplus's do not (compared with the factor removed).
- Mplus fixes the class sizes, as latents does; the free-size model misses
  by 0.1 to 10 SE (measured when the kit was built; the test keeps the
  fixed-size comparison). Its "Last Step" proportions equal latents'
  `effective_n / n`.
- One Latent GOLD step-one fit (`lpa_k3_lo_01`) is a variance spike (a class
  on one unit) and is left out of the end-to-end comparison.
- Weighted fits are refused (`latents_unsupported_weights`), so integer
  weights against duplicated rows cannot be tested; the refusal is asserted.

## Findings

See FINDINGS.md.
