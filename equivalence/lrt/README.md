# Likelihood-ratio tests: external references

`compare_models(..., test = "lrt")` (twice the log-likelihood difference of
consecutive nested models, its degrees of freedom, and the chi-square p-value)
compared with three independent programs. Run from the project root:

```sh
Rscript equivalence/lrt/make-kit.R   # references -> reference.rds
Rscript equivalence/lrt/compare.R    # latents against them -> comparison.csv
Rscript equivalence/run.R lrt        # the same comparison as a test
```

27 of 27 quantities agree (31 October 2026; R 4.5.2, mclust 6.1.2,
glca 1.4.2, Mplus 9 demo).

| Case | Reference | Models | Compared |
|---|---|---|---|
| Covariance structures | mclust 6.1.2 `Mclust(G = 2, modelNames = )` on `faithful`, EM tolerance 1e-12 | EEI in VVI in VVV; EEI in EEE in VVV | log likelihood of each alternative, LR statistic, df, p |
| Membership covariate | glca 1.4.2 `glca(item(...) ~ 1)` and `~ AGE`, `gofglca(test = "chisq")`, 3 classes, `gss08` complete cases, AGE standardised | `lca()` without and with `profile_covariates = "AGE"` | both log likelihoods, LR statistic from glca's log likelihoods and from its printed deviance, df, p |
| LTA measurement invariance | Mplus 9 demo, `mplus/lta-invariant.inp` (thresholds equal across occasions by label) and `mplus/lta-occasion.inp` (free at each occasion), 1200 persons, 3 binary items, 2 occasions, 2 classes, STARTS = 200 50, best likelihood replicated | `lta(measurement = "invariant")` and `"occasion"` | both log likelihoods, LR statistic, df, log10 p |

## Tolerances

- mclust: log likelihoods 1e-6, statistics 2e-6, p-values 1e-5 relative. With
  its default EM tolerance (1e-5, relative) mclust stops up to 2e-4 short of
  the maximum, so the kit tightens it; the agreement is then about 1e-11.
- glca: log likelihoods 1e-6 (glca's own convergence leaves 3e-7);
  `gofglca()` rounds the deviance to three decimals (0.437 against 0.43653
  from its log likelihoods), so the printed deviance and its p-value are
  compared within 5e-4.
- Mplus prints log likelihoods to three decimals: 5e-4 each, 2e-3 for the
  statistic, 0.01 for log10 p (the p-value is about 1e-23).

## Where the data come from

`mplus/lta.dat` is JStats `validation/r-reference/mplus-demo-lta-classvarying-measurement.dat`
(seeded synthetic data, generated there by `generate-mplus-demo-fixtures.R`).
The two Mplus inputs were written for this kit and run with
`/Applications/MplusDemo/mpdemo <name>.inp <name>.out` inside `mplus/`.

## Not covered externally

- The nesting rule itself (which pairs are refused) has no external
  counterpart; it is unit-tested in `tests/testthat/test-compare-lrt.R`, and
  the sabotage battery checks that those tests fail when every covariance
  structure is taken as nested (`lrt_covariance_nesting`).
- Calibration of the chi-square reference (size at 5% under the null) is not
  simulated here yet.
