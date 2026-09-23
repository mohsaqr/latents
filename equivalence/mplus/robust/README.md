# Genuine Mplus MLR comparison

`ESTIMATOR = MLR` on the same 1,200-individual, 60-group synthetic dataset used
by the original two-level comparison, so the point estimates are identical to
the retained ML run and only the standard errors differ.

| File | Contents |
|---|---|
| `robust.inp` | Mplus specification, identical to `../generated/varying.inp` except for the estimator |
| `robust.out` | Mplus VERSION 9 DEMO (Mac) output, including the H0 scaling correction factor |
| `robust-results.dat` | Saved estimates, robust standard errors, scaling factor and information criteria |
| `robust-posteriors.dat` | Saved class probabilities |
| `synthetic.dat` | Input data, copied from `../generated/` |
| `compare.R` | Reproduces the comparison from the project root |
| `comparison.csv` | Retained numerical result |

| Quantity | Largest absolute difference |
|---|---:|
| Robust standard errors (11 parameters) | 1.82e-7 |
| MLR scaling correction factor | 7.03e-7 |
| AIC, BIC, sample-size adjusted BIC | 4.00e-5 |

The scaling correction factor is computed natively as `tr(A^-1 B) / q`, where
`A` is the observed information and `B` the cross-product of the per-group
scores. Matching Mplus's reported value confirms that definition empirically.
Mplus saves results to eight significant digits, which sets the tolerances.

```sh
Rscript equivalence/mplus/robust/compare.R
```
