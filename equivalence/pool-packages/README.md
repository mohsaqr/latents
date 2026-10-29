# pool_imputations() against mitml and mitools

The harness compares `pool_imputations()` with mice's `pool.scalar()`. This kit
adds two implementations of Rubin's rules that share no code with mice:
`mitml::testEstimates()` and `mitools::MIcombine()`. Each is given latents' own
per-imputation estimates and variances (`get_results(x, "imputations")`), so
only the pooling arithmetic is compared: the pooled estimate, its standard
error, Rubin's degrees of freedom (complete-data df infinite), the relative
increase in variance (mitml only; mitools does not report it) and the
fraction of missing information.

```
Rscript equivalence/pool-packages/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-pool-packages.R` runs the same comparison under testthat.

**Result (2026-10-29): 36 quantities (765 values) compared, 36 agree.**

Cases (five imputations by `mice(method = "pmm", seed = 1)`): `iris` with 10%
of Sepal.Width and Petal.Length missing at random, two and three profiles;
`srl` with 30 values missing in each of three scales, two profiles, and with
test_anxiety as a membership covariate (`profile_covariates`). Parameters the
imputations leave unchanged (no between-imputation variance) are left out:
their relative increase is 0/0 in the references' arithmetic.

Tolerances, fixed before the comparison: 1e-10 of max(1, |value|); 1e-8 for
the degrees of freedom, which reach 1e28 where the imputations barely move a
parameter.
