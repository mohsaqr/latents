# Sampling weights: external and internal references

No mixture or latent-class software output with sampling weights exists in
the JStats fixture set (searched 2026-10-01: 445 Mplus inputs, none with
`WEIGHT` and `TYPE = MIXTURE`). The weighting convention is therefore pinned
two ways.

## 1. Mplus 9, weighted regression (`compare-mplus-complex.R`)

`JStats/validation/r-reference/mplus-demo-complex-weight-on`: 2,000 rows,
`WEIGHT IS wt; TYPE = COMPLEX; y ON x`. A one-profile, full-covariance
weighted `lpa()` of `(y, x)` estimates the weighted mean and covariance, from
which the regression and its conditional log likelihood follow exactly.

| quantity | Mplus | closed-form WLS | latents |
|---|---|---|---|
| slope | 0.892 | 0.8923652 | 0.8923652 |
| intercept | -0.108 | -0.1081710 | -0.1081710 |
| residual variance | 3.779 | 3.7796175 | 3.7796175 |
| log likelihood | -4167.500 | -4167.49988 | -4167.49988 |

latents equals the closed form to 2e-10. The log likelihood matches Mplus,
which confirms the scaling: weights rescaled to sum to n, as in latents
(unscaled weights would give a log likelihood about twice as large). Mplus's
residual variance is printed 6e-4 below the maximum, where the likelihood is
flat; the Mplus gate is 1e-3 for that reason.

## 2. Integer weights equal duplication (tests/testthat/test-weights.R)

For every engine (profile model under each covariance structure, categorical
and FIML, membership covariates, the group-class and cross-level families,
`lta()`, and `mixture_regression()` at all three nestings), integer weights
give the same estimates as the data with each unit repeated that many times,
and the log likelihood times `units / sum(weights)` equals the duplicated
fit's to 1e-8 or better. The weighted analytic scores equal the numeric
gradient of the weighted pseudo log likelihood (1e-5).
