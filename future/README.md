# Deferred features

Code removed from the shipped package so that everything shipped is defensible,
not so that everything is shipped. Nothing here is believed to be wrong; each
was cut because its evidence did not reach the bar the rest of the package
holds to. Each entry says what would have to be true to bring it back.

Removed in 0.12.0.

## `fit_random_intercept()` — `R/random-intercept.R`

A continuous Gaussian group random intercept with unit indicator loadings, as
an alternative to discrete group classes. Adaptive Gauss-Hermite quadrature,
with a check at a higher node count.

Why it is here, in the order the reasons matter:

1. **Independent validation reaches only the one-profile limit.** The Mplus
   comparison covers a single profile. Multiple-profile parameter recovery has
   never been checked against anything.
2. **No inference at all.** `parameter_inference()`, `vcov()` and `confint()`
   refused it with `multilpa_no_inference`. Variance parameters sit on a
   boundary and the integration adds its own error, neither of which had a
   treatment.
3. **It is a different model family.** A continuous group effect is not
   two-level latent profile analysis, which is what the package is named for.
4. It appeared in no vignette.

To restore: establish multiple-profile parameter recovery against an
independent implementation, and implement standard errors with an explicit
treatment of the variance boundary and of the quadrature error. Then move
`R/random-intercept.R` back, restore the blocks in `tests/` below, and put the
`Description:` clause and the README rows back.

Note that `multilpa_no_group_classes` and `multilpa_quadrature_check` were the
conditions this family raised. `multilpa_quadrature_check` was removed from the
catalogue because nothing else raises it. `multilpa_no_group_classes` is kept as
a guard in `R/diagnostics.R` and `R/get-data.R` but is unreachable while every
shipped family has discrete group classes — restoring this family makes it
reachable again.

## `lmr_lrt()` — `R/lmr.R`

The likelihood-ratio statistic and its Lo-Mendell-Rubin small-sample
adjustment. The adjustment factor itself was verified against two genuine Mplus
`TECH11` runs, reproducing the adjusted statistic to the precision Mplus prints
(`747.0633` against `747.063`, adjusted `723.7708` against `723.771`).

Why it is here: it returned a `p_value` column that was always `NA`, because
the Vuong-Lo-Mendell-Rubin reference distribution is not reproduced. The help
page said so plainly, but a function named `_lrt` that performs no test invites
the adjusted statistic to be reported as though it had one, and `bootstrap_lrt()`
already gives a calibrated p-value for the same comparison.

To restore: implement and independently validate the VLMR reference
distribution so a p-value can be returned. An ordinary chi-square approximation
is not an acceptable substitute — the null puts a class on the boundary of the
parameter space and leaves its parameters unidentified.

## Files

| Path | Was |
|---|---|
| `R/random-intercept.R` | `R/random-intercept.R` |
| `R/lmr.R` | `R/lmr.R` |
| `tests/test-random-intercept.R` | `tests/testthat/test-random-intercept.R` |
| `tests/removed-blocks.R` | fourteen `test_that()` blocks lifted out of nine shared test files; each is headed by the file it came from |

`future/` is listed in `.Rbuildignore`, so none of it enters the tarball.
