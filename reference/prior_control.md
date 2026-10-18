# Specify a conjugate prior for Gaussian mixture estimation

Returns a `latents_prior` object requesting maximum a posteriori (MAP)
estimation of the Gaussian profile means and covariances under the
conjugate prior of Fraley and Raftery (2007), as
[`mclust::priorControl()`](https://mclust-org.github.io/mclust/reference/priorControl.html)
does for
[`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html).
Pass the result as `multilpa(prior = )`. The prior keeps the covariance
estimates away from singularity.

## Usage

``` r
prior_control(shrinkage = 0.01, mean = NULL, dof = NULL, scale = NULL)

# S3 method for class 'latents_prior'
print(x, ...)
```

## Arguments

- shrinkage:

  Non-negative prior precision factor for the means, `kappa` (default
  `0.01`, as in mclust). Zero leaves the means unshrunk.

- mean:

  Optional prior mean, one value per continuous indicator, on the scale
  the profile means are estimated on: the indicators' own units, or
  their centred values under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)'s
  `centering`. Requires a positive `shrinkage`.

- dof:

  Optional prior degrees of freedom, a single positive number.

- scale:

  Optional prior scale: a single positive number for an axis-parallel or
  spherical structure, or a positive-definite `d` by `d` matrix for an
  ellipsoidal one.

- x:

  A `latents_prior` object.

- ...:

  Ignored.

## Value

An object of class `latents_prior`: a list with `shrinkage`, `mean`,
`dof` and `scale`, the last three `NULL` where the default is to be
computed from the data.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## Details

Every argument left `NULL` takes
[`mclust::defaultPrior()`](https://mclust-org.github.io/mclust/reference/defaultPrior.html)'s
value, computed from the indicators being fitted: the prior mean is
their column means, the degrees of freedom are `d + 2` for `d`
indicators, and the scale is `(1/G)^(2/d)` times the sample covariance
(ellipsoidal structures) or times the average sample variance
(axis-parallel and spherical structures), for `G` profiles.

The updates follow mclust's fitted-model convention, including its
structure-specific departures from a textbook conjugate update, so the
controls have different effects across covariance structures. For `EII`,
`scale` does not enter the variance update. For `EEI`, `shrinkage` and
`mean` enter the variance calculation but do not shrink the fitted
means. For `VEV`, `dof` does not enter the update. `EEI`, `EVI` and
`EEV` count the mean-prior contribution once rather than once per
profile.

## Conditions

`latents_bad_argument` for a negative or non-finite `shrinkage`, a
non-finite `mean`, a `mean` with a zero `shrinkage`, a non-positive or
non-finite `dof`, or a `scale` that is neither a positive number nor a
square, symmetric, positive-definite matrix.

## References

Fraley, C., & Raftery, A. E. (2007). Bayesian regularization for normal
mixture estimation and model-based clustering. *Journal of
Classification*, 24, 155–181.
[doi:10.1007/s00357-007-0004-5](https://doi.org/10.1007/s00357-007-0004-5)

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
whose `prior` argument takes this.

## Examples

``` r
prior_control()
#> Conjugate prior (Fraley & Raftery, 2007)
#>   shrinkage: 0.01
#>   mean:      mclust default (from the data)
#>   dof:       mclust default (from the data)
#>   scale:     mclust default (from the data)
prior_control(shrinkage = 0)
#> Conjugate prior (Fraley & Raftery, 2007)
#>   shrinkage: 0
#>   mean:      mclust default (from the data)
#>   dof:       mclust default (from the data)
#>   scale:     mclust default (from the data)
```
