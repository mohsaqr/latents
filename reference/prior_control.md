# Conjugate prior for Gaussian mixture estimation

Requests maximum a posteriori (MAP) estimation of the Gaussian profile
means and covariances under the conjugate prior of Fraley and Raftery
(2007), exactly as
[`mclust::priorControl()`](https://mclust-org.github.io/mclust/reference/priorControl.html)
does for
[`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html).
Pass the result as `multilpa(prior = )`. The prior keeps covariances
away from singularity, which is what lets a many-profile or
many-indicator model be fitted where maximum likelihood degenerates.

## Usage

``` r
prior_control(shrinkage = 0.01, mean = NULL, dof = NULL, scale = NULL)

# S3 method for class 'latents_prior'
print(x, ...)
```

## Arguments

- shrinkage:

  Non-negative prior precision factor for the means, `kappa`; mclust's
  default `0.01`. Zero leaves the means unshrunk.

- mean:

  Optional prior mean, one value per continuous indicator, in the
  indicators' own units. Naming one needs a positive `shrinkage`.

- dof:

  Optional prior degrees of freedom, a single number.

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

## Conditions

Raises `latents_bad_argument` for a negative or non-finite `shrinkage`,
a `mean` with a zero `shrinkage`, or a non-positive `dof` or `scale`.

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
