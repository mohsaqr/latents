# Specify prior settings for Gaussian profile estimation

`prior_control()` specifies prior settings that regularise the
estimation of Gaussian profile means and covariance parameters. The
returned object is supplied to the `prior` argument of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
The prior applies to the continuous-indicator measurement parameters;
class-membership probabilities are estimated without a prior.

The settings control mean shrinkage and the prior scale and degrees of
freedom for covariance estimation. Unspecified values are determined
from the data during fitting. The parameterisation follows the
conventions of mclust, with effects that depend on the chosen covariance
structure.

## Usage

``` r
prior_control(shrinkage = 0.01, mean = NULL, dof = NULL, scale = NULL)
```

## Arguments

- shrinkage:

  Non-negative prior precision factor for the means, `kappa` (default
  `0.01`, as in mclust). Zero leaves the means unshrunk.

- mean:

  Optional prior mean, one value per continuous indicator on the scale
  used to fit the profiles. Use original indicator units when no
  transformation is requested, and the transformed scale under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)'s
  `centering` or `scaling`. The separate internal centring used for
  numerical stability is handled automatically. Requires a positive
  `shrinkage`.

- dof:

  Optional prior degrees of freedom, a single positive number.

- scale:

  Optional prior scale: a single positive number for an axis-parallel or
  spherical structure, or a positive-definite `d` by `d` matrix for an
  ellipsoidal one.

## Value

An object of class `latents_prior`: a list with `shrinkage`, `mean`,
`dof` and `scale`, the last three `NULL` where the default is to be
computed from the data. Its
[`print()`](https://rdrr.io/r/base/print.html) method shows these
settings and returns the object invisibly.

## Details

The prior follows Fraley and Raftery (2007) and
[`mclust::priorControl()`](https://mclust-org.github.io/mclust/reference/priorControl.html).
Arguments set to `NULL` use the data-dependent defaults of
[`mclust::defaultPrior()`](https://mclust-org.github.io/mclust/reference/defaultPrior.html):
the prior mean is the vector of indicator means, the degrees of freedom
are `d + 2` for `d` indicators, and the scale is `(1/G)^(2/d)` times the
sample covariance (ellipsoidal structures) or times the average sample
variance (axis-parallel and spherical structures), for `G` profiles.

The updates follow mclust's fitted-model convention, including its
structure-specific departures from a textbook conjugate update, so the
controls have different effects across covariance structures. For `EII`,
`scale` does not enter the variance update. For `EEI`, `shrinkage` and
`mean` enter the variance calculation but do not shrink the fitted
means. For `VEV`, `dof` does not enter the update. `EEI`, `EVI` and
`EEV` count the mean-prior contribution once rather than once per
profile.

## Interpretation

Using a prior changes the estimation procedure from maximum likelihood
to regularised point estimation. It does not provide posterior samples
or Bayesian credible intervals. Ordinary Wald inference is unavailable
for these fits; supported bootstrap inference refits the model with the
prior settings.

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

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for the supported model specifications.

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
