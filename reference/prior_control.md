# Set a prior that regularises the means and covariances of Gaussian profile models

`prior_control()` specifies prior settings that regularise the
estimation of Gaussian profile means and covariance parameters. The
returned object is supplied to the `prior` argument of
[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) or
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
The prior applies to the continuous-indicator measurement parameters;
class-membership probabilities are estimated without a prior.

A prior is used when maximum likelihood gives a profile whose variance
collapses towards zero, or when a profile rests on few observations. It
pulls the variances towards a common scale and, with a positive
`shrinkage`, the profile means towards a common mean.

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
  `0.01`, as in mclust). A larger value pulls the profile means more
  strongly towards the prior `mean`; zero leaves the means unshrunk and
  keeps only the prior on the covariances.

- mean:

  Prior mean. `NULL` (default) uses the mean of each indicator.
  Otherwise one value per continuous indicator on the scale used to fit
  the profiles: original indicator units when no transformation is
  requested, and the transformed scale under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)'s
  `centering` or `scaling`. The separate internal centring used for
  numerical stability is handled automatically. Requires a positive
  `shrinkage`.

- dof:

  Prior degrees of freedom. `NULL` (default) uses `d + 2` for `d`
  indicators. Otherwise a single positive number; a larger value gives
  the prior `scale` more weight against the data.

- scale:

  Prior scale of the covariance. `NULL` (default) computes it from the
  data, as described under *Details*. Otherwise a single positive number
  for a structure without covariances (codes ending in `I`), or a
  positive-definite `d` by `d` matrix for a structure with covariances
  (the other codes).

## Value

An object of class `latents_prior`: a list with `shrinkage`, `mean`,
`dof` and `scale`, the last three `NULL` where the default is to be
computed from the data. Its
[`print()`](https://rdrr.io/r/base/print.html) method shows these
settings and returns the object invisibly. The object holds settings
only; nothing is estimated until it is passed to a fitting function.

## Details

The prior follows
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

## Where a prior applies

A prior is available for the covariance structures `"EII"`, `"VII"`,
`"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`, `"EEV"`, `"VEV"` and
`"VVV"`; mclust defines none for `"VEE"`, `"EVE"`, `"VVE"` and `"EVV"`.
The indicators must all be continuous and complete. A prior is not
combined with `categorical`, `ordinal` or `count` indicators, missing
values, `fixed`, membership covariates or sampling weights, and
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) has no
`prior` argument.
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
refuses these combinations with class `latents_unsupported_prior`
(`latents_unsupported_indicator` for ordinal or count indicators,
`latents_unsupported_weights` for weights).

## Conditions

`latents_bad_argument` for a negative or non-finite `shrinkage`, a
non-finite `mean`, a `mean` with a zero `shrinkage`, a non-positive or
non-finite `dof`, or a `scale` that is neither a positive number nor a
square, symmetric, positive-definite matrix. A `mean` or `scale` whose
size does not match the indicators or the structure is refused when the
model is fitted, also with `latents_bad_argument`.

## References

Fraley, C., & Raftery, A. E. (2007). Bayesian regularization for normal
mixture estimation and model-based clustering. *Journal of
Classification*, 24, 155-181.
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

# Used in a fit:
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           prior = prior_control(), n_starts = 2, seed = 1)
fit
#> Latent profile analysis: 2 profiles
#> 272 observations; varying diagonal residual covariance (VVI)
#> Posterior mode under mclust's conjugate prior; the log likelihood and criteria below are unpenalized
#> Log likelihood: -1229.215687 | AIC: 2476.431 | BIC: 2508.884
#> Converged: TRUE | iterations: 15 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.1746 (profiles, profile_1); Gaussian mass: 172.7; weight effective n: 174.2.
#> 
#> Class sizes   profile_1 63.5% (173)   profile_2 36.5% (99)
#> 
#> Means
#>  indicator profile_1 profile_2
#>  eruptions     4.301     2.074
#>    waiting    80.187    54.741
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
