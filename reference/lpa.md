# Fit a latent profile model

Latent profile analysis describes heterogeneity in continuous indicators
through a specified number of latent profiles. Each profile represents a
pattern of indicator means and within-profile variation. `lpa()`
estimates these patterns using a finite mixture of multivariate normal
distributions.

By default, indicators are independent within profiles and their
variances differ across profiles. Alternative covariance structures
allow associations between indicators within profiles and equality
constraints across profiles. The fitted model provides profile means,
variances, proportions, and posterior membership probabilities for each
observation, together with model-fit and classification summaries.

The model treats observations as independent. For observations nested
within groups, use
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
Additional measurement options support mixed continuous and discrete
indicators.

## Usage

``` r
lpa(data, vars, n_profiles, ...)
```

## Arguments

- data:

  A data frame with one row per observation.

- vars:

  Names of the indicator columns. Numeric columns are continuous
  indicators; categorical ones are named in `categorical`.

- n_profiles:

  Number of profiles; a positive integer.

- ...:

  Further arguments for
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md):
  the covariance structure (`model`, or `variance_model` and
  `covariance_model`, or `volume`, `shape` and `orientation`),
  `missing`, `categorical`, `n_starts`, `seed`, `boundary`, `prior`,
  `noise` and `profile_covariates`. Arguments `id` and `n_group_classes`
  are not accepted.

## Value

A fitted object as returned by
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
usually of class `multilpa`, or `multilpa_covariates` with membership
covariates. Its `call` records the original `lpa()` call.

## Conditions

`latents_bad_argument` when `id` or `n_group_classes` is passed; the
conditions
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
raises otherwise.

## References

Oberski, D. (2016). Mixture models: Latent profile and latent class
analysis. In J. Robertson & M. Kaptein (Eds.), *Modern Statistical
Methods for HCI* (pp. 275–287). Springer.

## See also

[`lca()`](https://pak.dynasite.org/latents/reference/lca.md) for
categorical indicators;
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for the two-level model;
[`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
to compare numbers of profiles.

## Examples

``` r
fit <- lpa(srl, c("cognitive_strategies", "intrinsic_value",
                  "self_efficacy", "self_regulation", "test_anxiety"),
           n_profiles = 2, variance_model = "equal",
           covariance_model = "full", n_starts = 3, seed = 1)
fit
#> Latent profile analysis: 2 profiles
#> 300 observations; equal full residual covariance (EEE)
#> Log likelihood: -600.924999 | AIC: 1253.850 | BIC: 1350.148
#> Converged: TRUE | iterations: 18 | best log likelihood reached by 3 of 3 starts (3 converged; emEM: 15 screened)
#> Smallest standardized Gaussian eigenvalue: 0.1025 (profiles, profile_1); Gaussian mass: 93.07; weight effective n: 111.8.
#> 
#> Class sizes   profile_1 31.0% (93)   profile_2 69.0% (207)
#> 
#> Means
#>             indicator profile_1 profile_2
#>  cognitive_strategies     5.215     5.928
#>       intrinsic_value     5.516     6.386
#>         self_efficacy     5.392     6.167
#>       self_regulation     5.108     5.740
#>          test_anxiety     3.845     3.635
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
plot(fit, what = "raincloud")
```
