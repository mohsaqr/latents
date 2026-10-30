# Fit a latent class model

Latent class analysis describes heterogeneity in categorical response
patterns through a specified number of latent classes. Each class is
characterized by the probabilities of the response categories for each
indicator. These probabilities describe how response patterns differ
between classes.

`lca()` estimates the class-specific response probabilities and class
proportions, together with posterior membership probabilities for each
observation. Indicators are assumed independent conditional on class
membership unless associations between selected indicators are
explicitly specified. The posterior probabilities describe uncertainty
in classification; an observation's assigned class is a summary of these
probabilities.

The model treats observations as independent. For observations nested
within groups, use
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md).
Ordinal and count indicators can also be fitted using the structured
measurement models described in
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

## Usage

``` r
lca(data, vars, n_classes, ...)
```

## Arguments

- data:

  A data frame with one row per observation.

- vars:

  Names of the discrete indicator columns. Categories are taken from the
  data. Columns named in `ordinal` or `count` (through `...`) are
  modelled as such; the rest are categorical.

- n_classes:

  Number of classes; a positive integer.

- ...:

  Further arguments for
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  such as `missing`, `n_starts`, `seed`, `min_probability` and
  `profile_covariates`, except `id` and `n_group_classes`.

## Value

A fitted object as returned by
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
usually of class `multilpa`, or `multilpa_covariates` with membership
covariates. Its `call` records the original `lca()` call. Classes are
labelled `profile_1`, `profile_2`, and so on in the returned tables.

## Conditions

`latents_bad_argument` when `id`, `n_group_classes` or `n_profiles` is
passed; the conditions
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
raises otherwise.

## References

Lazarsfeld, P. F., & Henry, N. W. (1968). *Latent Structure Analysis*.
Houghton Mifflin.

## See also

[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) for
continuous indicators;
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
for the two-level model;
[`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
to compare numbers of classes.

## Examples

``` r
fit <- lca(student_esm, c("time_with_friends", "on_social_media",
                          "tv_video_games", "sports"),
           n_classes = 2, n_starts = 3, seed = 1)
fit
#> Latent class analysis: 2 classes
#> 2582 observations; 4 categorical indicators
#> Classes are labelled profile_1, profile_2, ... in every table.
#> Log likelihood: -4560.239231 | AIC: 9138.478 | BIC: 9191.185
#> Converged: TRUE | iterations: 69 | best log likelihood reached by 3 of 3 starts (3 converged)
#> 
#> Class sizes   profile_1 35.6% (919)   profile_2 64.4% (1663)
#> 
#> Response probabilities
#>          indicator category profile_1 profile_2
#>  time_with_friends       no     0.782     0.897
#>  time_with_friends      yes     0.218     0.103
#>    on_social_media       no     1.000     0.552
#>    on_social_media      yes     0.000     0.448
#>     tv_video_games       no     1.000     0.670
#>     tv_video_games      yes     0.000     0.330
#>             sports       no     0.902     0.943
#>             sports      yes     0.098     0.057
#> 
#> Thresholds and standard errors: get_results(x, "responses").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
plot(fit, what = "heatmap")
```
