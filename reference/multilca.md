# Fit a two-level latent class model

Two-level latent class analysis for categorical indicators:
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
with every indicator in `vars` treated as categorical. Each profile is
described by an unrestricted probability for every category of every
item, and the group classes differ in how probable each profile is for
their observations. For models that combine categorical and continuous
indicators, call
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
and name the categorical ones in its `categorical` argument.

## Usage

``` r
multilca(data, vars, id, n_profiles, n_group_classes = 2L, ...)
```

## Arguments

- data:

  A data frame with one row per observation.

- vars:

  Names of the categorical indicator columns. Numeric, integer, logical,
  character and factor columns are accepted; factor levels keep their
  declared order and other types are sorted.

- id:

  Name of the group identifier column, or `NULL` for a single-level
  model.

- n_profiles:

  Number of observation-level latent classes (profiles).

- n_group_classes:

  Number of group-level latent classes.

- ...:

  Further arguments to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  such as `n_starts`, `seed`, `missing`, `tol` or `profile_covariates`.
  `categorical` is set to `vars` and cannot be supplied.

## Value

A fitted model of class `multilpa` (or `multilpa_covariates` when
membership covariates are given), exactly as
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
returns it with `categorical = vars`. The response probabilities are in
`get_results(fit, "responses")`, one row per profile, item and category.

## Conditions

`latents_bad_argument` when `categorical` is supplied, since every
indicator is categorical by definition here. Every condition of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
can also be raised.

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for continuous and mixed indicators, and
[`vignette("lca", package = "latents")`](https://pak.dynasite.org/latents/articles/lca.md).

## Examples

``` r
activities <- c("time_with_friends", "on_social_media", "tv_video_games")
fit <- multilca(subset(student_esm, day <= 1), vars = activities,
                id = "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 1, seed = 1)
get_results(fit, "responses")
#>    profile         indicator category probability threshold
#> 1        1 time_with_friends       no  0.97598029  3.704567
#> 2        2 time_with_friends       no  0.83034660  1.588086
#> 3        1 time_with_friends      yes  0.02401971        NA
#> 4        2 time_with_friends      yes  0.16965340        NA
#> 5        1   on_social_media       no  0.37231410 -0.522302
#> 6        2   on_social_media       no  0.98998036  4.593138
#> 7        1   on_social_media      yes  0.62768590        NA
#> 8        2   on_social_media      yes  0.01001964        NA
#> 9        1    tv_video_games       no  0.73677805  1.029289
#> 10       2    tv_video_games       no  0.86637631  1.869292
#> 11       1    tv_video_games      yes  0.26322195        NA
#> 12       2    tv_video_games      yes  0.13362369        NA
get_results(fit, "profile_probabilities")
#>   group_class profile  probability group_class_probability
#> 1           1       1 1.000000e+00               0.4532333
#> 2           1       2 3.822445e-11               0.4532333
#> 3           2       1 2.195735e-01               0.5467667
#> 4           2       2 7.804265e-01               0.5467667
```
