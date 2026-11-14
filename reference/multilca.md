# Fit a two-level latent class model

`multilca()` fits a latent class model for categorical data with
observations nested within groups. It describes response patterns
through latent classes and differences between groups through the
relative prevalence of those classes. Groups with similar class
compositions are represented by the same latent group class.

By default, the response probabilities defining the observation-level
classes are shared across groups. The model estimates these
probabilities together with class proportions and posterior membership
probabilities, allowing both observations and groups to be classified
while retaining information about classification uncertainty.

## Usage

``` r
multilca(data, vars, id, n_profiles, n_group_classes = 2L, ...)
```

## Arguments

- data:

  A data frame with one row per observation.

- vars:

  Names of the indicator columns. Numeric, integer, logical, character
  and factor columns are accepted; factor levels retain their declared
  order and other types are sorted.

- id:

  Name of the group identifier column, or `NULL` for a single-level
  model.

- n_profiles:

  Number of observation-level latent classes (profiles).

- n_group_classes:

  Number of group-level latent classes. Defaults to 2, or to 1 when
  `id = NULL`.

- ...:

  Further arguments to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  such as `n_starts`, `seed`, `missing`, `tol` or `profile_covariates`.
  `ordinal` and `count` name ordinal and count indicators; the remaining
  indicators are categorical. `categorical` is accepted but has no
  further effect; its names must be in `vars`.

## Value

A fitted `multilpa` object (`multilpa_covariates` when membership
covariates are specified), as returned by
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
Categorical response probabilities are in
`get_results(fit, "responses")`, one row per profile, item and category.

## Details

The function calls
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
treating all indicators in `vars` as categorical except those specified
through `ordinal` or `count` in `...`. Categorical indicators have
unrestricted class-specific response probabilities; ordinal and count
indicators use the corresponding structured measurement models. For a
combination of continuous and categorical indicators, use
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
with `categorical`.

## Conditions

`latents_bad_categorical` when `categorical` names a column not in
`vars`; otherwise the conditions of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

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
#> Warning: The fitted likelihood still carries a score worth 1.6% of a standard error; refit with a tighter `tol` before using Wald inference.
#>    profile         indicator category probability  threshold
#> 1        1 time_with_friends       no 0.829949372  1.5852685
#> 2        2 time_with_friends       no 0.975758476  3.6951481
#> 3        1 time_with_friends      yes 0.170050628         NA
#> 4        2 time_with_friends      yes 0.024241524         NA
#> 5        1   on_social_media       no 0.992389457  4.8705811
#> 6        2   on_social_media       no 0.372739579 -0.5204818
#> 7        1   on_social_media      yes 0.007610543         NA
#> 8        2   on_social_media      yes 0.627260421         NA
#> 9        1    tv_video_games       no 0.866790096  1.8728707
#> 10       2    tv_video_games       no 0.736932548  1.0300859
#> 11       1    tv_video_games      yes 0.133209904         NA
#> 12       2    tv_video_games      yes 0.263067452         NA
#>    probability_standard_error
#> 1                  0.04305267
#> 2                  0.01675757
#> 3                  0.04305267
#> 4                  0.01675757
#> 5                  0.12189263
#> 6                  0.06817299
#> 7                  0.12189263
#> 8                  0.06817299
#> 9                  0.04410404
#> 10                 0.03726423
#> 11                 0.04410404
#> 12                 0.03726423
get_results(fit, "profile_probabilities")
#>   group_class profile  probability group_class_probability
#> 1           1       1 7.775664e-01               0.5460114
#> 2           1       2 2.224336e-01               0.5460114
#> 3           2       1 2.380180e-07               0.4539886
#> 4           2       2 9.999998e-01               0.4539886
```
