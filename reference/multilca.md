# Fit a two-level latent class model: classes of categorical responses and classes of the groups they are nested in

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
  and factor columns are accepted, and every distinct value of a column
  is one category; factor levels retain their declared order and other
  types are sorted.

- id:

  Name of the group identifier column, or `NULL` for a single-level
  model ([`lca()`](https://pak.dynasite.org/latents/reference/lca.md)
  fits that directly).

- n_profiles:

  Number of observation-level latent classes (profiles).

- n_group_classes:

  Number of group-level latent classes. Defaults to 2, or to 1 when
  `id = NULL`. With 1, the groups do not differ in their class
  proportions.

- ...:

  Further arguments, passed to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  which documents all of them. The ones used most often:

  - `missing`: `"error"` (default) stops when an indicator has a missing
    value; `"fiml"` keeps every row and uses the responses it has.

  - `ordinal`, `count`: indicators in `vars` modelled as ordered
    categories or as counts; the remaining indicators are categorical.

  - `n_starts` (default 10), `seed`, `max_iter` and `tol`: the number of
    starts, the seed that makes them reproducible, and the limits of
    each start.

  - `min_probability`: the lower bound on a response probability
    (default `1e-10`).

  - `profile_covariates`, `group_covariates`: covariates that predict
    class and group-class membership.

  - `time`: the occasion column, used to order each group's rows in the
    sequence tables and
    [`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md).

  `categorical` is accepted but has no further effect; its names must be
  in `vars`.

## Value

A fitted model of class `multilpa`, or `multilpa_covariates` when
membership covariates are specified, as returned by
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
The classes are numbered 1, 2, and so on in the `profile` column of its
tables, and the group classes in `group_class`.
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
returns the tables by name, among them:

- `"responses"` (also `as.data.frame(fit)`): one row per class, item and
  category, with the response probability;

- `"profile_probabilities"`: the proportion of each class within each
  group class, and the proportion of each group class;

- `"assignments"`: the data with each row's most likely class, its
  group's most likely group class and the posterior probabilities;

- `"group_posteriors"`: the posterior probability of each group class
  for each group;

- `"classification"` and `"entropy"`: how clearly rows and groups are
  classified;

- `"information_criteria"` and `"starts"`.

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

Errors: `latents_bad_argument` when `id` is left out or an argument has
a value that is not accepted; `latents_bad_categorical` when
`categorical` names a column not in `vars`; `latents_bad_data` when a
column is missing from `data` or an indicator has missing values under
`missing = "error"`; `latents_unidentified` when more than one group
class is requested with one class or with one row per group. Warnings:
`latents_unconverged` when the reported start did not converge and
`latents_small_classes` when a class or group class is nearly empty.
Otherwise the conditions of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);
[latents-conditions](https://pak.dynasite.org/latents/reference/latents-conditions.md)
lists them all.

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for continuous and mixed indicators, and
[`vignette("lca", package = "latents")`](https://pak.dynasite.org/latents/articles/lca.md).

## Examples

``` r
activities <- c("time_with_friends", "on_social_media", "tv_video_games")
fit <- multilca(subset(student_esm, day <= 1), vars = activities,
                id = "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 2, seed = 1)
get_results(fit, "responses")
#>    profile         indicator category probability  threshold
#> 1        1 time_with_friends       no 0.975634724  3.6899293
#> 2        2 time_with_friends       no 0.829730459  1.5837182
#> 3        1 time_with_friends      yes 0.024365276         NA
#> 4        2 time_with_friends      yes 0.170269541         NA
#> 5        1   on_social_media       no 0.372974753 -0.5194761
#> 6        2   on_social_media       no 0.993730396  5.0657527
#> 7        1   on_social_media      yes 0.627025247         NA
#> 8        2   on_social_media      yes 0.006269604         NA
#> 9        1    tv_video_games       no 0.737018500  1.0305293
#> 10       2    tv_video_games       no 0.867019325  1.8748575
#> 11       1    tv_video_games      yes 0.262981500         NA
#> 12       2    tv_video_games      yes 0.132980675         NA
#>    probability_standard_error
#> 1                  0.01717744
#> 2                  0.04359206
#> 3                  0.01717744
#> 4                  0.04359206
#> 5                  0.06856696
#> 6                  0.12850275
#> 7                  0.06856696
#> 8                  0.12850275
#> 9                  0.03725565
#> 10                 0.04461679
#> 11                 0.03725565
#> 12                 0.04461679
get_results(fit, "profile_probabilities")
#>   group_class profile  probability group_class_probability
#> 1           1       1 2.240147e-01                 0.54559
#> 2           1       2 7.759853e-01                 0.54559
#> 3           2       1 1.000000e+00                 0.45441
#> 4           2       2 1.688848e-09                 0.45441
```
