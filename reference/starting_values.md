# Get the estimates of a fitted profile or class model as starting values for another fit

`starting_values()` returns the parameter estimates of a fitted profile
or class model in the form
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
accepts as `start`, so that a second fit begins from the first one's
solution. The measurement parameters can be taken alone, for a refit
that estimates the class probabilities anew. Starting values are where
the estimation begins, not constraints: to hold parameters fixed, use
the `fixed` argument of the fitting function.

## Usage

``` r
starting_values(
  x,
  covariance = c("auto", "drop", "keep"),
  what = c("all", "measurement")
)
```

## Arguments

- x:

  The source of the values, one of:

  - a fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    without membership covariates;

  - a list with `profile_probabilities`, `group_probabilities` and the
    measurement blocks the model uses: `means` with `variances` or
    `covariances` for continuous indicators, `response_probabilities`
    for categorical ones (a model with only categorical indicators needs
    no `means`). Other elements are ignored.

  A fit with membership covariates and a model of another family carry
  no such blocks and raise `latents_bad_start`, as does a transition fit
  from [`lta()`](https://pak.dynasite.org/latents/reference/lta.md),
  which takes no `start`.

- covariance:

  How the within-profile covariances are returned:

  - `"auto"` (default): `covariances` when `x` has them (a fit with a
    full covariance structure), otherwise `variances`;

  - `"drop"`: always `variances`, the diagonal form;

  - `"keep"`: `covariances`, raising `latents_bad_start` when `x` has
    none.

- what:

  Which blocks to return:

  - `"all"` (default): the measurement blocks and the mixing blocks
    (`profile_probabilities` and `group_probabilities`);

  - `"measurement"`: the measurement blocks only. The mixing blocks are
    sized for the model that produced them, so leave them out for a
    refit with another number of group classes, together with `fixed` in
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

## Value

An object of class `multilpa_start`: the list that
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
accepts as `start`, which can be passed on as it is. Its elements:

- `means`: a matrix, profiles by continuous indicators;

- `variances`: a matrix of the same shape; or, when covariances are
  returned, `covariances` in its place: an array, indicators by
  indicators by profiles;

- `profile_probabilities`: a matrix, group classes by profiles, and
  `group_probabilities`: a vector with one value per group class. Both
  are left out under `what = "measurement"`;

- `response_probabilities`: for a model with categorical indicators, a
  list with one matrix per indicator, profiles by categories.

The Gaussian blocks are unnamed and matched by position.
`response_probabilities` keeps the indicator and category names, and
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
matches its blocks by them.

Printing the object summarises the blocks.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the means and variances as a table, one row per profile and indicator,
and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
the tables `"profiles"`, `"responses"`, `"profile_probabilities"` and
`"covariances"`.

Conditions raised:

- `latents_bad_start`: `x` lacks a block that is needed, is a transition
  fit, or has no covariances under `covariance = "keep"`.
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  raises the same class for a `start` whose categorical indicators or
  categories do not match the model;

- `latents_unsupported_noise`: `x` is a fit with a noise component;

- `latents_bad_argument`: `x` is not a list.

## Details

This is a utility for controlled refitting. `starting_values()` checks
that the blocks it needs are present; the fitting function checks their
dimensions and values.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the values as tidy tables;
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
for two-stage estimation with fixed measurement parameters.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
starting_values(fit)
#> Starting values: 2 profile(s), 2 continuous indicator(s)
#> Residual covariance: diagonal
#> Categorical indicators: 0
#> Mixing blocks: 1 group class(es)
#> Pass this to multilpa(start = ) as it is.
as.data.frame(starting_values(fit))
#>   profile indicator        mean  variance standard_deviation
#> 1       1         1  0.31636435 0.7483133          0.8650510
#> 2       1         2 -0.05291897 0.9149761          0.9565438
#> 3       2         1 -1.11359384 0.1018719          0.3191738
#> 4       2         2  0.79093314 0.2857893          0.5345926
refit <- multilpa(example_data, c("score_a", "score_b"), "school",
                  n_profiles = 2, n_group_classes = 1, n_starts = 1,
                  start = starting_values(fit))
logLik(refit)
#> 'log Lik.' -323.0197 (df=9)
```
