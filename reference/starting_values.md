# Extract starting values for a latent profile fit

Returns a `multilpa_start` object, the starting-value list
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
accepts as `start`, taken from a fitted model or from any list that
carries the parameter blocks. The blocks are validated when the list is
built.

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

  A fitted `multilpa` model, or a list carrying `profile_probabilities`,
  `group_probabilities`, and the measurement blocks the model uses:
  `means` with `variances` or `covariances` for continuous indicators,
  `response_probabilities` for categorical ones. An all-categorical
  model needs no Gaussian block. Other elements are ignored. A missing
  block raises `latents_bad_start`; a fit with a noise component raises
  `latents_unsupported_noise`.

- covariance:

  `"auto"` (the default) keeps `covariances` when the object carries
  them, `"drop"` always returns the diagonal parameterization, and
  `"keep"` requires `covariances` and raises `latents_bad_start` when
  they are absent.

- what:

  `"all"` (the default) returns the measurement and mixing blocks.
  `"measurement"` omits `profile_probabilities` and
  `group_probabilities`, which are sized for the model that produced
  them; use it with `fixed` in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  for a staged fit that changes the number of group classes.

## Value

An object of class `multilpa_start`: the list
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
accepts as `start`, with elements `means`, `variances`,
`profile_probabilities` and `group_probabilities`. When covariances are
returned, `covariances` replaces `variances`; when the object carries a
categorical measurement model, `response_probabilities` is added.
`what = "measurement"` drops `profile_probabilities` and
`group_probabilities`. The Gaussian blocks are unnamed and indexed by
position; [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
returns them as a labelled table. `response_probabilities` keeps each
indicator's name and each block's category names, and
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
matches the blocks by these labels, so a fit whose categorical
indicators or categories do not match raises `latents_bad_start`. The
object can be passed directly as `start`.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the values as tidy tables.

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
#> 1       1         1  0.31636979 0.7483074          0.8650476
#> 2       1         2 -0.05291982 0.9149785          0.9565451
#> 3       2         1 -1.11359585 0.1018727          0.3191750
#> 4       2         2  0.79091589 0.2857973          0.5346002
refit <- multilpa(example_data, c("score_a", "score_b"), "school",
                    n_profiles = 2, n_group_classes = 1, n_starts = 1,
                    start = starting_values(fit))
logLik(refit)
#> 'log Lik.' -323.0197 (df=9)
```
