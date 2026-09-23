# Build starting values for a multilevel latent profile fit

Returns the starting-value list
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
accepts, taken from a fitted model or from any object that already
carries the parameter blocks, such as a retained reference solution.
This exists so that callers never assemble the list by hand, and so that
starting values are validated where they are built rather than deep
inside the fitting loop.

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
  `group_probabilities`, and whichever measurement blocks the model
  uses: `means` with `variances` or `covariances` for continuous
  indicators, `response_probabilities` for categorical ones. An
  all-categorical model needs no Gaussian block. Any other elements are
  ignored.

- covariance:

  `"auto"` keeps `covariances` when the object carries them, `"drop"`
  always returns the diagonal parameterization, and `"keep"` requires
  `covariances` and fails when they are absent.

- what:

  `"all"` returns the measurement and mixing blocks. `"measurement"`
  omits `profile_probabilities` and `group_probabilities`, which are
  sized for the model that produced them and must not be carried into a
  staged fit that changes the number of group classes. Pass the result
  as `start` alongside `fixed` in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

## Value

An object of class `multilpa_start`: the list
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
accepts as `start`, with elements `means`, `variances`,
`profile_probabilities`, and `group_probabilities`, plus `covariances`
when the full-covariance parameterization is returned and
`response_probabilities` when the object carries a categorical
measurement model. `what = "measurement"` drops `profile_probabilities`
and `group_probabilities`. The Gaussian blocks are unnamed and indexed
by position, matching what
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
expects of `start`; read them with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), which
labels the positions, rather than out of the list.
`response_probabilities` is the exception: the list keeps each
indicator's name and each block keeps its category names, because those
labels are what lets
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
attach a block to the item it was estimated for rather than to whichever
item happens to sit in that position. A fit that names its categorical
indicators in another order, or that has other categories, therefore
refuses the start with `latents_bad_start` instead of reporting a
likelihood for a measurement model it never held. The list payload is
unchanged by the class, so a `multilpa_start` can be passed straight
back as `start`.

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
#> 1       1         1  0.31651781 0.7482672          0.8650244
#> 2       1         2 -0.05304516 0.9149716          0.9565415
#> 3       2         1 -1.11301874 0.1021424          0.3195973
#> 4       2         2  0.79087174 0.2858249          0.5346260
refit <- multilpa(example_data, c("score_a", "score_b"), "school",
                    n_profiles = 2, n_group_classes = 1, n_starts = 1,
                    start = starting_values(fit))
logLik(refit)
#> 'log Lik.' -323.0197 (df=9)
```
