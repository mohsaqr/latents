# Print a set of starting values

Print a set of starting values

## Usage

``` r
# S3 method for class 'multilpa_start'
print(x, ...)
```

## Arguments

- x:

  A `multilpa_start` object from
  [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md).

- ...:

  Reserved for compatibility with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The input, invisibly.

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
print(starting_values(fit))
#> Starting values: 2 profile(s), 2 continuous indicator(s)
#> Residual covariance: diagonal
#> Categorical indicators: 0
#> Mixing blocks: 1 group class(es)
#> Pass this to multilpa(start = ) as it is.
```
