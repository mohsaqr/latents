# Extract the number of independent groups

Extract the number of independent groups

## Usage

``` r
# S3 method for class 'multilpa'
nobs(object, ...)
```

## Arguments

- object:

  An `multilpa` model.

- ...:

  Reserved for compatibility with
  [`nobs()`](https://rdrr.io/r/stats/nobs.html).

## Value

A single integer: the number of observed groups, which are the
independent units of the two-level likelihood. For the individual count
alongside every other sample-size-dependent quantity, call
`get_results(x, "information_criteria")`.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
nobs(fit)
#> [1] 12
```
