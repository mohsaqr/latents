# Extract the multilevel model log likelihood

Extract the multilevel model log likelihood

## Usage

``` r
# S3 method for class 'multilpa'
logLik(object, ...)
```

## Arguments

- object:

  An `multilpa` model.

- ...:

  Reserved for compatibility with
  [`logLik()`](https://rdrr.io/r/stats/logLik.html).

## Value

A `logLik` object with parameter count `df` and the number of observed
groups as `nobs`. Thus
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) uses group-count BIC.
For the individual-count alternative, and every other criterion, call
`get_results(x, "information_criteria")`, which reports both conventions
side by side.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
logLik(fit)
#> 'log Lik.' -323.0197 (df=9)
AIC(fit)
#> [1] 664.0394
BIC(fit)
#> [1] 668.4036
```
