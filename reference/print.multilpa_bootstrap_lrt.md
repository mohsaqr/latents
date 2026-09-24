# Print a parametric bootstrap likelihood-ratio comparison

Print a parametric bootstrap likelihood-ratio comparison

## Usage

``` r
# S3 method for class 'multilpa_bootstrap_lrt'
print(x, ...)
```

## Arguments

- x:

  An `multilpa_bootstrap_lrt` result.

- ...:

  Reserved for compatibility with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The input, invisibly. Called for the side effect of printing the two
models compared, the observed statistic, the p-value with its Monte
Carlo standard error and the number of valid replicates, and the blocks
held fixed where there are any.

## Examples

``` r
set.seed(1)
example_data <- data.frame(
  school = rep(seq_len(10), each = 10),
  score = rnorm(100, rep(c(-2, 2), each = 50))
)
smaller <- multilpa(example_data, "score", "school", n_profiles = 1,
                    n_group_classes = 1, n_starts = 2, seed = 1)
larger <- multilpa(example_data, "score", "school", n_profiles = 2,
                   n_group_classes = 1, n_starts = 2, seed = 1)
# `iter` is small so the example runs quickly; use many more for inference.
comparison <- bootstrap_lrt(smaller, larger, iter = 9, n_starts = 1,
                            max_iter = 2000, tol = 1e-6, seed = 1)
print(comparison)
#> Parametric bootstrap likelihood-ratio comparison
#> Null: 1 profiles, 1 group classes; alternative: 2 profiles, 1 group classes
#> Observed statistic: 44.198110
#> p-value: 0.1 (Monte Carlo SE 0.0949) from 9 of 9 valid replicates
```
