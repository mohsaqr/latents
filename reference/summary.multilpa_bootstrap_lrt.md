# Summarise a parametric bootstrap likelihood-ratio comparison

Summarise a parametric bootstrap likelihood-ratio comparison

## Usage

``` r
# S3 method for class 'multilpa_bootstrap_lrt'
summary(object, ...)
```

## Arguments

- object:

  An `multilpa_bootstrap_lrt` result.

- ...:

  Reserved for compatibility with
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `summary_multilpa_bootstrap_lrt`, with a `print`
method and an
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) accessor.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the one-row test table by default and one row per replicate with
`what = "replicates"`.

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
comparison <- bootstrap_lrt(smaller, larger, iter = 9, n_starts = 2,
                            max_iter = 2000, tol = 1e-6, seed = 1)
summary(comparison)
#> Parametric bootstrap likelihood-ratio comparison
#> 0 replicate(s) reached a parameter boundary; 0 raised an error.
#> The reference distribution is simulated, not chi-square.
#> 
#> -- test ------------------------------------------------------------
#>  null_profiles null_group_classes alternative_profiles
#>              1                  1                    2
#>  alternative_group_classes statistic p_value monte_carlo_se iter n_valid fixed
#>                          1      44.2     0.1        0.09487    9       9  <NA>
#> 
#> -- replicates ------------------------------------------------------
#>  replicate statistic valid boundary null_replications alternative_replications
#>          1    6.3351  TRUE    FALSE                 2                        2
#>          2    1.0521  TRUE    FALSE                 2                        2
#>          3    0.9975  TRUE    FALSE                 2                        1
#>          4    0.2128  TRUE    FALSE                 2                        2
#>          5    9.2209  TRUE    FALSE                 2                        2
#>          6    1.2794  TRUE    FALSE                 2                        2
#>          7    1.7191  TRUE    FALSE                 2                        2
#>          8    3.6272  TRUE    FALSE                 2                        2
#>          9    6.5931  TRUE    FALSE                 2                        2
#>  warnings error
#>      <NA>  <NA>
#>      <NA>  <NA>
#>      <NA>  <NA>
#>      <NA>  <NA>
#>      <NA>  <NA>
#>      <NA>  <NA>
#>      <NA>  <NA>
#>      <NA>  <NA>
#>      <NA>  <NA>
#> 
#> 2 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
