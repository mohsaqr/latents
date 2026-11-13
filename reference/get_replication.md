# How many starts reached the best log likelihood

The number of starts, how many were screened, converged and failed, and
how many replicated the best log likelihood within tolerance. For a fit,
the table is the one
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
computes with its defaults; for a grid, the tables of every candidate
are stacked with a first column `candidate` (a candidate whose
estimation failed has no row); for a diagnostics result, its own table.

## Usage

``` r
get_replication(x)
```

## Arguments

- x:

  A fitted model of any family, an enumeration grid, or the result of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).

## Value

A data frame, one row per fit (and candidate):
`n_starts, n_screened, n_converged, n_failed, n_best_replicated, share, tolerance, best_start, log_likelihood`.
Any other object raises `latents_not_a_fit`.

## See also

[`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md)
draws it for a grid;
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
get_replication(fit)
#>   n_starts n_screened n_converged n_failed n_best_replicated share tolerance
#> 1        2          0           2        0                 2     1  0.000544
#>   best_start log_likelihood
#> 1          1      -1147.806
```
