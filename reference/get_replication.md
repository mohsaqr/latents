# Get how many random starts converged, failed and reached the best log likelihood

`get_replication()` returns how many of a model's random starts ran,
converged and failed, and how many arrived at the best log likelihood. A
best solution that only one start reached may not be the highest one.

## Usage

``` r
get_replication(x)
```

## Arguments

- x:

  One of:

  - a fitted model of any family: from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).
    The table is the one
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    computes with its defaults;

  - an enumeration grid of any kind, from
    [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
    [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
    [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
    or
    [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).
    The tables of the fitted candidates are stacked, with a first column
    `candidate` that names the candidate of each row. A candidate whose
    estimation failed has no rows;

  - the result of
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).
    The table is the one the result holds.

## Value

A data frame with one row per fit (one row per candidate for a grid) and
the columns:

- `n_starts`: the number of starts that were run to the end;

- `n_screened`: the number of starts that were only screened in a short
  run; 0 when the fit screened none;

- `n_converged`: the number of starts that converged;

- `n_failed`: the number of starts that failed or gave no log
  likelihood;

- `n_best_replicated`: the number of starts whose log likelihood is
  within `tolerance` of the best;

- `share`: `n_best_replicated` divided by `n_starts`;

- `tolerance`: the difference in log likelihood up to which two starts
  count as the same solution;

- `best_start`: the number of the start that gave the fit;

- `log_likelihood`: the log likelihood of the fit.

A count the fit does not record is `NA`.

Errors: `latents_not_a_fit` when `x` is none of these, for example a
comparison, pooled imputations or a data frame.

## See also

[`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md)
draws the replication of every candidate of a grid;
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md);
[`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
refits a model from other seeds.

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
get_replication(fit)
#>   n_starts n_screened n_converged n_failed n_best_replicated share tolerance
#> 1        2          0           2        0                 2     1  0.000544
#>   best_start log_likelihood
#> 1          1      -1147.806

# One row per candidate of a grid
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
get_replication(grid)
#>   candidate n_starts n_screened n_converged n_failed n_best_replicated share
#> 1 profiles1        2          8           2        0                 2     1
#> 2 profiles2        2          8           2        0                 2     1
#> 3 profiles3        2          8           2        0                 2     1
#>   tolerance best_start log_likelihood
#> 1  0.000544          1      -1516.706
#> 2  0.000544          1      -1157.680
#> 3  0.000544          1      -1133.455
```
