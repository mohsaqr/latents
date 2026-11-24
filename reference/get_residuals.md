# Get the bivariate residuals: the association between pairs of indicators that each class leaves unexplained

`get_residuals()` returns the bivariate residuals, a check of local
independence: for each pair of indicators within each class, the
association that the model leaves unexplained.

## Usage

``` r
get_residuals(x)
```

## Arguments

- x:

  One of:

  - a profile, class or transition fit as described above. The table is
    the one
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    computes with its defaults, which assesses each class separately;

  - an enumeration grid of such fits. The tables of the fitted
    candidates are stacked, with a first column `candidate`;

  - the result of
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).
    The table is the one the result holds, so
    `diagnostics(fit, by = "overall")` gives the residuals pooled over
    the classes.

## Value

A data frame with one row per class and pair of indicators (and per
candidate for a grid) and the columns:

- `profile`: the class, as `profile_1`, `profile_2`, ...; `overall` for
  residuals pooled over the classes;

- `indicator_1`, `indicator_2`: the pair;

- `kind`: the kind of pair, as listed above;

- `observed`, `expected`: the observed association and the one the model
  implies;

- `residual`: the residual association;

- `effective_n`: the effective number of units behind it;

- `statistic`, `df`, `p_value`: the test, `NA` where there is none;

- `p_adjusted`: the p-value adjusted across the pairs.

Errors: `latents_bad_argument` when the fit is of a family without
bivariate residuals (a mixture regression, a growth mixture, an additive
or cross-level fit, or a general transition fit of class
`multilpa_lta`), and when no candidate of a grid has them; the message
names the tables the object does have. `latents_unsupported_noise` when
the fit has a noise component. `latents_not_a_fit` when `x` is not a
fit, a grid or a diagnostics result.

## Details

The statistic depends on the two indicators (column `kind`):

- `"gaussian"`, two continuous indicators: the residual correlation,
  with a normal test;

- `"categorical"`, two categorical or ordinal indicators: a chi-square
  test of the observed two-way table against the one the model expects;

- `"mixed"` or `"count"`, a continuous with a categorical indicator or a
  pair with a count indicator: a correlation, without a test.

The p-values are adjusted across the pairs by the Benjamini-Hochberg
method. A fit with sampling weights reports the associations without
tests.
[`?get_results`](https://pak.dynasite.org/latents/reference/get_results.md)
describes the statistics in full.

Only profile, class and transition fits have these residuals: fits of
[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
[`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
and
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
with or without covariates, and transition fits of
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) of class
`multilpa_transitions`.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the statistics;
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
and
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md)
for the check that reads this table.

## Examples

``` r
fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
           n_classes = 2, n_starts = 2, seed = 1)
get_residuals(fit)
#>     profile       indicator_1     indicator_2        kind     observed
#> 1 profile_2 time_with_friends          sports categorical 7.373151e-02
#> 2 profile_1 time_with_friends          sports categorical 1.695572e-02
#> 3 profile_1   on_social_media          sports categorical 7.743269e-04
#> 4 profile_1 time_with_friends on_social_media categorical 4.846471e-04
#> 5 profile_2   on_social_media          sports categorical 4.268375e-05
#> 6 profile_2 time_with_friends on_social_media categorical 2.674778e-05
#>       expected     residual effective_n    statistic df    p_value p_adjusted
#> 1 1.774695e-16 7.373151e-02    740.9498 4.028016e+00  1 0.04475055  0.2685033
#> 2 1.414172e-16 1.695572e-02   1841.0502 5.292958e-01  1 0.46690308  0.9994096
#> 3 2.210547e-18 7.743269e-04   1841.0502 1.115168e-03  1 0.97336029  0.9994096
#> 4 7.739697e-17 4.846471e-04   1841.0502 4.427395e-04  1 0.98321265  0.9994096
#> 5 4.163988e-17 4.268375e-05    740.9498 1.366122e-06  1 0.99906742  0.9994096
#> 6 5.955492e-17 2.674778e-05    740.9498 5.476002e-07  1 0.99940957  0.9994096
```
