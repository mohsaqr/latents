# Bivariate residuals of the indicators

Each pair of indicators within each class, checked for association the
model leaves unexplained: a residual correlation with a normal test for
two Gaussian indicators, a chi-square of observed against expected
counts for two categorical ones, and a descriptive correlation without a
test for a mixed pair or one with a count indicator, with
Benjamini-Hochberg adjusted p-values, for the families that define them.
For a fit, the table is the one
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
computes with its defaults; for a grid, the tables of every candidate
are stacked with a first column `candidate` (a candidate whose
estimation failed has no row); for a diagnostics result, its own table.

## Usage

``` r
get_residuals(x)
```

## Arguments

- x:

  A fitted model of any family, an enumeration grid, or the result of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).

## Value

A data frame, one row per class and indicator pair (and candidate):
`profile, indicator_1, indicator_2, kind, observed, expected, residual, effective_n, statistic, df, p_value, p_adjusted`.
A family without bivariate residuals (regression, growth, additive and
cross-level fits, a general transition fit with occasion-dependent
structure, or a fit with a noise component) raises
`latents_bad_argument`; a grid none of whose candidates has them does
too. Any other object raises `latents_not_a_fit`.

## See also

[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
for the test;
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md).

## Examples

``` r
fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
           n_classes = 2, n_starts = 2, seed = 1)
get_residuals(fit)
#>     profile       indicator_1     indicator_2        kind     observed
#> 1 profile_2 time_with_friends          sports categorical 7.373151e-02
#> 2 profile_1 time_with_friends          sports categorical 1.695573e-02
#> 3 profile_1   on_social_media          sports categorical 7.743261e-04
#> 4 profile_1 time_with_friends on_social_media categorical 4.846465e-04
#> 5 profile_2   on_social_media          sports categorical 4.268376e-05
#> 6 profile_2 time_with_friends on_social_media categorical 2.674778e-05
#>       expected     residual effective_n    statistic df    p_value p_adjusted
#> 1 0.000000e+00 7.373151e-02    740.9499 4.028016e+00  1 0.04475055  0.2685033
#> 2 3.454893e-17 1.695573e-02   1841.0501 5.292958e-01  1 0.46690306  0.9994096
#> 3 5.529567e-17 7.743261e-04   1841.0501 1.115166e-03  1 0.97336032  0.9994096
#> 4 1.349080e-16 4.846465e-04   1841.0501 4.427385e-04  1 0.98321267  0.9994096
#> 5 1.577941e-16 4.268376e-05    740.9499 1.366123e-06  1 0.99906742  0.9994096
#> 6 7.263467e-19 2.674778e-05    740.9499 5.476005e-07  1 0.99940957  0.9994096
```
