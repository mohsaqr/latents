# Coerce a class-enumeration summary to its primary table

Plain coercion, as the base generic means it: one object, one data
frame. A summary carries every table the object it describes can
produce, and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
names them.

## Usage

``` r
# S3 method for class 'summary_multilpa_enumeration'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `summary_multilpa_enumeration`.

- row.names:

  Passed to [`data.frame()`](https://rdrr.io/r/base/data.frame.html);
  `NULL` gives default row names.

- optional:

  Ignored, present for generic compatibility.

- ...:

  Must be empty. An argument here raises `latents_bad_argument` naming
  it, rather than being dropped.

## Value

A base `data.frame`: one row per candidate model in the grid.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for every other table this summary holds.

## Examples

``` r
candidates <- enumerate_classes(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 1:2, n_group_classes = 1, n_starts = 2,
  seed = 1
)
as.data.frame(summary(candidates))
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic
#> 1          1               1   VVI     -10007.463           10 20034.93
#> 2          2               1   VVI      -8615.233           21 17272.47
#>        kic bic_groups bic_individual sabic_groups sabic_individual caic_groups
#> 1 20047.93   20061.56       20087.52     20029.97         20055.76    20071.56
#> 2 17296.47   17328.40       17382.92     17262.05         17316.21    17349.40
#>   caic_individual awe_groups awe_individual icl_groups icl_individual
#> 1        20097.52   20138.19       20190.12   20061.56       20087.52
#> 2        17403.92   17489.33       17752.78   17328.40       17537.33
#>   clc_groups clc_individual profile_entropy group_entropy converged boundary
#> 1   20014.93       20014.93              NA            NA      TRUE    FALSE
#> 2   17230.47       17384.87       0.9216738            NA      TRUE    FALSE
#>   n_best_replicated warnings error
#> 1                 2           <NA>
#> 2                 2           <NA>
get_results(summary(candidates), what = "criteria")
#>    criterion  convention n_profiles n_group_classes model    value
#> 1        aic        <NA>          2               1   VVI 17272.47
#> 2        kic        <NA>          2               1   VVI 17296.47
#> 3        bic      groups          2               1   VVI 17328.40
#> 4        bic individuals          2               1   VVI 17382.92
#> 5      sabic      groups          2               1   VVI 17262.05
#> 6      sabic individuals          2               1   VVI 17316.21
#> 7       caic      groups          2               1   VVI 17349.40
#> 8       caic individuals          2               1   VVI 17403.92
#> 9        awe      groups          2               1   VVI 17489.33
#> 10       awe individuals          2               1   VVI 17752.78
#> 11       icl      groups          2               1   VVI 17328.40
#> 12       icl individuals          2               1   VVI 17537.33
#> 13       clc      groups          2               1   VVI 17230.47
#> 14       clc individuals          2               1   VVI 17384.87
```
