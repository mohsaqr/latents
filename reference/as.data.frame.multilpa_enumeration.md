# Coerce a class-enumeration grid to its primary table

Plain coercion, as the base generic means it: one object, one data
frame. The other tables are named rather than positional, so they belong
to
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
which takes `what` and refuses a name this object has not.

## Usage

``` r
# S3 method for class 'multilpa_enumeration'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `multilpa_enumeration`.

- row.names:

  Passed to [`data.frame()`](https://rdrr.io/r/base/data.frame.html);
  `NULL` gives default row names.

- optional:

  Ignored, present for generic compatibility.

- ...:

  Must be empty. An argument here raises `latents_bad_argument` naming
  it, rather than being dropped, because `what =` used to live on this
  generic and silently returning the primary table instead of the one
  that was asked for is the one outcome worth refusing.

## Value

A base `data.frame`: one row per candidate model in the grid.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for every other table this object holds.

## Examples

``` r
candidates <- enumerate_classes(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 1:2, n_group_classes = 1, n_starts = 2,
  seed = 1
)
as.data.frame(candidates)
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
```
