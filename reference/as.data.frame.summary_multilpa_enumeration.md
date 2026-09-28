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
#> 1          1               1   EEI     -10007.463           10 20034.93
#> 2          2               1   EEI      -8654.537           16 17341.07
#> 3          1               1   VVI     -10007.463           10 20034.93
#> 4          2               1   VVI      -8615.233           21 17272.47
#> 5          1               1   EEE      -8710.142           20 17460.28
#> 6          2               1   EEE      -8536.422           26 17124.84
#> 7          1               1   VVV      -8710.142           20 17460.28
#> 8          2               1   VVV      -8489.610           41 17061.22
#>        kic bic_groups bic_individual sabic_groups sabic_individual caic_groups
#> 1 20047.93   20061.56       20087.52     20029.97         20055.76    20071.56
#> 2 17360.07   17383.69       17425.23     17333.14         17374.40    17399.69
#> 3 20047.93   20061.56       20087.52     20029.97         20055.76    20071.56
#> 4 17296.47   17328.40       17382.92     17262.05         17316.21    17349.40
#> 5 17483.28   17513.55       17565.48     17450.37         17501.95    17533.55
#> 6 17153.84   17194.09       17261.60     17111.95         17179.01    17220.09
#> 7 17483.28   17513.55       17565.48     17450.37         17501.95    17533.55
#> 8 17105.22   17170.42       17276.87     17040.89         17146.63    17211.42
#>   caic_individual awe_groups awe_individual icl_groups icl_individual
#> 1        20097.52   20138.19       20190.12   20061.56       20087.52
#> 2        17441.23   17506.30       17755.41   17383.69       17591.25
#> 3        20097.52   20138.19       20190.12   20061.56       20087.52
#> 4        17403.92   17489.33       17752.78   17328.40       17537.32
#> 5        17585.48   17666.82       17770.68   17513.55       17565.48
#> 6        17287.60   17393.34       17793.05   17194.09       17526.29
#> 7        17585.48   17666.82       17770.68   17513.55       17565.48
#> 8        17317.87   17484.62       17944.67   17170.42       17524.02
#>   clc_groups clc_individual profile_entropy group_entropy converged boundary
#> 1   20014.93       20014.93              NA            NA      TRUE    FALSE
#> 2   17309.07       17475.10       0.9157810            NA      TRUE    FALSE
#> 3   20014.93       20014.93              NA            NA      TRUE    FALSE
#> 4   17230.47       17384.87       0.9216765            NA      TRUE    FALSE
#> 5   17420.28       17420.28              NA            NA      TRUE    FALSE
#> 6   17072.84       17337.54       0.8657262            NA      TRUE    FALSE
#> 7   17420.28       17420.28              NA            NA      TRUE    FALSE
#> 8   16979.22       17226.37       0.8746271            NA      TRUE    FALSE
#>   n_best_replicated warnings error
#> 1                 2           <NA>
#> 2                 2           <NA>
#> 3                 2           <NA>
#> 4                 2           <NA>
#> 5                 2           <NA>
#> 6                 2           <NA>
#> 7                 2           <NA>
#> 8                 2           <NA>
get_results(summary(candidates), what = "criteria")
#>    criterion  convention n_profiles n_group_classes model    value
#> 1        aic        <NA>          2               1   VVV 17061.22
#> 2        kic        <NA>          2               1   VVV 17105.22
#> 3        bic      groups          2               1   VVV 17170.42
#> 4        bic individuals          2               1   EEE 17261.60
#> 5      sabic      groups          2               1   VVV 17040.89
#> 6      sabic individuals          2               1   VVV 17146.63
#> 7       caic      groups          2               1   VVV 17211.42
#> 8       caic individuals          2               1   EEE 17287.60
#> 9        awe      groups          2               1   EEE 17393.34
#> 10       awe individuals          2               1   VVI 17752.78
#> 11       icl      groups          2               1   VVV 17170.42
#> 12       icl individuals          2               1   VVV 17524.02
#> 13       clc      groups          2               1   VVV 16979.22
#> 14       clc individuals          2               1   VVV 17226.37
```
