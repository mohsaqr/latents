# Every classification diagnostic, in one call

[`summary()`](https://rdrr.io/r/base/summary.html) reports what a model
estimated. This reports whether to believe it: how sharply the posterior
separates the classes, how big each class actually is, how confidently
each unit was assigned, and whether the within-profile independence the
model assumes survives contact with the data.

## Usage

``` r
diagnostics(x, data = NULL, plots = FALSE, by = c("profile", "overall"), ...)

# S3 method for class 'multilpa_diagnostics'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'multilpa_diagnostics'
print(x, ...)

# S3 method for class 'multilpa_diagnostics'
plot(x, ...)
```

## Arguments

- x:

  A fitted model of this package.

- data:

  Optional. The data the model was fitted to. A fit carries the columns
  it was built from, so this is only needed to override them.

- plots:

  `TRUE` also draws the classification plots, as a side effect, before
  returning. Equivalent to calling
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the
  result.

- by:

  `"profile"`, the default, assesses each profile's bivariate residuals
  separately; `"overall"` pools them. It is the one argument of a
  gathered table this function forwards, because it is the one that
  changes what a gathered table means rather than which fit it is taken
  from. The `level` the classification tables use is not an argument
  here: it follows from whether the fit has discrete group classes.

- ...:

  For `diagnostics()`, nothing further is accepted. An argument this
  function cannot forward raises an error of class
  `latents_bad_argument` naming it, rather than being dropped on the way
  to a table that then means something other than what was asked for.
  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), style
  overrides, as in
  [`plot.multilpa()`](https://pak.dynasite.org/latents/reference/plot.multilpa.md).

- row.names, optional:

  Passed to the base generic; `row.names` is applied to the returned
  table.

## Value

An object of class `multilpa_diagnostics`. Read its tables with
`get_results(result, what = )`, which offers `"entropy"`,
`"classification"`, `"average_posteriors"`, `"residuals"` and `"all"`,
and never with `$`. A model family that has no bivariate residuals
leaves that table out of `"all"`, and asking for it by name raises
`latents_no_group_classes`.

[`print()`](https://rdrr.io/r/base/print.html) returns the object
invisibly, having printed one line per diagnostic: relative entropy,
smallest class, lowest average posterior and largest residual, at each
level the fit has.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns the
object invisibly, having drawn the case-level entropy and posterior
panels. [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
returns the entropy table, the primary one.

## Details

It gathers the four classification tables of
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
— `"entropy"`, `"classification"`, `"average_posteriors"` and
`"residuals"` — and prints one line of reading per diagnostic rather
than four differently shaped tables. The tables themselves come back
from
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
on the result, or on the fit.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for these tables and every other one,
[`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
for the before-the-fit counterpart, and
[`summary()`](https://rdrr.io/r/base/summary.html) for what the model
estimated rather than whether to trust it.

## Examples

``` r
fit <- multilpa(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
  seed = 1
)
quality <- diagnostics(fit)
quality
#> Classification quality: 2 profiles, 2 group classes
#> 
#>   Relative entropy     individuals 0.938    groups 0.940
#>   Smallest class       individuals 587.9 (41.3%)    groups 34.5 (32.5%)
#>   Lowest avg posterior individuals 0.982    groups 0.969
#>   Largest residual     0.234  (lectures, attendance; profile_1)
#> 
#> Tables: get_results(x, what = "entropy" | "classification" | 
#>         "average_posteriors" | "residuals" | "all"). plot(x) draws them.
get_results(quality, what = "classification")
#>         level class n_modal proportion_modal estimated_n estimated_proportion
#> 1 individuals     1     836        0.5879044   834.08182            0.5865554
#> 2 individuals     2     586        0.4120956   587.91818            0.4134446
#> 3      groups     1      35        0.3301887    34.47767            0.3252610
#> 4      groups     2      71        0.6698113    71.52233            0.6747390
#>   average_posterior odds_correct_classification
#> 1         0.9849254                    46.05384
#> 2         0.9817676                    76.39357
#> 3         0.9689220                    64.67542
#> 4         0.9920366                    60.05195
get_results(diagnostics(fit, by = "overall"), what = "residuals")
#>    profile indicator_1 indicator_2     kind     observed expected     residual
#> 1  overall  forum_read  attendance gaussian  0.215876827        0  0.215876827
#> 2  overall    lectures  attendance gaussian  0.192335674        0  0.192335674
#> 3  overall      browse  attendance gaussian  0.180768516        0  0.180768516
#> 4  overall  forum_post  attendance gaussian  0.172386082        0  0.172386082
#> 5  overall      browse  forum_read gaussian  0.062656428        0  0.062656428
#> 6  overall      browse  forum_post gaussian -0.025112589        0 -0.025112589
#> 7  overall      browse    lectures gaussian -0.021262512        0 -0.021262512
#> 8  overall  forum_read  forum_post gaussian  0.019856866        0  0.019856866
#> 9  overall    lectures  forum_post gaussian  0.017086939        0  0.017086939
#> 10 overall    lectures  forum_read gaussian  0.001918653        0  0.001918653
#>    effective_n  statistic df      p_value   p_adjusted
#> 1         1422  8.2619762 NA 1.432803e-16 1.432803e-16
#> 2         1422  7.3365887 NA 2.191066e-13 2.191066e-13
#> 3         1422  6.8851414 NA 5.773015e-12 5.773015e-12
#> 4         1422  6.5592139 NA 5.409216e-11 5.409216e-11
#> 5         1422  2.3633395 NA 1.811107e-02 1.811107e-02
#> 6         1422 -0.9461805 NA 3.440565e-01 3.440565e-01
#> 7         1422 -0.8010714 NA 4.230903e-01 4.230903e-01
#> 8         1422  0.7480989 NA 4.544005e-01 4.544005e-01
#> 9         1422  0.6437211 NA 5.197563e-01 5.197563e-01
#> 10        1422  0.0722750 NA 9.423831e-01 9.423831e-01
```
