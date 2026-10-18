# Classification diagnostics of a fitted model

Gathers the four classification tables of
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
(`"entropy"`, `"classification"`, `"average_posteriors"` and
`"residuals"`) for a fitted model, prints one line per diagnostic, and
by default draws the classification plots. The diagnostics are relative
entropy, estimated class sizes, average posterior probabilities of the
assigned classes, and the bivariate residuals that check the
within-profile independence the model assumes.

## Usage

``` r
diagnostics(x, data = NULL, plots = TRUE, by = c("profile", "overall"), ...)

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

  `TRUE` (the default) also draws the four classification plots:
  effective profile sizes, posterior probabilities of the assigned
  profiles, case entropy, and the average posterior probability matrix.
  Without the ggplot2 package a message says so and the tables are still
  returned. `FALSE` returns the tables only.

- by:

  `"profile"` (the default) assesses each profile's bivariate residuals
  separately; `"overall"` pools them. It is the only argument forwarded
  to the gathered tables. The `level` of the classification tables is
  set from whether the fit has discrete group classes.

- ...:

  Not used. For `diagnostics()` and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), any further
  argument raises `latents_bad_argument`; style the returned plots with
  ggplot2.

- row.names, optional:

  Part of the base generic; `row.names` is applied to the returned
  table.

## Value

An object of class `multilpa_diagnostics`. Its tables are read with
`get_results(result, what = )`, where `what` is `"entropy"`,
`"classification"`, `"average_posteriors"`, `"residuals"` or `"all"`. A
model family without bivariate residuals leaves that table out of
`"all"`, and asking for it by name raises `latents_no_group_classes`; a
fit with a noise component has none either, and raises
`latents_unsupported_noise`.

[`print()`](https://rdrr.io/r/base/print.html) prints one line per
diagnostic (relative entropy, smallest class, lowest average posterior
and largest residual, at each level the fit has) and returns the object
invisibly. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
returns the four classification plots as a `latents_plots` list of
ggplot objects, which draws all four when printed.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the entropy table.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for these and the other tables,
[`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
for the variables before fitting, and
[`summary()`](https://rdrr.io/r/base/summary.html) for the estimates.

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
#> 1 individuals     1     836        0.5879044   834.09578            0.5865652
#> 2 individuals     2     586        0.4120956   587.90422            0.4134348
#> 3      groups     1      71        0.6698113    71.51853            0.6747031
#> 4      groups     2      35        0.3301887    34.48147            0.3252969
#>   average_posterior odds_correct_classification
#> 1         0.9849349                    46.08135
#> 2         0.9817573                    76.35264
#> 3         0.9920142                    59.89166
#> 4         0.9689851                    64.80071
get_results(diagnostics(fit, by = "overall"), what = "residuals")








#>    profile indicator_1 indicator_2     kind     observed expected     residual
#> 1  overall  forum_read  attendance gaussian  0.215876301        0  0.215876301
#> 2  overall    lectures  attendance gaussian  0.192343503        0  0.192343503
#> 3  overall      browse  attendance gaussian  0.180765960        0  0.180765960
#> 4  overall  forum_post  attendance gaussian  0.172391933        0  0.172391933
#> 5  overall      browse  forum_read gaussian  0.062652530        0  0.062652530
#> 6  overall      browse  forum_post gaussian -0.025110724        0 -0.025110724
#> 7  overall      browse    lectures gaussian -0.021258479        0 -0.021258479
#> 8  overall  forum_read  forum_post gaussian  0.019861247        0  0.019861247
#> 9  overall    lectures  forum_post gaussian  0.017095209        0  0.017095209
#> 10 overall    lectures  forum_read gaussian  0.001925371        0  0.001925371
#>    effective_n   statistic df      p_value   p_adjusted
#> 1         1422  8.26195538 NA 1.433053e-16 1.433053e-16
#> 2         1422  7.33689497 NA 2.186060e-13 2.186060e-13
#> 3         1422  6.88504191 NA 5.777052e-12 5.777052e-12
#> 4         1422  6.55944109 NA 5.400982e-11 5.400982e-11
#> 5         1422  2.36319212 NA 1.811827e-02 1.811827e-02
#> 6         1422 -0.94611024 NA 3.440924e-01 3.440924e-01
#> 7         1422 -0.80091942 NA 4.231783e-01 4.231783e-01
#> 8         1422  0.74826396 NA 4.543010e-01 4.543010e-01
#> 9         1422  0.64403271 NA 5.195542e-01 5.195542e-01
#> 10        1422  0.07252808 NA 9.421817e-01 9.421817e-01
```
