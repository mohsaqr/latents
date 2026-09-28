# Everything about a fit, in one call

Prints the estimates, then the classification diagnostics, then draws
every plot the fit supports. A first look, not a substitute for the
verbs: each section is what
[`summary()`](https://rdrr.io/r/base/summary.html),
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) return, and
the tables are reachable from those rather than from here.

## Usage

``` r
report(
  x,
  data = NULL,
  plots = TRUE,
  by = c("profile", "overall"),
  rows = 10L,
  ...
)
```

## Arguments

- x:

  A fitted model of this package.

- data:

  Optional. The data the model was fitted to; a fit carries the columns
  it was built from.

- plots:

  `TRUE`, the default, draws the plots. `FALSE` prints only.

- by:

  Passed to
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md),
  and from there to the bivariate residuals: `"profile"`, the default,
  assesses each profile separately, `"overall"` pools them.

- rows:

  How many rows of each of the summary's tables to print, passed to
  `print(summary(x))`. A first look at a fit with thousands of
  observations would otherwise be mostly posteriors.

- ...:

  Nothing further is accepted. An argument this function cannot forward
  raises an error of class `latents_bad_argument` naming it, before
  anything has been printed, rather than being dropped.

## Value

The fitted model, invisibly. Called for the printing and drawing.

## What it prints

[`summary()`](https://rdrr.io/r/base/summary.html), which is every table
the fit can produce, then
[`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md),
then the condensed reading of
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md),
then every plot view the fit supports. Each section is what that verb
returns, and the tables are reached from
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
rather than from here.

## See also

[`summary()`](https://rdrr.io/r/base/summary.html),
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md),
[`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md),
[`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md).

## Examples

``` r
fit <- multilpa(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
  seed = 1
)
report(fit, plots = FALSE)
#> Multilevel LPA: 2 profiles and 2 group classes
#> Individuals: 1422; groups: 106; parameters: 23; converged: TRUE
#> Log likelihood: -8439.816424; AIC: 16925.633
#> BIC (groups): 16986.892; BIC (individuals): 17046.609
#> Best likelihood replicated in 4/4 starts (absolute tolerance 0.00844).
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator    mean variance standard_deviation mean_standard_error
#>        1     browse  0.5396   0.5899             0.7680             0.02710
#>        1   lectures  0.4274   0.8447             0.9191             0.03249
#>        1 forum_read  0.6198   0.4815             0.6939             0.02473
#>        1 forum_post  0.5305   0.6888             0.8299             0.02953
#>        1 attendance  0.6490   0.3841             0.6198             0.02226
#>        2     browse -0.7656   0.5285             0.7270             0.03121
#>        2   lectures -0.6065   0.5384             0.7337             0.03093
#>        2 forum_read -0.8792   0.3631             0.6026             0.02622
#>        2 forum_post -0.7527   0.4212             0.6490             0.02768
#>        2 attendance -0.9207   0.3736             0.6113             0.02643
#>  variance_standard_error
#>                  0.02948
#>                  0.04210
#>                  0.02440
#>                  0.03482
#>                  0.01956
#>                  0.03222
#>                  0.03209
#>                  0.02279
#>                  0.02608
#>                  0.02279
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse     0.5899
#>        1   lectures      browse     0.0000
#>        1 forum_read      browse     0.0000
#>        1 forum_post      browse     0.0000
#>        1 attendance      browse     0.0000
#>        1     browse    lectures     0.0000
#>        1   lectures    lectures     0.8447
#>        1 forum_read    lectures     0.0000
#>        1 forum_post    lectures     0.0000
#>        1 attendance    lectures     0.0000
#>    ... 40 more rows.  get_results(x, what = "covariances")
#> 
#> -- profile_probabilities -------------------------------------------
#>  group_class profile probability group_class_probability
#>            1       1      0.7859                  0.6747
#>            1       2      0.2141                  0.6747
#>            2       1      0.1705                  0.3253
#>            2       2      0.8295                  0.3253
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1          834.10               0.5866
#>  individuals     2          587.90               0.4134
#>       groups     1           71.52               0.6747
#>       groups     2           34.48               0.3253
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 9.996e-01  TRUE
#>    2     1       1 9.851e-01  TRUE
#>    3     1       1 2.421e-06 FALSE
#>    4     1       1 3.863e-06 FALSE
#>    5     1       1 5.749e-06 FALSE
#>    6     1       1 2.278e-05 FALSE
#>    7     1       1 2.301e-05 FALSE
#>    8     1       1 1.107e-05 FALSE
#>    9     1       1 2.173e-06 FALSE
#>   10     1       1 9.510e-06 FALSE
#>    ... 2834 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13         -61.52           1 1.473e-05 FALSE
#>      2         13         -71.94           1 6.391e-08 FALSE
#>      3         14         -76.69           1 1.000e+00  TRUE
#>      4         13         -67.56           1 1.000e+00  TRUE
#>      5         14         -78.80           1 1.248e-06 FALSE
#>      6         12         -80.90           1 9.802e-01  TRUE
#>      7         15         -89.36           1 1.000e+00  TRUE
#>      8         14         -85.22           1 1.067e-03 FALSE
#>      9         14         -83.80           1 1.000e+00  TRUE
#>     10         15         -82.06           1 1.565e-06 FALSE
#>    ... 202 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student browse lectures forum_read forum_post attendance profile group_class
#>        1   0.73    -0.08       0.30       0.67       0.72       1           2
#>        1   0.68     0.53      -0.02      -0.55       0.70       1           2
#>        1  -1.51    -0.51      -0.91      -0.80      -0.97       2           2
#>        1   0.64    -1.02      -1.62      -1.20      -1.26       2           2
#>        1  -0.94    -0.37      -0.71      -0.28      -1.52       2           2
#>        1  -0.36    -0.46      -0.77      -0.43      -1.34       2           2
#>        1   0.38    -1.85      -0.84      -1.09      -1.12       2           2
#>        1  -1.44     0.66      -0.67      -1.67      -1.00       2           2
#>        1  -0.76    -0.43      -1.04      -0.69      -1.37       2           2
#>        1   0.26    -0.95      -0.51      -1.89      -1.45       2           2
#>  uncertainty posterior_profile_1 posterior_profile_2
#>    4.436e-04           9.996e-01           0.0004436
#>    1.485e-02           9.851e-01           0.0148536
#>    2.421e-06           2.421e-06           0.9999976
#>    3.863e-06           3.863e-06           0.9999961
#>    5.749e-06           5.749e-06           0.9999943
#>    2.278e-05           2.278e-05           0.9999772
#>    2.301e-05           2.301e-05           0.9999770
#>    1.107e-05           1.107e-05           0.9999889
#>    2.173e-06           2.173e-06           0.9999978
#>    9.510e-06           9.510e-06           0.9999905
#>    ... 1412 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     836           0.5879      834.10               0.5866
#>  individuals     2     586           0.4121      587.90               0.4134
#>       groups     1      71           0.6698       71.52               0.6747
#>       groups     2      35           0.3302       34.48               0.3253
#>  average_posterior odds_correct_classification
#>             0.9849                       46.08
#>             0.9818                       76.35
#>             0.9920                       59.89
#>             0.9690                       64.80
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        836          0.984935
#>  individuals              1     2        836          0.015065
#>  individuals              2     1        586          0.018243
#>  individuals              2     2        586          0.981757
#>       groups              1     1         71          0.992014
#>       groups              1     2         71          0.007986
#>       groups              2     1         35          0.031015
#>       groups              2     2         35          0.968985
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.98718
#>  individuals          1              2     0.01282
#>  individuals          2              1     0.02142
#>  individuals          2              2     0.97858
#>       groups          1              1     0.98482
#>       groups          1              2     0.01518
#>       groups          2              1     0.01644
#>       groups          2              2     0.98356
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.01327
#>  individuals    1              1     2 -0.01327
#>  individuals    2              1     1  1.01327
#>  individuals    2              1     2 -0.01327
#>  individuals    3              2     1 -0.02218
#>  individuals    3              2     2  1.02218
#>  individuals    4              2     1 -0.02218
#>  individuals    4              2     2  1.02218
#>  individuals    5              2     1 -0.02218
#>  individuals    5              2     2  1.02218
#>    ... 3046 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2    1422      60.940           0.9382
#>       groups         2     106       4.435           0.9396
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_1    lectures  attendance gaussian  0.23393        0  0.23393
#>  profile_2  forum_read  attendance gaussian  0.22869        0  0.22869
#>  profile_1  forum_read  attendance gaussian  0.20880        0  0.20880
#>  profile_1  forum_post  attendance gaussian  0.20314        0  0.20314
#>  profile_2      browse  attendance gaussian  0.19248        0  0.19248
#>  profile_1      browse  attendance gaussian  0.17312        0  0.17312
#>  profile_2    lectures  attendance gaussian  0.12018        0  0.12018
#>  profile_2  forum_post  attendance gaussian  0.11883        0  0.11883
#>  profile_2      browse  forum_read gaussian  0.10925        0  0.10925
#>  profile_2  forum_read  forum_post gaussian  0.06676        0  0.06676
#>  effective_n statistic df   p_value p_adjusted
#>        834.1     6.871 NA 6.374e-12  6.374e-12
#>        587.9     5.630 NA 1.799e-08  1.799e-08
#>        834.1     6.109 NA 1.001e-09  1.001e-09
#>        834.1     5.939 NA 2.867e-09  2.867e-09
#>        587.9     4.714 NA 2.429e-06  2.429e-06
#>        834.1     5.042 NA 4.619e-07  4.619e-07
#>        587.9     2.921 NA 3.493e-03  3.493e-03
#>        587.9     2.888 NA 3.883e-03  3.883e-03
#>        587.9     2.653 NA 7.983e-03  7.983e-03
#>        587.9     1.617 NA 1.059e-01  1.059e-01
#>    ... 10 more rows.  get_results(x, what = "residuals")
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>           -8440           23 16926 16952      16987          17047        16914
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             16974       17010           17070      17172          17404
#>  icl_groups icl_individual clc_groups clc_individual
#>       16996          17168      16889          17002
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>            1422          1422      106          2               2      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                   VVI           23                            23          -8440
#>    aic bic_groups bic_individual converged iterations boundary small_classes
#>  16926      16987          17047      TRUE          9    FALSE         FALSE
#>  best_start n_best_replicated
#>           1                 4
#> 
#> -- stages ----------------------------------------------------------
#>  stage group_classes fixed log_likelihood parameters
#>  joint             2  <NA>          -8440         23
#>  parameters_with_measurement converged
#>                           23      TRUE
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error boundary
#>      1          -8440      TRUE          9  <NA>    FALSE
#>      2          -8440      TRUE         12  <NA>    FALSE
#>      3          -8440      TRUE          9  <NA>    FALSE
#>      4          -8440      TRUE         12  <NA>    FALSE
#> 
#> -- data ------------------------------------------------------------
#>  student browse lectures forum_read forum_post attendance
#>        1   0.73    -0.08       0.30       0.67       0.72
#>        1   0.68     0.53      -0.02      -0.55       0.70
#>        1  -1.51    -0.51      -0.91      -0.80      -0.97
#>        1   0.64    -1.02      -1.62      -1.20      -1.26
#>        1  -0.94    -0.37      -0.71      -0.28      -1.52
#>        1  -0.36    -0.46      -0.77      -0.43      -1.34
#>        1   0.38    -1.85      -0.84      -1.09      -1.12
#>        1  -1.44     0.66      -0.67      -1.67      -1.00
#>        1  -0.76    -0.43      -1.04      -0.69      -1.37
#>        1   0.26    -0.95      -0.51      -1.89      -1.45
#>    ... 1412 more rows.  get_results(x, what = "data")
#> 
#> 19 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
#> 
#>     variable    n n_missing          mean        sd   min  max n_distinct
#> 1     browse 1422         0 -2.812940e-05 0.9890804 -3.14 2.49        398
#> 2   lectures 1422         0 -4.219409e-05 0.9889178 -2.91 3.22        399
#> 3 forum_read 1422         0  7.032349e-05 0.9890236 -2.73 2.81        400
#> 4 forum_post 1422         0 -2.109705e-05 0.9890038 -2.61 3.10        397
#> 5 attendance 1422         0  7.735584e-05 0.9889426 -2.45 2.41        401
#>   n_groups        icc
#> 1      106 0.19229259
#> 2      106 0.08779169
#> 3      106 0.22403384
#> 4      106 0.15034618
#> 5      106 0.22145308
#> 
#> Classification quality: 2 profiles, 2 group classes
#> 
#>   Relative entropy     individuals 0.938    groups 0.940
#>   Smallest class       individuals 587.9 (41.3%)    groups 34.5 (32.5%)
#>   Lowest avg posterior individuals 0.982    groups 0.969
#>   Largest residual     0.234  (lectures, attendance; profile_1)
#> 
#> Tables: get_results(x, what = "entropy" | "classification" | 
#>         "average_posteriors" | "residuals" | "all"). plot(x) draws them.
```
