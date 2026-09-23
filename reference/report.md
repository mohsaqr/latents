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
#> Log likelihood: -8439.816436; AIC: 16925.633
#> BIC (groups): 16986.892; BIC (individuals): 17046.609
#> Best likelihood replicated in 4/4 starts (absolute tolerance 0.00844).
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator    mean variance standard_deviation
#>        1     browse  0.5396   0.5899             0.7680
#>        1   lectures  0.4274   0.8447             0.9191
#>        1 forum_read  0.6198   0.4815             0.6939
#>        1 forum_post  0.5305   0.6888             0.8299
#>        1 attendance  0.6491   0.3841             0.6198
#>        2     browse -0.7656   0.5285             0.7270
#>        2   lectures -0.6065   0.5384             0.7337
#>        2 forum_read -0.8792   0.3631             0.6026
#>        2 forum_post -0.7527   0.4212             0.6490
#>        2 attendance -0.9206   0.3736             0.6113
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
#>            1       1      0.1704                  0.3252
#>            1       2      0.8296                  0.3252
#>            2       1      0.7859                  0.6748
#>            2       2      0.2141                  0.6748
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1          834.08               0.5866
#>  individuals     2          587.92               0.4134
#>       groups     1           34.48               0.3253
#>       groups     2           71.52               0.6747
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 9.996e-01  TRUE
#>    2     1       1 9.851e-01  TRUE
#>    3     1       1 2.418e-06 FALSE
#>    4     1       1 3.858e-06 FALSE
#>    5     1       1 5.743e-06 FALSE
#>    6     1       1 2.275e-05 FALSE
#>    7     1       1 2.298e-05 FALSE
#>    8     1       1 1.105e-05 FALSE
#>    9     1       1 2.171e-06 FALSE
#>   10     1       1 9.498e-06 FALSE
#>    ... 2834 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13         -61.52           1 1.000e+00  TRUE
#>      2         13         -71.94           1 1.000e+00  TRUE
#>      3         14         -76.69           1 5.976e-10 FALSE
#>      4         13         -67.56           1 1.321e-09 FALSE
#>      5         14         -78.80           1 1.000e+00  TRUE
#>      6         12         -80.90           1 1.979e-02 FALSE
#>      7         15         -89.36           1 1.118e-09 FALSE
#>      8         14         -85.22           1 9.989e-01  TRUE
#>      9         14         -83.80           1 4.287e-10 FALSE
#>     10         15         -82.06           1 1.000e+00  TRUE
#>    ... 202 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student browse lectures forum_read forum_post attendance profile group_class
#>        1   0.73    -0.08       0.30       0.67       0.72       1           1
#>        1   0.68     0.53      -0.02      -0.55       0.70       1           1
#>        1  -1.51    -0.51      -0.91      -0.80      -0.97       2           1
#>        1   0.64    -1.02      -1.62      -1.20      -1.26       2           1
#>        1  -0.94    -0.37      -0.71      -0.28      -1.52       2           1
#>        1  -0.36    -0.46      -0.77      -0.43      -1.34       2           1
#>        1   0.38    -1.85      -0.84      -1.09      -1.12       2           1
#>        1  -1.44     0.66      -0.67      -1.67      -1.00       2           1
#>        1  -0.76    -0.43      -1.04      -0.69      -1.37       2           1
#>        1   0.26    -0.95      -0.51      -1.89      -1.45       2           1
#>  uncertainty posterior_profile_1 posterior_profile_2
#>    4.442e-04           9.996e-01           0.0004442
#>    1.487e-02           9.851e-01           0.0148709
#>    2.418e-06           2.418e-06           0.9999976
#>    3.858e-06           3.858e-06           0.9999961
#>    5.743e-06           5.743e-06           0.9999943
#>    2.275e-05           2.275e-05           0.9999772
#>    2.298e-05           2.298e-05           0.9999770
#>    1.105e-05           1.105e-05           0.9999889
#>    2.171e-06           2.171e-06           0.9999978
#>    9.498e-06           9.498e-06           0.9999905
#>    ... 1412 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     836           0.5879      834.08               0.5866
#>  individuals     2     586           0.4121      587.92               0.4134
#>       groups     1      35           0.3302       34.48               0.3253
#>       groups     2      71           0.6698       71.52               0.6747
#>  average_posterior odds_correct_classification
#>             0.9849                       46.05
#>             0.9818                       76.39
#>             0.9689                       64.68
#>             0.9920                       60.05
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        836          0.984925
#>  individuals              1     2        836          0.015075
#>  individuals              2     1        586          0.018232
#>  individuals              2     2        586          0.981768
#>       groups              1     1         35          0.968922
#>       groups              1     2         35          0.031078
#>       groups              2     1         71          0.007963
#>       groups              2     2         71          0.992037
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.98719
#>  individuals          1              2     0.01281
#>  individuals          2              1     0.02144
#>  individuals          2              2     0.97856
#>       groups          1              1     0.98360
#>       groups          1              2     0.01640
#>       groups          2              1     0.01521
#>       groups          2              2     0.98479
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.01326
#>  individuals    1              1     2 -0.01326
#>  individuals    2              1     1  1.01326
#>  individuals    2              1     2 -0.01326
#>  individuals    3              2     1 -0.02220
#>  individuals    3              2     2  1.02220
#>  individuals    4              2     1 -0.02220
#>  individuals    4              2     2  1.02220
#>  individuals    5              2     1 -0.02220
#>  individuals    5              2     2  1.02220
#>    ... 3046 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2    1422      60.943           0.9382
#>       groups         2     106       4.433           0.9397
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_1    lectures  attendance gaussian  0.23392        0  0.23392
#>  profile_2  forum_read  attendance gaussian  0.22871        0  0.22871
#>  profile_1  forum_read  attendance gaussian  0.20879        0  0.20879
#>  profile_1  forum_post  attendance gaussian  0.20313        0  0.20313
#>  profile_2      browse  attendance gaussian  0.19250        0  0.19250
#>  profile_1      browse  attendance gaussian  0.17311        0  0.17311
#>  profile_2    lectures  attendance gaussian  0.12018        0  0.12018
#>  profile_2  forum_post  attendance gaussian  0.11884        0  0.11884
#>  profile_2      browse  forum_read gaussian  0.10928        0  0.10928
#>  profile_2  forum_read  forum_post gaussian  0.06677        0  0.06677
#>  effective_n statistic df   p_value p_adjusted
#>        834.1     6.871 NA 6.393e-12  6.393e-12
#>        587.9     5.631 NA 1.793e-08  1.793e-08
#>        834.1     6.109 NA 1.003e-09  1.003e-09
#>        834.1     5.939 NA 2.875e-09  2.875e-09
#>        587.9     4.715 NA 2.423e-06  2.423e-06
#>        834.1     5.041 NA 4.626e-07  4.626e-07
#>        587.9     2.921 NA 3.491e-03  3.491e-03
#>        587.9     2.888 NA 3.879e-03  3.879e-03
#>        587.9     2.653 NA 7.968e-03  7.968e-03
#>        587.9     1.617 NA 1.058e-01  1.058e-01
#>    ... 10 more rows.  get_results(x, what = "residuals")
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>           -8440           23 16926 16952      16987          17047        16914
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             16974       17010           17070      17172          17404
#>  icl_groups icl_individual clc_groups clc_individual
#>       16996          17168      16888          17002
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>            1422          1422      106          2               2      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                   VVI           23                            23          -8440
#>    aic bic_groups bic_individual converged iterations boundary small_classes
#>  16926      16987          17047      TRUE         10    FALSE         FALSE
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
#>      1          -8440      TRUE         10  <NA>    FALSE
#>      2          -8440      TRUE         10  <NA>    FALSE
#>      3          -8440      TRUE         10  <NA>    FALSE
#>      4          -8440      TRUE         13  <NA>    FALSE
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
