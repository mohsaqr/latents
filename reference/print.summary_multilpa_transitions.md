# Print a latent transition summary

Print a latent transition summary

## Usage

``` r
# S3 method for class 'summary_multilpa_transitions'
print(x, digits = 4L, rows = 10L, ...)
```

## Arguments

- x:

  A `summary_multilpa_transitions` object.

- digits:

  Number of printed significant digits.

- rows:

  How many rows of each table to print. A longer table is shown to that
  depth, with its remaining row count and the
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  call that returns it whole.

- ...:

  Additional arguments passed to matrix printing.

## Value

The summary, invisibly; called for what it prints.

## Examples

``` r
fit <- lta(subset(course_engagement, student <= 40),
           c("browse", "lectures", "forum_read"), "student",
           n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
print(summary(fit), digits = 3)
#> Latent transition model: 2 profiles and 1 group class
#> Observations: 547; groups: 40; up to 15 occasions (unbalanced, observed grid)
#> Parameters: 15; converged: TRUE
#> Log likelihood: -2025.406100; AIC: 4080.812
#> BIC (groups): 4106.145; BIC (individuals): 4145.379
#> Best likelihood replicated in 2/2 starts (absolute tolerance 0.00203).
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator   mean variance standard_deviation mean_standard_error
#>        1     browse  0.555    0.549              0.741              0.0402
#>        1   lectures  0.437    0.825              0.908              0.0489
#>        1 forum_read  0.609    0.499              0.706              0.0391
#>        2     browse -0.871    0.460              0.678              0.0519
#>        2   lectures -0.626    0.574              0.757              0.0571
#>        2 forum_read -0.932    0.326              0.571              0.0436
#>  variance_standard_error
#>                   0.0422
#>                   0.0625
#>                   0.0396
#>                   0.0494
#>                   0.0611
#>                   0.0354
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse      0.549
#>        1   lectures      browse      0.000
#>        1 forum_read      browse      0.000
#>        1     browse    lectures      0.000
#>        1   lectures    lectures      0.825
#>        1 forum_read    lectures      0.000
#>        1     browse  forum_read      0.000
#>        1   lectures  forum_read      0.000
#>        1 forum_read  forum_read      0.499
#>        2     browse      browse      0.460
#>    ... 8 more rows.  get_results(x, what = "covariances")
#> 
#> -- transitions -----------------------------------------------------
#>  group_class from to probability expected_count stable estimated
#>            1    1  1       0.891          297.4   TRUE      TRUE
#>            1    1  2       0.109           36.5  FALSE      TRUE
#>            1    2  1       0.190           32.8  FALSE      TRUE
#>            1    2  2       0.810          140.2   TRUE      TRUE
#>  group_class_probability
#>                        1
#>                        1
#>                        1
#>                        1
#> 
#> -- initial ---------------------------------------------------------
#>  group_class profile probability prevalence group_class_probability
#>            1       1       0.742      0.658                       1
#>            1       2       0.258      0.342                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1             360                0.658
#>  individuals     2             187                0.342
#>       groups     1              40                1.000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1  0.998908  TRUE
#>    2     1       1  0.985600  TRUE
#>    3     1       1  0.001292 FALSE
#>    4     1       1  0.001097 FALSE
#>    5     1       1  0.000372 FALSE
#>    6     1       1  0.001270 FALSE
#>    7     1       1  0.002108 FALSE
#>    8     1       1  0.000705 FALSE
#>    9     1       1  0.000244 FALSE
#>   10     1       1  0.009924 FALSE
#>    ... 1084 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13          -39.8           1         1  TRUE
#>      2         13          -45.2           1         1  TRUE
#>      3         14          -47.0           1         1  TRUE
#>      4         13          -42.3           1         1  TRUE
#>      5         14          -45.7           1         1  TRUE
#>      6         12          -51.3           1         1  TRUE
#>      7         15          -58.1           1         1  TRUE
#>      8         14          -52.3           1         1  TRUE
#>      9         14          -56.5           1         1  TRUE
#>     10         15          -57.9           1         1  TRUE
#>    ... 30 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student sequence browse lectures forum_read profile group_class uncertainty
#>        1        1   0.73    -0.08       0.30       1           1    0.001092
#>        1        2   0.68     0.53      -0.02       1           1    0.014400
#>        1        3  -1.51    -0.51      -0.91       2           1    0.001292
#>        1        4   0.64    -1.02      -1.62       2           1    0.001097
#>        1        5  -0.94    -0.37      -0.71       2           1    0.000372
#>        1        6  -0.36    -0.46      -0.77       2           1    0.001270
#>        1        7   0.38    -1.85      -0.84       2           1    0.002108
#>        1        8  -1.44     0.66      -0.67       2           1    0.000705
#>        1        9  -0.76    -0.43      -1.04       2           1    0.000244
#>        1       10   0.26    -0.95      -0.51       2           1    0.009924
#>  posterior_profile_1 posterior_profile_2
#>             0.998908             0.00109
#>             0.985600             0.01440
#>             0.001292             0.99871
#>             0.001097             0.99890
#>             0.000372             0.99963
#>             0.001270             0.99873
#>             0.002108             0.99789
#>             0.000705             0.99929
#>             0.000244             0.99976
#>             0.009924             0.99008
#>    ... 537 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     358            0.654         360                0.658
#>  individuals     2     189            0.346         187                0.342
#>       groups     1      40            1.000          40                1.000
#>  average_posterior odds_correct_classification
#>              0.988                        44.3
#>              0.968                        58.0
#>              1.000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        358            0.9884
#>  individuals              1     2        358            0.0116
#>  individuals              2     1        189            0.0321
#>  individuals              2     2        189            0.9679
#>       groups              1     1         40            1.0000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1      0.9831
#>  individuals          1              2      0.0169
#>  individuals          2              1      0.0222
#>  individuals          2              2      0.9778
#>       groups          1              1      1.0000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class  weight
#>  individuals    1              1     1  1.0175
#>  individuals    1              1     2 -0.0175
#>  individuals    2              1     1  1.0175
#>  individuals    2              1     2 -0.0175
#>  individuals    3              2     1 -0.0231
#>  individuals    3              2     2  1.0231
#>  individuals    4              2     1 -0.0231
#>  individuals    4              2     2  1.0231
#>  individuals    5              2     1 -0.0231
#>  individuals    5              2     2  1.0231
#>    ... 1124 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     547        29.1            0.923
#>       groups         1      40         0.0               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2      browse    lectures gaussian  -0.1508        0  -0.1508
#>  profile_2      browse  forum_read gaussian   0.1218        0   0.1218
#>  profile_1    lectures  forum_read gaussian   0.0708        0   0.0708
#>  profile_1      browse  forum_read gaussian   0.0431        0   0.0431
#>  profile_2    lectures  forum_read gaussian  -0.0384        0  -0.0384
#>  profile_1      browse    lectures gaussian   0.0370        0   0.0370
#>  effective_n statistic df p_value p_adjusted
#>          187    -2.062 NA  0.0392     0.0392
#>          187     1.660 NA  0.0969     0.0969
#>          360     1.341 NA  0.1800     0.1800
#>          360     0.814 NA  0.4154     0.4154
#>          187    -0.522 NA  0.6018     0.6018
#>          360     0.700 NA  0.4840     0.4840
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters  aic  kic bic_groups bic_individual sabic_groups
#>           -2025           15 4081 4099       4106           4145         4059
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>              4098        4121            4160       4206           4343
#>  icl_groups icl_individual clc_groups clc_individual
#>        4106           4204       4051           4109
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             547           547       40          2               1      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                  <NA>           15                            15          -2025
#>   aic bic_groups bic_individual converged iterations boundary small_classes
#>  4081       4106           4145      TRUE          8    FALSE         FALSE
#>  best_start n_best_replicated
#>           1                 2
#> 
#> -- sequences -------------------------------------------------------
#>  group group_class time profile
#>      1           1    1       1
#>      1           1    2       1
#>      1           1    3       2
#>      1           1    4       2
#>      1           1    5       2
#>      1           1    6       2
#>      1           1    7       2
#>      1           1    8       2
#>      1           1    9       2
#>      1           1   10       2
#>    ... 537 more rows.  get_results(x, what = "sequences")
#> 
#> -- sequence_summary ------------------------------------------------
#>  group_class groups observations mean_length median_length shortest longest
#>            1     40          547        13.7            14       11      15
#>  complete gaps
#>        10    0
#> 
#> -- sequence_lengths ------------------------------------------------
#>  group group_class occasions observations complete
#>      1           1        13           13    FALSE
#>      2           1        13           13    FALSE
#>      3           1        14           14    FALSE
#>      4           1        13           13    FALSE
#>      5           1        14           14    FALSE
#>      6           1        12           12    FALSE
#>      7           1        15           15     TRUE
#>      8           1        14           14    FALSE
#>      9           1        14           14    FALSE
#>     10           1        15           15     TRUE
#>    ... 30 more rows.  get_results(x, what = "sequence_lengths")
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error
#>      1          -2025      TRUE          8  <NA>
#>      2          -2025      TRUE          9  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  student sequence browse lectures forum_read
#>        1        1   0.73    -0.08       0.30
#>        1        2   0.68     0.53      -0.02
#>        1        3  -1.51    -0.51      -0.91
#>        1        4   0.64    -1.02      -1.62
#>        1        5  -0.94    -0.37      -0.71
#>        1        6  -0.36    -0.46      -0.77
#>        1        7   0.38    -1.85      -0.84
#>        1        8  -1.44     0.66      -0.67
#>        1        9  -0.76    -0.43      -1.04
#>        1       10   0.26    -0.95      -0.51
#>    ... 537 more rows.  get_results(x, what = "data")
#> 
#> 22 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
