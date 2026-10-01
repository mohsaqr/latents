# Print a covariate LPA summary

Print a covariate LPA summary

## Usage

``` r
# S3 method for class 'summary_multilpa_covariates'
print(x, digits = 4L, rows = 10L, ...)
```

## Arguments

- x:

  A `summary_multilpa_covariates` object.

- digits:

  Number of printed significant digits.

- rows:

  How many rows of each table to print. A longer table is shown to that
  depth, with its remaining row count and the
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  call that returns it whole.

- ...:

  Passed to the underlying `data.frame` printing.

## Value

The summary, invisibly. Called for the side effect of printing the model
line, the measurement model, the membership regressions, the effective
memberships at both levels, the likelihood and information criteria, any
warnings, and the restart diagnostics.

## Examples

``` r
fit <- multilpa(subset(course_engagement, student <= 40),
                c("browse", "lectures", "forum_read"), "student",
                n_profiles = 2, n_group_classes = 1,
                profile_covariates = "previous_grade", n_starts = 2, seed = 1)
print(summary(fit), digits = 3)
#> Multilevel LPA with covariates: 2 profiles, 1 group classes
#> Individuals: 547; groups: 40; parameters: 14; converged: TRUE
#> 1 profile covariate(s); 0 group covariate(s)
#> Log likelihood: -2099.951166; AIC: 4227.902
#> BIC (groups): 4251.547; BIC (individuals): 4288.165
#> Membership coefficients carry no standard errors here; parameter_inference() has them.
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator   mean variance standard_deviation mean_standard_error
#>        1     browse  0.567    0.542              0.736              0.0420
#>        1   lectures  0.458    0.814              0.902              0.0509
#>        1 forum_read  0.631    0.474              0.688              0.0405
#>        2     browse -0.840    0.489              0.699              0.0574
#>        2   lectures -0.625    0.565              0.752              0.0583
#>        2 forum_read -0.916    0.339              0.582              0.0488
#>  variance_standard_error
#>                   0.0438
#>                   0.0641
#>                   0.0393
#>                   0.0572
#>                   0.0625
#>                   0.0417
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse      0.542
#>        1   lectures      browse      0.000
#>        1 forum_read      browse      0.000
#>        1     browse    lectures      0.000
#>        1   lectures    lectures      0.814
#>        1 forum_read    lectures      0.000
#>        1     browse  forum_read      0.000
#>        1   lectures  forum_read      0.000
#>        1 forum_read  forum_read      0.474
#>        2     browse      browse      0.489
#>    ... 8 more rows.  get_results(x, what = "covariances")
#> 
#> -- coefficients ----------------------------------------------------
#>    level   outcome           term parameter estimate
#>  profile profile_1    (Intercept)     logit    0.614
#>  profile profile_1 previous_grade     logit    0.651
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1             353                0.645
#>  individuals     2             194                0.355
#>       groups     1              40                1.000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1   0.99459  TRUE
#>    2     1       1   0.98232  TRUE
#>    3     1       1   0.00613 FALSE
#>    4     1       1   0.01348 FALSE
#>    5     1       1   0.01092 FALSE
#>    6     1       1   0.04160 FALSE
#>    7     1       1   0.08375 FALSE
#>    8     1       1   0.02065 FALSE
#>    9     1       1   0.00324 FALSE
#>   10     1       1   0.11097 FALSE
#>    ... 1084 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13          -46.0           1         1  TRUE
#>      2         13          -52.0           1         1  TRUE
#>      3         14          -49.0           1         1  TRUE
#>      4         13          -46.3           1         1  TRUE
#>      5         14          -52.8           1         1  TRUE
#>      6         12          -50.4           1         1  TRUE
#>      7         15          -59.4           1         1  TRUE
#>      8         14          -56.4           1         1  TRUE
#>      9         14          -60.2           1         1  TRUE
#>     10         15          -61.1           1         1  TRUE
#>    ... 30 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student browse lectures forum_read profile group_class uncertainty
#>        1   0.73    -0.08       0.30       1           1     0.00541
#>        1   0.68     0.53      -0.02       1           1     0.01768
#>        1  -1.51    -0.51      -0.91       2           1     0.00613
#>        1   0.64    -1.02      -1.62       2           1     0.01348
#>        1  -0.94    -0.37      -0.71       2           1     0.01092
#>        1  -0.36    -0.46      -0.77       2           1     0.04160
#>        1   0.38    -1.85      -0.84       2           1     0.08375
#>        1  -1.44     0.66      -0.67       2           1     0.02065
#>        1  -0.76    -0.43      -1.04       2           1     0.00324
#>        1   0.26    -0.95      -0.51       2           1     0.11097
#>  posterior_profile_1 posterior_profile_2
#>              0.99459             0.00541
#>              0.98232             0.01768
#>              0.00613             0.99387
#>              0.01348             0.98652
#>              0.01092             0.98908
#>              0.04160             0.95840
#>              0.08375             0.91625
#>              0.02065             0.97935
#>              0.00324             0.99676
#>              0.11097             0.88903
#>    ... 537 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     352            0.644         353                0.645
#>  individuals     2     195            0.356         194                0.355
#>       groups     1      40            1.000          40                1.000
#>  average_posterior odds_correct_classification
#>              0.973                        20.1
#>              0.948                        32.8
#>              1.000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        352            0.9734
#>  individuals              1     2        352            0.0266
#>  individuals              2     1        195            0.0525
#>  individuals              2     2        195            0.9475
#>       groups              1     1         40            1.0000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1      0.9710
#>  individuals          1              2      0.0290
#>  individuals          2              1      0.0483
#>  individuals          2              2      0.9517
#>       groups          1              1      1.0000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class  weight
#>  individuals    1              1     1  1.0314
#>  individuals    1              1     2 -0.0314
#>  individuals    2              1     1  1.0314
#>  individuals    2              1     2 -0.0314
#>  individuals    3              2     1 -0.0523
#>  individuals    3              2     2  1.0523
#>  individuals    4              2     1 -0.0523
#>  individuals    4              2     2  1.0523
#>  individuals    5              2     1 -0.0523
#>  individuals    5              2     2  1.0523
#>    ... 1124 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     547        55.2            0.854
#>       groups         1      40         0.0               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2      browse    lectures gaussian  -0.1352        0  -0.1352
#>  profile_2      browse  forum_read gaussian   0.1226        0   0.1226
#>  profile_2    lectures  forum_read gaussian  -0.0557        0  -0.0557
#>  profile_1    lectures  forum_read gaussian   0.0455        0   0.0455
#>  profile_1      browse  forum_read gaussian   0.0283        0   0.0283
#>  profile_1      browse    lectures gaussian   0.0181        0   0.0181
#>  effective_n statistic df p_value p_adjusted
#>          194    -1.880 NA  0.0601     0.0601
#>          194     1.704 NA  0.0883     0.0883
#>          194    -0.771 NA  0.4406     0.4406
#>          353     0.852 NA  0.3939     0.3939
#>          353     0.529 NA  0.5967     0.5967
#>          353     0.339 NA  0.7348     0.7348
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters  aic  kic bic_groups bic_individual sabic_groups
#>           -2100           14 4228 4245       4252           4288         4208
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>              4244        4266            4302       4345           4529
#>  icl_groups icl_individual clc_groups clc_individual
#>        4252           4399       4200           4310
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_groups n_profiles n_group_classes n_profile_covariates
#>             547       40          2               1                    1
#>  n_group_covariates variance_model covariance_model n_parameters log_likelihood
#>                   0        varying         diagonal           14          -2100
#>   aic bic_groups bic_individual converged boundary extreme_logits n_starts
#>  4228       4252           4288      TRUE    FALSE          FALSE        2
#>  weights
#>     <NA>
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error
#>      1          -2100      TRUE          7  <NA>
#>      2          -2100      TRUE         15  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  student browse lectures forum_read
#>        1   0.73    -0.08       0.30
#>        1   0.68     0.53      -0.02
#>        1  -1.51    -0.51      -0.91
#>        1   0.64    -1.02      -1.62
#>        1  -0.94    -0.37      -0.71
#>        1  -0.36    -0.46      -0.77
#>        1   0.38    -1.85      -0.84
#>        1  -1.44     0.66      -0.67
#>        1  -0.76    -0.43      -1.04
#>        1   0.26    -0.95      -0.51
#>    ... 537 more rows.  get_results(x, what = "data")
#> 
#> 18 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
