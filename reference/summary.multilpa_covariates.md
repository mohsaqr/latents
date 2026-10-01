# Summarize a covariate LPA fit

Collects the model-level fit, the measurement model, the membership
regressions and the restart diagnostics into one object, so that none of
them has to be read out of the fit by hand.

## Usage

``` r
# S3 method for class 'multilpa_covariates'
summary(object, ...)
```

## Arguments

- object:

  A covariate LPA fit from `multilpa(profile_covariates = )`.

- ...:

  Reserved for compatibility with
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `summary_multilpa_covariates`, with a `print` method
and an [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
accessor. [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
returns the one-row model summary by default; `what = "profiles"`,
`"coefficients"` and `"starts"` return the measurement model, the
membership coefficients and the restart diagnostics. The membership
coefficients carry no standard errors here;
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports those.

## Examples

``` r
fit <- multilpa(subset(course_engagement, student <= 40),
                c("browse", "lectures", "forum_read"), "student",
                n_profiles = 2, n_group_classes = 1,
                profile_covariates = "previous_grade", n_starts = 2, seed = 1)
summary(fit)
#> Multilevel LPA with covariates: 2 profiles, 1 group classes
#> Individuals: 547; groups: 40; parameters: 14; converged: TRUE
#> 1 profile covariate(s); 0 group covariate(s)
#> Log likelihood: -2099.951166; AIC: 4227.902
#> BIC (groups): 4251.547; BIC (individuals): 4288.165
#> Membership coefficients carry no standard errors here; parameter_inference() has them.
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator    mean variance standard_deviation mean_standard_error
#>        1     browse  0.5669   0.5419             0.7361             0.04201
#>        1   lectures  0.4576   0.8137             0.9020             0.05091
#>        1 forum_read  0.6311   0.4740             0.6885             0.04049
#>        2     browse -0.8400   0.4887             0.6991             0.05736
#>        2   lectures -0.6249   0.5655             0.7520             0.05835
#>        2 forum_read -0.9163   0.3388             0.5820             0.04877
#>  variance_standard_error
#>                  0.04383
#>                  0.06406
#>                  0.03928
#>                  0.05724
#>                  0.06253
#>                  0.04172
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse     0.5419
#>        1   lectures      browse     0.0000
#>        1 forum_read      browse     0.0000
#>        1     browse    lectures     0.0000
#>        1   lectures    lectures     0.8137
#>        1 forum_read    lectures     0.0000
#>        1     browse  forum_read     0.0000
#>        1   lectures  forum_read     0.0000
#>        1 forum_read  forum_read     0.4740
#>        2     browse      browse     0.4887
#>    ... 8 more rows.  get_results(x, what = "covariances")
#> 
#> -- coefficients ----------------------------------------------------
#>    level   outcome           term parameter estimate
#>  profile profile_1    (Intercept)     logit   0.6143
#>  profile profile_1 previous_grade     logit   0.6507
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1           352.9               0.6451
#>  individuals     2           194.1               0.3549
#>       groups     1            40.0               1.0000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1  0.994590  TRUE
#>    2     1       1  0.982315  TRUE
#>    3     1       1  0.006128 FALSE
#>    4     1       1  0.013476 FALSE
#>    5     1       1  0.010921 FALSE
#>    6     1       1  0.041605 FALSE
#>    7     1       1  0.083749 FALSE
#>    8     1       1  0.020653 FALSE
#>    9     1       1  0.003240 FALSE
#>   10     1       1  0.110966 FALSE
#>    ... 1084 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13         -46.02           1         1  TRUE
#>      2         13         -52.04           1         1  TRUE
#>      3         14         -49.04           1         1  TRUE
#>      4         13         -46.34           1         1  TRUE
#>      5         14         -52.81           1         1  TRUE
#>      6         12         -50.38           1         1  TRUE
#>      7         15         -59.39           1         1  TRUE
#>      8         14         -56.36           1         1  TRUE
#>      9         14         -60.18           1         1  TRUE
#>     10         15         -61.09           1         1  TRUE
#>    ... 30 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student browse lectures forum_read profile group_class uncertainty
#>        1   0.73    -0.08       0.30       1           1    0.005410
#>        1   0.68     0.53      -0.02       1           1    0.017685
#>        1  -1.51    -0.51      -0.91       2           1    0.006128
#>        1   0.64    -1.02      -1.62       2           1    0.013476
#>        1  -0.94    -0.37      -0.71       2           1    0.010921
#>        1  -0.36    -0.46      -0.77       2           1    0.041605
#>        1   0.38    -1.85      -0.84       2           1    0.083749
#>        1  -1.44     0.66      -0.67       2           1    0.020653
#>        1  -0.76    -0.43      -1.04       2           1    0.003240
#>        1   0.26    -0.95      -0.51       2           1    0.110966
#>  posterior_profile_1 posterior_profile_2
#>             0.994590             0.00541
#>             0.982315             0.01768
#>             0.006128             0.99387
#>             0.013476             0.98652
#>             0.010921             0.98908
#>             0.041605             0.95840
#>             0.083749             0.91625
#>             0.020653             0.97935
#>             0.003240             0.99676
#>             0.110966             0.88903
#>    ... 537 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     352           0.6435       352.9               0.6451
#>  individuals     2     195           0.3565       194.1               0.3549
#>       groups     1      40           1.0000        40.0               1.0000
#>  average_posterior odds_correct_classification
#>             0.9734                       20.12
#>             0.9475                       32.81
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        352           0.97338
#>  individuals              1     2        352           0.02662
#>  individuals              2     1        195           0.05248
#>  individuals              2     2        195           0.94752
#>       groups              1     1         40           1.00000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.97100
#>  individuals          1              2     0.02900
#>  individuals          2              1     0.04827
#>  individuals          2              2     0.95173
#>       groups          1              1     1.00000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.03143
#>  individuals    1              1     2 -0.03143
#>  individuals    2              1     1  1.03143
#>  individuals    2              1     2 -0.03143
#>  individuals    3              2     1 -0.05231
#>  individuals    3              2     2  1.05231
#>  individuals    4              2     1 -0.05231
#>  individuals    4              2     2  1.05231
#>  individuals    5              2     1 -0.05231
#>  individuals    5              2     2  1.05231
#>    ... 1124 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     547       55.18           0.8545
#>       groups         1      40        0.00               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2      browse    lectures gaussian -0.13518        0 -0.13518
#>  profile_2      browse  forum_read gaussian  0.12265        0  0.12265
#>  profile_2    lectures  forum_read gaussian -0.05573        0 -0.05573
#>  profile_1    lectures  forum_read gaussian  0.04555        0  0.04555
#>  profile_1      browse  forum_read gaussian  0.02828        0  0.02828
#>  profile_1      browse    lectures gaussian  0.01811        0  0.01811
#>  effective_n statistic df p_value p_adjusted
#>        194.1   -1.8804 NA 0.06006    0.06006
#>        194.1    1.7042 NA 0.08835    0.08835
#>        194.1   -0.7712 NA 0.44056    0.44056
#>        352.9    0.8525 NA 0.39394    0.39394
#>        352.9    0.5291 NA 0.59671    0.59671
#>        352.9    0.3388 NA 0.73479    0.73479
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
get_results(fit, what = "coefficients")
#>     level   outcome           term parameter  estimate
#> 1 profile profile_1    (Intercept)     logit 0.6143165
#> 2 profile profile_1 previous_grade     logit 0.6506674
```
