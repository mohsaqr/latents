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
set.seed(1)
example_data <- data.frame(group = rep(seq_len(20), each = 10),
                           z = rnorm(200))
example_data$y <- rnorm(200,
  ifelse(runif(200) < plogis(example_data$z), -3, 3))
fit <- multilpa(example_data, "y", "group", n_profiles = 2,
                n_group_classes = 1, profile_covariates = "z",
                n_starts = 2, seed = 1)
summary(fit)
#> Multilevel LPA with covariates: 2 profiles, 1 group classes
#> Individuals: 200; groups: 20; parameters: 6; converged: TRUE
#> 1 profile covariate(s); 0 group covariate(s)
#> Log likelihood: -425.206918; AIC: 862.414
#> BIC (groups): 868.388; BIC (individuals): 882.204
#> Membership coefficients carry no standard errors here; parameter_inference() has them.
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator   mean variance standard_deviation
#>        1         y  3.009    1.235              1.111
#>        2         y -3.005    1.114              1.055
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile indicator indicator_2 covariance
#>        1         y           y      1.235
#>        2         y           y      1.114
#> 
#> -- coefficients ----------------------------------------------------
#>    level   outcome          term parameter  estimate
#>  profile profile_1 group_class_1     logit  0.004333
#>  profile profile_1             z     logit -0.826689
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1           99.15               0.4958
#>  individuals     2          100.85               0.5042
#>       groups     1           20.00               1.0000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 1.000e+00  TRUE
#>    2     1       1 2.358e-09 FALSE
#>    3     1       1 1.000e+00  TRUE
#>    4     1       1 1.000e+00  TRUE
#>    5     1       1 1.000e+00  TRUE
#>    6     1       1 1.000e+00  TRUE
#>    7     1       1 4.294e-07 FALSE
#>    8     1       1 1.000e+00  TRUE
#>    9     1       1 1.901e-09 FALSE
#>   10     1       1 1.000e+00  TRUE
#>    ... 390 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         10         -22.07           1         1  TRUE
#>      2         10         -18.94           1         1  TRUE
#>      3         10         -19.53           1         1  TRUE
#>      4         10         -20.69           1         1  TRUE
#>      5         10         -22.77           1         1  TRUE
#>      6         10         -20.35           1         1  TRUE
#>      7         10         -20.38           1         1  TRUE
#>      8         10         -18.84           1         1  TRUE
#>      9         10         -21.99           1         1  TRUE
#>     10         10         -19.17           1         1  TRUE
#>    ... 10 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  group      y profile group_class uncertainty posterior_profile_1
#>      1  3.894       1           1   4.514e-10           1.000e+00
#>      1 -4.047       2           1   2.358e-09           2.358e-09
#>      1  4.971       1           1   9.857e-13           1.000e+00
#>      1  2.616       1           1   2.879e-06           1.000e+00
#>      1  4.654       1           1   1.506e-11           1.000e+00
#>      1  4.512       1           1   1.277e-11           1.000e+00
#>      1 -2.917       2           1   4.294e-07           4.294e-07
#>      1  3.567       1           1   8.299e-09           1.000e+00
#>      1 -4.025       2           1   1.901e-09           1.901e-09
#>      1  3.323       1           1   1.321e-08           1.000e+00
#>  posterior_profile_2
#>            4.514e-10
#>            1.000e+00
#>            9.856e-13
#>            2.879e-06
#>            1.506e-11
#>            1.277e-11
#>            1.000e+00
#>            8.299e-09
#>            1.000e+00
#>            1.321e-08
#>    ... 190 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1      99            0.495       99.15               0.4958
#>  individuals     2     101            0.505      100.85               0.5042
#>       groups     1      20            1.000       20.00               1.0000
#>  average_posterior odds_correct_classification
#>             0.9998                      5148.3
#>             0.9983                       570.9
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1         99         0.9998025
#>  individuals              1     2         99         0.0001975
#>  individuals              2     1        101         0.0017193
#>  individuals              2     2        101         0.9982807
#>       groups              1     1         20         1.0000000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1   0.9982487
#>  individuals          1              2   0.0017513
#>  individuals          2              1   0.0001939
#>  individuals          2              2   0.9998061
#>       groups          1              1   1.0000000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class     weight
#>  individuals    1              1     1  1.0017547
#>  individuals    1              1     2 -0.0017547
#>  individuals    2              2     1 -0.0001943
#>  individuals    2              2     2  1.0001943
#>  individuals    3              1     1  1.0017547
#>  individuals    3              1     2 -0.0017547
#>  individuals    4              1     1  1.0017547
#>  individuals    4              1     2 -0.0017547
#>  individuals    5              1     1  1.0017547
#>  individuals    5              1     2 -0.0017547
#>    ... 410 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     200      0.8458           0.9939
#>       groups         1      20      0.0000               NA
#> 
#> -- residuals -------------------------------------------------------
#>    (no rows)
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>          -425.2            6 862.4 871.4      868.4          882.2        849.9
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             863.2       874.4           888.2      904.4          933.7
#>  icl_groups icl_individual clc_groups clc_individual
#>       868.4          883.9      850.4          852.1
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_groups n_profiles n_group_classes n_profile_covariates
#>             200       20          2               1                    1
#>  n_group_covariates variance_model covariance_model n_parameters log_likelihood
#>                   0        varying         diagonal            6         -425.2
#>    aic bic_groups bic_individual converged boundary extreme_logits n_starts
#>  862.4      868.4          882.2      TRUE    FALSE          FALSE        2
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error
#>      1         -425.2      TRUE          5  <NA>
#>      2         -425.2      TRUE          7  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  group      y
#>      1  3.894
#>      1 -4.047
#>      1  4.971
#>      1  2.616
#>      1  4.654
#>      1  4.512
#>      1 -2.917
#>      1  3.567
#>      1 -4.025
#>      1  3.323
#>    ... 190 more rows.  get_results(x, what = "data")
#> 
#> 18 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
get_results(fit, what = "coefficients")
#>     level   outcome          term parameter     estimate
#> 1 profile profile_1 group_class_1     logit  0.004333305
#> 2 profile profile_1             z     logit -0.826689379
```
