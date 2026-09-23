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
set.seed(1)
example_data <- data.frame(group = rep(seq_len(20), each = 10),
                           z = rnorm(200))
example_data$y <- rnorm(200,
  ifelse(runif(200) < plogis(example_data$z), -3, 3))
fit <- multilpa(example_data, "y", "group", n_profiles = 2,
                n_group_classes = 1, profile_covariates = "z",
                n_starts = 2, seed = 1)
print(summary(fit), digits = 3)
#> Multilevel LPA with covariates: 2 profiles, 1 group classes
#> Individuals: 200; groups: 20; parameters: 6; converged: TRUE
#> 1 profile covariate(s); 0 group covariate(s)
#> Log likelihood: -425.206918; AIC: 862.414
#> BIC (groups): 868.388; BIC (individuals): 882.204
#> Membership coefficients carry no standard errors here; parameter_inference() has them.
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator  mean variance standard_deviation
#>        1         y  3.01     1.24               1.11
#>        2         y -3.01     1.11               1.06
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile indicator indicator_2 covariance
#>        1         y           y       1.24
#>        2         y           y       1.11
#> 
#> -- coefficients ----------------------------------------------------
#>    level   outcome          term parameter estimate
#>  profile profile_1 group_class_1     logit  0.00433
#>  profile profile_1             z     logit -0.82669
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1            99.2                0.496
#>  individuals     2           100.8                0.504
#>       groups     1            20.0                1.000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1  1.00e+00  TRUE
#>    2     1       1  2.36e-09 FALSE
#>    3     1       1  1.00e+00  TRUE
#>    4     1       1  1.00e+00  TRUE
#>    5     1       1  1.00e+00  TRUE
#>    6     1       1  1.00e+00  TRUE
#>    7     1       1  4.29e-07 FALSE
#>    8     1       1  1.00e+00  TRUE
#>    9     1       1  1.90e-09 FALSE
#>   10     1       1  1.00e+00  TRUE
#>    ... 390 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         10          -22.1           1         1  TRUE
#>      2         10          -18.9           1         1  TRUE
#>      3         10          -19.5           1         1  TRUE
#>      4         10          -20.7           1         1  TRUE
#>      5         10          -22.8           1         1  TRUE
#>      6         10          -20.4           1         1  TRUE
#>      7         10          -20.4           1         1  TRUE
#>      8         10          -18.8           1         1  TRUE
#>      9         10          -22.0           1         1  TRUE
#>     10         10          -19.2           1         1  TRUE
#>    ... 10 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  group     y profile group_class uncertainty posterior_profile_1
#>      1  3.89       1           1    4.51e-10            1.00e+00
#>      1 -4.05       2           1    2.36e-09            2.36e-09
#>      1  4.97       1           1    9.86e-13            1.00e+00
#>      1  2.62       1           1    2.88e-06            1.00e+00
#>      1  4.65       1           1    1.51e-11            1.00e+00
#>      1  4.51       1           1    1.28e-11            1.00e+00
#>      1 -2.92       2           1    4.29e-07            4.29e-07
#>      1  3.57       1           1    8.30e-09            1.00e+00
#>      1 -4.02       2           1    1.90e-09            1.90e-09
#>      1  3.32       1           1    1.32e-08            1.00e+00
#>  posterior_profile_2
#>             4.51e-10
#>             1.00e+00
#>             9.86e-13
#>             2.88e-06
#>             1.51e-11
#>             1.28e-11
#>             1.00e+00
#>             8.30e-09
#>             1.00e+00
#>             1.32e-08
#>    ... 190 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1      99            0.495        99.2                0.496
#>  individuals     2     101            0.505       100.8                0.504
#>       groups     1      20            1.000        20.0                1.000
#>  average_posterior odds_correct_classification
#>              1.000                        5148
#>              0.998                         571
#>              1.000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1         99          0.999802
#>  individuals              1     2         99          0.000198
#>  individuals              2     1        101          0.001719
#>  individuals              2     2        101          0.998281
#>       groups              1     1         20          1.000000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1    0.998249
#>  individuals          1              2    0.001751
#>  individuals          2              1    0.000194
#>  individuals          2              2    0.999806
#>       groups          1              1    1.000000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class    weight
#>  individuals    1              1     1  1.001755
#>  individuals    1              1     2 -0.001755
#>  individuals    2              2     1 -0.000194
#>  individuals    2              2     2  1.000194
#>  individuals    3              1     1  1.001755
#>  individuals    3              1     2 -0.001755
#>  individuals    4              1     1  1.001755
#>  individuals    4              1     2 -0.001755
#>  individuals    5              1     1  1.001755
#>  individuals    5              1     2 -0.001755
#>    ... 410 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     200       0.846            0.994
#>       groups         1      20       0.000               NA
#> 
#> -- residuals -------------------------------------------------------
#>    (no rows)
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters aic kic bic_groups bic_individual sabic_groups
#>            -425            6 862 871        868            882          850
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>               863         874             888        904            934
#>  icl_groups icl_individual clc_groups clc_individual
#>         868            884        850            852
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_groups n_profiles n_group_classes n_profile_covariates
#>             200       20          2               1                    1
#>  n_group_covariates variance_model covariance_model n_parameters log_likelihood
#>                   0        varying         diagonal            6           -425
#>  aic bic_groups bic_individual converged boundary extreme_logits n_starts
#>  862        868            882      TRUE    FALSE          FALSE        2
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error
#>      1           -425      TRUE          5  <NA>
#>      2           -425      TRUE          7  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  group     y
#>      1  3.89
#>      1 -4.05
#>      1  4.97
#>      1  2.62
#>      1  4.65
#>      1  4.51
#>      1 -2.92
#>      1  3.57
#>      1 -4.02
#>      1  3.32
#>    ... 190 more rows.  get_results(x, what = "data")
#> 
#> 18 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
