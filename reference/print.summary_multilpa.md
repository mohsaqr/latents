# Print a multilevel LPA summary

A human-facing report of the fit. The same content is available as tidy
tables from
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).

## Usage

``` r
# S3 method for class 'summary_multilpa'
print(x, digits = 4L, rows = 10L, ...)
```

## Arguments

- x:

  A `summary_multilpa` object.

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

The summary, invisibly. Called for the side effect of printing the
estimates block by block, the effective class memberships at both
levels, the likelihood and information criteria, any convergence or
boundary warnings, and the restart diagnostics.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the same content as data.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
print(summary(fit), digits = 3)
#> Multilevel LPA: 2 profiles and 1 group classes
#> Individuals: 120; groups: 12; parameters: 9; converged: TRUE
#> Log likelihood: -323.019688; AIC: 664.039
#> BIC (groups): 668.404; BIC (individuals): 689.127
#> Best likelihood replicated in 1/2 starts (absolute tolerance 0.000324).
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator    mean variance standard_deviation mean_standard_error
#>        1   score_a  0.3164    0.748              0.865               0.103
#>        1   score_b -0.0529    0.915              0.957               0.100
#>        2   score_a -1.1136    0.102              0.319               0.177
#>        2   score_b  0.7909    0.286              0.535               0.223
#>  variance_standard_error
#>                   0.1179
#>                   0.1286
#>                   0.0812
#>                   0.1448
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile indicator indicator_2 covariance
#>        1   score_a     score_a      0.748
#>        1   score_b     score_a      0.000
#>        1   score_a     score_b      0.000
#>        1   score_b     score_b      0.915
#>        2   score_a     score_a      0.102
#>        2   score_b     score_a      0.000
#>        2   score_a     score_b      0.000
#>        2   score_b     score_b      0.286
#> 
#> -- profile_probabilities -------------------------------------------
#>  group_class profile probability group_class_probability
#>            1       1       0.887                       1
#>            1       2       0.113                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1           106.4                0.887
#>  individuals     2            13.6                0.113
#>       groups     1            12.0                1.000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1     1.000  TRUE
#>    2     1       1     0.198 FALSE
#>    3     1       1     0.593  TRUE
#>    4     1       1     0.985  TRUE
#>    5     1       1     0.298 FALSE
#>    6     1       1     0.296 FALSE
#>    7     1       1     1.000  TRUE
#>    8     1       1     0.998  TRUE
#>    9     1       1     1.000  TRUE
#>   10     1       1     1.000  TRUE
#>    ... 230 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         10          -29.6           1         1  TRUE
#>      2         10          -30.1           1         1  TRUE
#>      3         10          -23.8           1         1  TRUE
#>      4         10          -24.5           1         1  TRUE
#>      5         10          -24.2           1         1  TRUE
#>      6         10          -27.8           1         1  TRUE
#>      7         10          -29.2           1         1  TRUE
#>      8         10          -30.4           1         1  TRUE
#>      9         10          -28.5           1         1  TRUE
#>     10         10          -26.4           1         1  TRUE
#>    ... 2 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  school score_a score_b profile group_class uncertainty posterior_profile_1
#>       1   2.287 -1.5547       1           1    0.00e+00               1.000
#>       1  -1.197  1.5699       2           1    1.98e-01               0.198
#>       1  -0.694  0.6884       1           1    4.07e-01               0.593
#>       1  -0.412 -0.1776       1           1    1.52e-02               0.985
#>       1  -0.971  0.7292       2           1    2.98e-01               0.298
#>       1  -0.947  1.5333       2           1    2.96e-01               0.296
#>       1   0.748  0.5066       1           1    2.96e-08               1.000
#>       1  -0.117  0.0333       1           1    1.97e-03               0.998
#>       1   0.153 -1.4676       1           1    9.60e-08               1.000
#>       1   2.190  1.0192       1           1    0.00e+00               1.000
#>  posterior_profile_2
#>             4.19e-28
#>             8.02e-01
#>             4.07e-01
#>             1.52e-02
#>             7.02e-01
#>             7.04e-01
#>             2.96e-08
#>             1.97e-03
#>             9.60e-08
#>             6.05e-23
#>    ... 110 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     104            0.867       106.4                0.887
#>  individuals     2      16            0.133        13.6                0.113
#>       groups     1      12            1.000        12.0                1.000
#>  average_posterior odds_correct_classification
#>              0.976                        5.13
#>              0.693                       17.61
#>              1.000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        104            0.9757
#>  individuals              1     2        104            0.0243
#>  individuals              2     1         16            0.3073
#>  individuals              2     2         16            0.6927
#>       groups              1     1         12            1.0000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1      0.9538
#>  individuals          1              2      0.0462
#>  individuals          2              1      0.1860
#>  individuals          2              2      0.8140
#>       groups          1              1      1.0000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class  weight
#>  individuals    1              1     1  1.0602
#>  individuals    1              1     2 -0.0602
#>  individuals    2              2     1 -0.2422
#>  individuals    2              2     2  1.2422
#>  individuals    3              1     1  1.0602
#>  individuals    3              1     2 -0.0602
#>  individuals    4              1     1  1.0602
#>  individuals    4              1     2 -0.0602
#>  individuals    5              2     1 -0.2422
#>  individuals    5              2     2  1.2422
#>    ... 242 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     120        15.7            0.811
#>       groups         1      12         0.0               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2     score_a     score_b gaussian  0.22312        0  0.22312
#>  profile_1     score_a     score_b gaussian  0.00869        0  0.00869
#>  effective_n statistic df p_value p_adjusted
#>         13.6    0.7394 NA    0.46       0.46
#>        106.4    0.0884 NA    0.93       0.93
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters aic kic bic_groups bic_individual sabic_groups
#>            -323            9 664 676        668            689          641
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>               661         677             698        718            791
#>  icl_groups icl_individual clc_groups clc_individual
#>         668            721        646            678
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             120           120       12          2               1      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                   VVI            9                             9           -323
#>  aic bic_groups bic_individual converged iterations boundary small_classes
#>  664        668            689      TRUE         30    FALSE         FALSE
#>  best_start n_best_replicated weights
#>           2                 1    <NA>
#> 
#> -- stages ----------------------------------------------------------
#>  stage group_classes fixed log_likelihood parameters
#>  joint             1  <NA>           -323          9
#>  parameters_with_measurement converged
#>                            9      TRUE
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error boundary
#>      1           -325      TRUE         48  <NA>    FALSE
#>      2           -323      TRUE         30  <NA>    FALSE
#> 
#> -- data ------------------------------------------------------------
#>  school score_a score_b
#>       1   2.287 -1.5547
#>       1  -1.197  1.5699
#>       1  -0.694  0.6884
#>       1  -0.412 -0.1776
#>       1  -0.971  0.7292
#>       1  -0.947  1.5333
#>       1   0.748  0.5066
#>       1  -0.117  0.0333
#>       1   0.153 -1.4676
#>       1   2.190  1.0192
#>    ... 110 more rows.  get_results(x, what = "data")
#> 
#> 19 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
