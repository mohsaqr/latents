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
#> Log likelihood: -323.019699; AIC: 664.039
#> BIC (groups): 668.404; BIC (individuals): 689.127
#> Best likelihood replicated in 2/2 starts (absolute tolerance 0.000324).
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator   mean variance standard_deviation
#>        1   score_a  0.317    0.748              0.865
#>        1   score_b -0.053    0.915              0.957
#>        2   score_a -1.113    0.102              0.320
#>        2   score_b  0.791    0.286              0.535
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
#>            1       1       0.886                       1
#>            1       2       0.114                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1           106.4                0.886
#>  individuals     2            13.6                0.114
#>       groups     1            12.0                1.000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1     1.000  TRUE
#>    2     1       1     0.198 FALSE
#>    3     1       1     0.592  TRUE
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
#>       1  -0.694  0.6884       1           1    4.08e-01               0.592
#>       1  -0.412 -0.1776       1           1    1.54e-02               0.985
#>       1  -0.971  0.7292       2           1    2.98e-01               0.298
#>       1  -0.947  1.5333       2           1    2.96e-01               0.296
#>       1   0.748  0.5066       1           1    3.13e-08               1.000
#>       1  -0.117  0.0333       1           1    2.01e-03               0.998
#>       1   0.153 -1.4676       1           1    9.88e-08               1.000
#>       1   2.190  1.0192       1           1    0.00e+00               1.000
#>  posterior_profile_2
#>             4.96e-28
#>             8.02e-01
#>             4.08e-01
#>             1.54e-02
#>             7.02e-01
#>             7.04e-01
#>             3.13e-08
#>             2.01e-03
#>             9.88e-08
#>             7.10e-23
#>    ... 110 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     104            0.867       106.4                0.886
#>  individuals     2      16            0.133        13.6                0.114
#>       groups     1      12            1.000        12.0                1.000
#>  average_posterior odds_correct_classification
#>              0.976                        5.11
#>              0.693                       17.61
#>              1.000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        104            0.9755
#>  individuals              1     2        104            0.0245
#>  individuals              2     1         16            0.3071
#>  individuals              2     2         16            0.6929
#>       groups              1     1         12            1.0000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1      0.9538
#>  individuals          1              2      0.0462
#>  individuals          2              1      0.1866
#>  individuals          2              2      0.8134
#>       groups          1              1      1.0000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class  weight
#>  individuals    1              1     1  1.0602
#>  individuals    1              1     2 -0.0602
#>  individuals    2              2     1 -0.2432
#>  individuals    2              2     2  1.2432
#>  individuals    3              1     1  1.0602
#>  individuals    3              1     2 -0.0602
#>  individuals    4              1     1  1.0602
#>  individuals    4              1     2 -0.0602
#>  individuals    5              2     1 -0.2432
#>  individuals    5              2     2  1.2432
#>    ... 242 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     120        15.8             0.81
#>       groups         1      12         0.0               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2     score_a     score_b gaussian  0.22260        0  0.22260
#>  profile_1     score_a     score_b gaussian  0.00882        0  0.00882
#>  effective_n statistic df p_value p_adjusted
#>         13.6    0.7381 NA   0.460      0.460
#>        106.4    0.0897 NA   0.929      0.929
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
#>  664        668            689      TRUE        130    FALSE         FALSE
#>  best_start n_best_replicated
#>           2                 2
#> 
#> -- stages ----------------------------------------------------------
#>  stage group_classes fixed log_likelihood parameters
#>  joint             1  <NA>           -323          9
#>  parameters_with_measurement converged
#>                            9      TRUE
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error boundary
#>      1           -323      TRUE        114  <NA>    FALSE
#>      2           -323      TRUE        130  <NA>    FALSE
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
