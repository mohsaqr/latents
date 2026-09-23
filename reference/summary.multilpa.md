# Summarize a fitted multilevel latent profile model

The summary collects the estimates, the effective class counts, the
information criteria and the restart diagnostics of a fit. Every field
is built explicitly, and a field the fit does not carry takes its
documented default rather than being silently absent, so the summary of
a diagonal fit and the summary of a staged full-covariance fit have
exactly the same names. Read the numbers with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), which
returns them as tidy tables; the printed form is a human-facing report.

## Usage

``` r
# S3 method for class 'multilpa'
summary(object, ...)
```

## Arguments

- object:

  An `multilpa` model.

- ...:

  Reserved for compatibility with
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

A `summary_multilpa` object: a named list, never carrying an `NA` name,
with the fit's dimensions (`n_observations`, `n_informative`,
`n_groups`, `n_profiles`, `n_group_classes`), its specification
(`variance_model`, `covariance_model`, `missing`, `continuous`, `fixed`,
`staged`), its estimates (`means`, `variances`, `standard_deviations`,
`covariances`, `response_probabilities`, `profile_probabilities`,
`group_probabilities`), the effective class counts at both levels, the
likelihood, parameter counts and information criteria, and the restart
diagnostics. `covariances` is `NULL` under the diagonal parameterization
and `response_probabilities` is `NULL` when no indicator is categorical;
both keep their names in either case. Use
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
rather than reading the fields.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the tidy tables.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
summary(fit)
#> Multilevel LPA: 2 profiles and 1 group classes
#> Individuals: 120; groups: 12; parameters: 9; converged: TRUE
#> Log likelihood: -323.019699; AIC: 664.039
#> BIC (groups): 668.404; BIC (individuals): 689.127
#> Best likelihood replicated in 2/2 starts (absolute tolerance 0.000324).
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator     mean variance standard_deviation
#>        1   score_a  0.31652   0.7483             0.8650
#>        1   score_b -0.05305   0.9150             0.9565
#>        2   score_a -1.11302   0.1021             0.3196
#>        2   score_b  0.79087   0.2858             0.5346
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile indicator indicator_2 covariance
#>        1   score_a     score_a     0.7483
#>        1   score_b     score_a     0.0000
#>        1   score_a     score_b     0.0000
#>        1   score_b     score_b     0.9150
#>        2   score_a     score_a     0.1021
#>        2   score_b     score_a     0.0000
#>        2   score_a     score_b     0.0000
#>        2   score_b     score_b     0.2858
#> 
#> -- profile_probabilities -------------------------------------------
#>  group_class profile probability group_class_probability
#>            1       1      0.8864                       1
#>            1       2      0.1136                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1          106.37               0.8864
#>  individuals     2           13.63               0.1136
#>       groups     1           12.00               1.0000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1    1.0000  TRUE
#>    2     1       1    0.1984 FALSE
#>    3     1       1    0.5916  TRUE
#>    4     1       1    0.9846  TRUE
#>    5     1       1    0.2976 FALSE
#>    6     1       1    0.2959 FALSE
#>    7     1       1    1.0000  TRUE
#>    8     1       1    0.9980  TRUE
#>    9     1       1    1.0000  TRUE
#>   10     1       1    1.0000  TRUE
#>    ... 230 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         10         -29.64           1         1  TRUE
#>      2         10         -30.05           1         1  TRUE
#>      3         10         -23.78           1         1  TRUE
#>      4         10         -24.45           1         1  TRUE
#>      5         10         -24.16           1         1  TRUE
#>      6         10         -27.79           1         1  TRUE
#>      7         10         -29.24           1         1  TRUE
#>      8         10         -30.43           1         1  TRUE
#>      9         10         -28.55           1         1  TRUE
#>     10         10         -26.44           1         1  TRUE
#>    ... 2 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  school score_a  score_b profile group_class uncertainty posterior_profile_1
#>       1  2.2872 -1.55472       1           1   0.000e+00              1.0000
#>       1 -1.1968  1.56989       2           1   1.984e-01              0.1984
#>       1 -0.6943  0.68845       1           1   4.084e-01              0.5916
#>       1 -0.4123 -0.17760       1           1   1.540e-02              0.9846
#>       1 -0.9707  0.72920       2           1   2.976e-01              0.2976
#>       1 -0.9473  1.53325       2           1   2.959e-01              0.2959
#>       1  0.7481  0.50658       1           1   3.132e-08              1.0000
#>       1 -0.1170  0.03333       1           1   2.009e-03              0.9980
#>       1  0.1527 -1.46755       1           1   9.884e-08              1.0000
#>       1  2.1900  1.01916       1           1   0.000e+00              1.0000
#>  posterior_profile_2
#>            4.963e-28
#>            8.016e-01
#>            4.084e-01
#>            1.540e-02
#>            7.024e-01
#>            7.041e-01
#>            3.132e-08
#>            2.009e-03
#>            9.884e-08
#>            7.100e-23
#>    ... 110 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     104           0.8667      106.37               0.8864
#>  individuals     2      16           0.1333       13.63               0.1136
#>       groups     1      12           1.0000       12.00               1.0000
#>  average_posterior odds_correct_classification
#>             0.9755                       5.112
#>             0.6929                      17.609
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        104           0.97555
#>  individuals              1     2        104           0.02445
#>  individuals              2     1         16           0.30710
#>  individuals              2     2         16           0.69290
#>       groups              1     1         12           1.00000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.95381
#>  individuals          1              2     0.04619
#>  individuals          2              1     0.18658
#>  individuals          2              2     0.81342
#>       groups          1              1     1.00000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.06021
#>  individuals    1              1     2 -0.06021
#>  individuals    2              2     1 -0.24318
#>  individuals    2              2     2  1.24318
#>  individuals    3              1     1  1.06021
#>  individuals    3              1     2 -0.06021
#>  individuals    4              1     1  1.06021
#>  individuals    4              1     2 -0.06021
#>  individuals    5              2     1 -0.24318
#>  individuals    5              2     2  1.24318
#>    ... 242 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     120       15.76           0.8105
#>       groups         1      12        0.00               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2     score_a     score_b gaussian 0.222599        0 0.222599
#>  profile_1     score_a     score_b gaussian 0.008819        0 0.008819
#>  effective_n statistic df p_value p_adjusted
#>        13.63   0.73809 NA  0.4605     0.4605
#>       106.37   0.08967 NA  0.9286     0.9286
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters aic kic bic_groups bic_individual sabic_groups
#>            -323            9 664 676      668.4          689.1        641.2
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             660.7       677.4           698.1      717.8          790.7
#>  icl_groups icl_individual clc_groups clc_individual
#>       668.4          720.7        646          677.6
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             120           120       12          2               1      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                   VVI            9                             9           -323
#>  aic bic_groups bic_individual converged iterations boundary small_classes
#>  664      668.4          689.1      TRUE        130    FALSE         FALSE
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
#>  school score_a  score_b
#>       1  2.2872 -1.55472
#>       1 -1.1968  1.56989
#>       1 -0.6943  0.68845
#>       1 -0.4123 -0.17760
#>       1 -0.9707  0.72920
#>       1 -0.9473  1.53325
#>       1  0.7481  0.50658
#>       1 -0.1170  0.03333
#>       1  0.1527 -1.46755
#>       1  2.1900  1.01916
#>    ... 110 more rows.  get_results(x, what = "data")
#> 
#> 19 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
