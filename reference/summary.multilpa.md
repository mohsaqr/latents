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
#> Log likelihood: -323.019688; AIC: 664.039
#> BIC (groups): 668.404; BIC (individuals): 689.127
#> Best likelihood replicated in 1/2 starts (absolute tolerance 0.000324).
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator     mean variance standard_deviation mean_standard_error
#>        1   score_a  0.31637   0.7483             0.8650              0.1034
#>        1   score_b -0.05292   0.9150             0.9565              0.1003
#>        2   score_a -1.11360   0.1019             0.3192              0.1775
#>        2   score_b  0.79092   0.2858             0.5346              0.2233
#>  variance_standard_error
#>                  0.11790
#>                  0.12865
#>                  0.08122
#>                  0.14476
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
#>        2   score_a     score_a     0.1019
#>        2   score_b     score_a     0.0000
#>        2   score_a     score_b     0.0000
#>        2   score_b     score_b     0.2858
#> 
#> -- profile_probabilities -------------------------------------------
#>  group_class profile probability group_class_probability
#>            1       1      0.8865                       1
#>            1       2      0.1135                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1          106.39               0.8865
#>  individuals     2           13.61               0.1135
#>       groups     1           12.00               1.0000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1    1.0000  TRUE
#>    2     1       1    0.1984 FALSE
#>    3     1       1    0.5928  TRUE
#>    4     1       1    0.9848  TRUE
#>    5     1       1    0.2979 FALSE
#>    6     1       1    0.2963 FALSE
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
#>       1 -0.6943  0.68845       1           1   4.072e-01              0.5928
#>       1 -0.4123 -0.17760       1           1   1.524e-02              0.9848
#>       1 -0.9707  0.72920       2           1   2.979e-01              0.2979
#>       1 -0.9473  1.53325       2           1   2.963e-01              0.2963
#>       1  0.7481  0.50658       1           1   2.963e-08              1.0000
#>       1 -0.1170  0.03333       1           1   1.972e-03              0.9980
#>       1  0.1527 -1.46755       1           1   9.602e-08              1.0000
#>       1  2.1900  1.01916       1           1   0.000e+00              1.0000
#>  posterior_profile_2
#>            4.187e-28
#>            8.016e-01
#>            4.072e-01
#>            1.524e-02
#>            7.021e-01
#>            7.037e-01
#>            2.963e-08
#>            1.972e-03
#>            9.602e-08
#>            6.050e-23
#>    ... 110 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     104           0.8667      106.39               0.8865
#>  individuals     2      16           0.1333       13.61               0.1135
#>       groups     1      12           1.0000       12.00               1.0000
#>  average_posterior odds_correct_classification
#>             0.9757                       5.128
#>             0.6927                      17.611
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        104           0.97565
#>  individuals              1     2        104           0.02435
#>  individuals              2     1         16           0.30734
#>  individuals              2     2         16           0.69266
#>       groups              1     1         12           1.00000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.95378
#>  individuals          1              2     0.04622
#>  individuals          2              1     0.18599
#>  individuals          2              2     0.81401
#>       groups          1              1     1.00000
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
#>  individuals         2     120       15.74           0.8107
#>       groups         1      12        0.00               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2     score_a     score_b gaussian 0.223117        0 0.223117
#>  profile_1     score_a     score_b gaussian 0.008692        0 0.008692
#>  effective_n statistic df p_value p_adjusted
#>        13.61   0.73936 NA  0.4597     0.4597
#>       106.39   0.08838 NA  0.9296     0.9296
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters aic kic bic_groups bic_individual sabic_groups
#>            -323            9 664 676      668.4          689.1        641.2
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             660.7       677.4           698.1      717.8          790.7
#>  icl_groups icl_individual clc_groups clc_individual
#>       668.4          720.6        646          677.5
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             120           120       12          2               1      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                   VVI            9                             9           -323
#>  aic bic_groups bic_individual converged iterations boundary small_classes
#>  664      668.4          689.1      TRUE         30    FALSE         FALSE
#>  best_start n_best_replicated
#>           2                 1
#> 
#> -- stages ----------------------------------------------------------
#>  stage group_classes fixed log_likelihood parameters
#>  joint             1  <NA>           -323          9
#>  parameters_with_measurement converged
#>                            9      TRUE
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error boundary
#>      1         -324.9      TRUE         48  <NA>    FALSE
#>      2         -323.0      TRUE         30  <NA>    FALSE
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
