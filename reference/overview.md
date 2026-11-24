# Print and plot everything about any object of the package: its summary, its checks and every plot

`overview()` prints and plots everything about any object the package
returns: the summary, the diagnostic checks where the object has them,
and every plot the object can draw. What is shown depends on the object:

- a profile, class or transition fit of
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  or [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) (class
  `multilpa_transitions`): the output of
  [`report()`](https://pak.dynasite.org/latents/reference/report.md),
  that is, the summary, the descriptives, the checks and every plot
  view;

- a fit of another family (a mixture regression, a growth mixture, a
  general transition fit, an additive or a cross-level fit): its
  summary, the checks of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md),
  and the plots of
  [`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md);

- pooled imputations: the summary and the plots of
  [`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md);

- an enumeration grid: its summary and the plots of
  [`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md)
  (criteria, entropy, AvePP, criterion changes, heat map, smallest
  class, trust and tree);

- a comparison, a bootstrap test, a bootstrap table of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  or a sensitivity analysis: the object or its summary printed, and its
  plot;

- a
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  result: the result printed, and the plots of the fit it checked.

## Usage

``` r
overview(x, plots = TRUE)
```

## Arguments

- x:

  Any object returned by the package: a fitted model, pooled
  imputations, an enumeration grid, a comparison, a bootstrap result, a
  sensitivity analysis or a diagnostics result.

- plots:

  Whether to draw the plots:

  - `TRUE` (default): prints, then draws every plot. The plots need the
    ggplot2 package; without it only the printed part is shown;

  - `FALSE`: prints only.

## Value

`x`, invisibly. The function is called for what it prints and draws.

Errors: `latents_not_a_fit` when `x` is not a fit, grid or other result
of the package, for example a data frame; `latents_bad_argument` when
`plots` is not `TRUE` or `FALSE`.

## See also

[`report()`](https://pak.dynasite.org/latents/reference/report.md) for
the fits it covers;
[`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md)
for the plots alone;
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md)
for the checks as a table.

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 2, seed = 1)
overview(fit, plots = FALSE)
#> Model fit
#> 
#> Family          gaussian
#> Nesting         observation
#> Classes         2
#> Group classes   1
#> Observations    900
#> Parameters      7
#> Log likelihood  -3144.01
#> AIC             6302.01
#> BIC             6335.63
#> SABIC           6313.40
#> ICL             7030.42
#> Entropy         0.443
#> Smallest class  43.5%
#> Converged       yes
#> 
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.37  [32.96, 35.78]  <.001
#> Class 1  hours          4.56  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.91  [52.77, 57.05]  <.001
#> Class 2  hours          0.81  [ 0.47,  1.14]  <.001
#> 
#> Classes
#> 
#> Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
#> -------  -----  --------  -------------  --------------  -----------
#> Class 1  0.565    508.10            573           0.789         5.24
#> Class 2  0.435    391.90            327           0.829         6.89
#> 
#> Class membership (log odds, 95% CI)
#> 
#> Class    Term       Estimate         95% CI     p
#> -------  ---------  --------  -------------  ----
#> Class 2  Intercept     -0.26  [-0.54, 0.02]  .068
#> 
#> 
#> Diagnostics: Mixture of gaussian regressions, 2 classes (one class per row)
#> No check needs attention.
#> 
#> Estimation
#>   ok     best start converged                                yes                      
#>   info   best log likelihood reached by                      2 of 2 starts            
#>   ok     starts failed or unfinished                         0 of 2 (1 screened)      
#> 
#> Class size
#>   info   smallest class                         individuals  43.5% (n = 392, class 2) 
#>   info   modal vs estimated proportions         individuals  largest gap 0.072        
#>   ok     classes with no assigned unit          individuals  0 of 2                   
#> 
#> Classification
#>   info   relative entropy                       individuals  0.443                    
#>   info   lowest average posterior               individuals  0.789 (class 1)          
#>   info   lowest odds of correct classification  individuals  2.9 (class 1)            
#> 
#> Local independence
#>   n/a    bivariate residuals                                                          
#>          not defined for this model family
#> 
#> Information
#>   info   log likelihood, parameters, BIC                     -3144.01, 7, 6335.63     
#> 
#> Every check with its reference: get_results(x, "checks"). Other tables: get_results(x, what = ).
```
