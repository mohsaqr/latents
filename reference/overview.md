# Everything about any object, printed and drawn

`overview()` takes any object the package returns and shows all of it:
the summary, the diagnostic checks and every plot the object can draw.

- A profile, class or transition fit:
  [`report()`](https://pak.dynasite.org/latents/reference/report.md)
  (summary, descriptives, checks, every plot view).

- Any other fit (mixture regression, growth mixture, general transition,
  additive, cross-level, pooled): its summary, its
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  checks where defined, and
  [`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md).

- An enumeration grid: its summary and
  [`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md)
  (criteria, entropy, AvePP, smallest class, BIC changes, trust, tree).

- A comparison, bootstrap test, sensitivity analysis or diagnostics
  result: the object printed and its plot.

## Usage

``` r
overview(x, plots = TRUE)
```

## Arguments

- x:

  Any object returned by the package.

- plots:

  `TRUE` (the default) draws the plots; `FALSE` prints only.

## Value

`x`, invisibly. Called for its printed and drawn output.

An object the package did not return raises `latents_not_a_fit`.

## See also

[`report()`](https://pak.dynasite.org/latents/reference/report.md),
[`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md),
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md).

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
