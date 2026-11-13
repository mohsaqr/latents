# Every diagnostic check of a fit

Each check
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
runs, with its status, criterion and reference: estimation, class sizes,
classification, local independence and information. For a fit, the table
is the one
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
computes with its defaults; for a grid, the tables of every candidate
are stacked with a first column `candidate` (a candidate whose
estimation failed has no row); for a diagnostics result, its own table.

## Usage

``` r
get_checks(x)
```

## Arguments

- x:

  A fitted model of any family, an enumeration grid, or the result of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).

## Value

A data frame, one row per check (and candidate):
`domain, check, level, value, display, criterion, status, detail, reference`.
Any other object raises `latents_not_a_fit`.

## See also

[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
for the checks and their criteria;
[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md),
[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md),
[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md),
[`get_replication()`](https://pak.dynasite.org/latents/reference/get_replication.md),
[`get_residuals()`](https://pak.dynasite.org/latents/reference/get_residuals.md)
for the tables behind them;
[`overview()`](https://pak.dynasite.org/latents/reference/overview.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
get_checks(fit)
#>                domain                                 check       level
#> 1          Estimation                  best start converged       model
#> 2          Estimation        best log likelihood reached by       model
#> 3          Estimation           starts failed or unfinished       model
#> 4          Estimation                  estimates at a bound       model
#> 5          Class size                        smallest class individuals
#> 6          Class size        modal vs estimated proportions individuals
#> 7          Class size         classes with no assigned unit individuals
#> 8      Classification                      relative entropy individuals
#> 9      Classification              lowest average posterior individuals
#> 10     Classification lowest odds of correct classification individuals
#> 11 Local independence        pairs dependent within a class       model
#> 12        Information       log likelihood, parameters, BIC       model
#>            value                 display            criterion status
#> 1   1.000000e+00                     yes                          ok
#> 2   2.000000e+00           2 of 2 starts                        info
#> 3   0.000000e+00                  0 of 2                          ok
#> 4   0.000000e+00                      no                          ok
#> 5   3.565167e-01 35.7% (n = 97, class 2)                        info
#> 6   1.009108e-04       largest gap 0.000                        info
#> 7   0.000000e+00                  0 of 2                          ok
#> 8   9.988000e-01                   0.999                        info
#> 9   9.996149e-01         1.000 (class 2)                        info
#> 10  4.684632e+03        4684.6 (class 2)                        info
#> 11  2.000000e+00     2 of 2 tested pairs BH-adjusted p < 0.05  check
#> 12 -1.147806e+03    -1147.81, 9, 2346.06                        info
#>                                                                                      detail
#> 1                                                                                          
#> 2                                                                                          
#> 3                                                                                          
#> 4                                                                                          
#> 5                                                                                          
#> 6                                                                                          
#> 7                                                                                          
#> 8                                                                                          
#> 9                                                                                          
#> 10                                                                                         
#> 11 smallest adjusted p: eruptions x waiting in profile_1, residual 0.374, adjusted p <1e-04
#> 12                                   meaningful only against other models: compare_models()
#>                              reference
#> 1                                     
#> 2                         Masyn (2013)
#> 3                                     
#> 4                                     
#> 5                                     
#> 6                                     
#> 7                                     
#> 8  Clark & Muthen (2009); Masyn (2013)
#> 9                         Nagin (2005)
#> 10                        Nagin (2005)
#> 11           Vermunt & Magidson (2005)
#> 12                                    
```
