# Print latents results

Print methods for the objects the package returns: fitted models,
enumerations, bootstrap tests, covariate and transition models, pooled
imputations, their summaries, result tables and plot lists. Each prints
a compact overview;
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
returns the tables themselves as data frames.

## Usage

``` r
# S3 method for class 'multilpa_enumeration'
print(x, ...)

# S3 method for class 'multilpa_start'
print(x, ...)

# S3 method for class 'summary_multilpa_additive'
print(x, digits = 4L, ...)

# S3 method for class 'multilpa_additive'
print(x, ...)

# S3 method for class 'multilpa_covariates'
print(x, rows = 20L, ...)

# S3 method for class 'summary_multilpa_covariates'
print(x, digits = 4L, rows = 10L, ...)

# S3 method for class 'summary_multilpa_enumeration'
print(x, digits = 4L, rows = 10L, ...)

# S3 method for class 'multilpa_bootstrap_lrt'
print(x, ...)

# S3 method for class 'summary_multilpa_bootstrap_lrt'
print(x, digits = 4L, rows = 10L, ...)

# S3 method for class 'latents_growth_mixture'
print(x, digits = 4L, ...)

# S3 method for class 'summary_latents_growth_mixture'
print(x, digits = 4L, ...)

# S3 method for class 'latents_table'
print(x, n = 20L, ...)

# S3 method for class 'multilpa'
print(x, rows = 20L, ...)

# S3 method for class 'summary_multilpa'
print(x, digits = 4L, rows = 10L, ...)

# S3 method for class 'latents_mixture_regression'
print(x, digits = 4L, ...)

# S3 method for class 'summary_latents_mixture_regression'
print(x, digits = 4L, ...)

# S3 method for class 'latents_plots'
print(x, ...)

# S3 method for class 'latents_pooled'
print(x, digits = 4L, rows = 20L, ...)

# S3 method for class 'multilpa_transitions'
print(x, rows = 20L, ...)

# S3 method for class 'summary_multilpa_transitions'
print(x, digits = 4L, rows = 10L, ...)
```

## Arguments

- x:

  An object returned by a latents function, or the
  [`summary()`](https://rdrr.io/r/base/summary.html) of one.

- ...:

  Passed to the underlying data frame printing, or ignored.

- digits:

  Number of significant digits printed.

- rows:

  How many rows of each table to print. A longer table is shown to that
  depth, with its remaining row count and the
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  call that returns it whole.

- n:

  The most rows of a result table printed.

## Value

`x`, invisibly. Called for the side effect of printing.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the tables,
[latents-summary](https://pak.dynasite.org/latents/reference/latents-summary.md)
and
[latents-as-data-frame](https://pak.dynasite.org/latents/reference/latents-as-data-frame.md).

## Examples

``` r
fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 2, seed = 1)
print(fit, rows = 3)
#> Two-level latent profile analysis: 2 profiles, 2 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -5390.952704 | AIC: 10811.905 | BIC (groups): 10851.857
#> Converged: TRUE | iterations: 12 | best start: 1/2
#> 
#>  profile     browse   lectures forum_read    count proportion
#>        1  0.5388561  0.4194852  0.6130248 842.6325  0.5925686
#>        2 -0.7837817 -0.6102031 -0.8914110 579.3675  0.4074314
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
