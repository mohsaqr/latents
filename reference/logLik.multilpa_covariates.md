# Extract a covariate LPA log likelihood

Extract a covariate LPA log likelihood

## Usage

``` r
# S3 method for class 'multilpa_covariates'
logLik(object, ...)
```

## Arguments

- object:

  A covariate LPA fit.

- ...:

  Reserved.

## Value

A `logLik` object carrying the maximized log likelihood, the free
parameter count as `df`, and the number of observed groups as `nobs`, so
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) uses the group-count
BIC.

## Examples

``` r
fit <- multilpa(subset(course_engagement, student <= 40),
                c("browse", "lectures", "forum_read"), "student",
                n_profiles = 2, n_group_classes = 1,
                profile_covariates = "previous_grade", n_starts = 2, seed = 1)
logLik(fit)
#> 'log Lik.' -2099.951 (df=14)
```
