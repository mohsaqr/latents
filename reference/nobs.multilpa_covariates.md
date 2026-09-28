# Count independent groups in a covariate LPA fit

Count independent groups in a covariate LPA fit

## Usage

``` r
# S3 method for class 'multilpa_covariates'
nobs(object, ...)
```

## Arguments

- object:

  A covariate LPA fit.

- ...:

  Reserved.

## Value

A single integer: the number of observed groups, which are the
independent units of this likelihood.

## Examples

``` r
fit <- multilpa(subset(course_engagement, student <= 40),
                c("browse", "lectures", "forum_read"), "student",
                n_profiles = 2, n_group_classes = 1,
                profile_covariates = "previous_grade", n_starts = 2, seed = 1)
nobs(fit)
#> [1] 40
```
