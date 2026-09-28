# Print a covariate LPA fit

Print a covariate LPA fit

## Usage

``` r
# S3 method for class 'multilpa_covariates'
print(x, rows = 20L, ...)
```

## Arguments

- x:

  A covariate LPA fit.

- rows:

  How many rows of the printed table to show before truncating.

- ...:

  Reserved.

## Value

The model, invisibly. Called for the side effect of printing the class
counts, the covariate counts, the log likelihood with the information
criteria, and the convergence diagnostics.

## Examples

``` r
fit <- multilpa(subset(course_engagement, student <= 40),
                c("browse", "lectures", "forum_read"), "student",
                n_profiles = 2, n_group_classes = 1,
                profile_covariates = "previous_grade", n_starts = 2, seed = 1)
print(fit)
#> Multilevel LPA with covariates: 2 profiles, 1 group classes
#> Log likelihood -2099.951166; AIC 4227.902; BIC (groups) 4251.547; converged TRUE
#> 
#>  profile     browse   lectures forum_read    count proportion
#>        1  0.5669052  0.4576051  0.6311310 352.8628  0.6450874
#>        2 -0.8400039 -0.6248737 -0.9162981 194.1372  0.3549126
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
