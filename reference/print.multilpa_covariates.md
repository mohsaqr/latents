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
set.seed(1)
example_data <- data.frame(group = rep(seq_len(20), each = 10),
                           z = rnorm(200))
example_data$y <- rnorm(200,
  ifelse(runif(200) < plogis(example_data$z), -3, 3))
fit <- multilpa(example_data, "y", "group", n_profiles = 2,
                n_group_classes = 1, profile_covariates = "z",
                n_starts = 2, seed = 1)
print(fit)
#> Multilevel LPA with covariates: 2 profiles, 1 group classes
#> Log likelihood -425.206918; AIC 862.414; BIC (groups) 868.388; converged TRUE
#> 
#>  profile         y    count proportion
#>        1  3.008702  99.1541  0.4957705
#>        2 -3.005361 100.8459  0.5042295
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
