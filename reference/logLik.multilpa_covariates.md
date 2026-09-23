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
set.seed(1)
example_data <- data.frame(group = rep(seq_len(20), each = 10),
                           z = rnorm(200))
example_data$y <- rnorm(200,
  ifelse(runif(200) < plogis(example_data$z), -3, 3))
fit <- multilpa(example_data, "y", "group", n_profiles = 2,
                n_group_classes = 1, profile_covariates = "z",
                n_starts = 2, seed = 1)
logLik(fit)
#> 'log Lik.' -425.2069 (df=6)
```
