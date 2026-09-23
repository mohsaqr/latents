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
set.seed(1)
example_data <- data.frame(group = rep(seq_len(20), each = 10),
                           z = rnorm(200))
example_data$y <- rnorm(200,
  ifelse(runif(200) < plogis(example_data$z), -3, 3))
fit <- multilpa(example_data, "y", "group", n_profiles = 2,
                n_group_classes = 1, profile_covariates = "z",
                n_starts = 2, seed = 1)
nobs(fit)
#> [1] 20
```
