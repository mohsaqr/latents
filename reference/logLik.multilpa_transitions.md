# Log likelihood of a fitted latent transition model

Log likelihood of a fitted latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
logLik(object, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.

- ...:

  Ignored.

## Value

A `logLik` object carrying the maximized observed-data log likelihood,
the free parameter count as `df`, and the number of groups as `nobs`.
Groups are the independent units, because a group's occasions are
dependent by construction in this model.

## Examples

``` r
fit <- lta(subset(course_engagement, student <= 40),
           c("browse", "lectures", "forum_read"), "student",
           n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
logLik(fit)
#> 'log Lik.' -2025.406 (df=15)
```
