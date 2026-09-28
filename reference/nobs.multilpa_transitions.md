# Number of independent units in a fitted latent transition model

Number of independent units in a fitted latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.

- ...:

  Ignored.

## Value

A single integer: the number of observed groups. Occasions within a
group are dependent by construction, so they are not independent
observations.

## Examples

``` r
fit <- lta(subset(course_engagement, student <= 40),
           c("browse", "lectures", "forum_read"), "student",
           n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
nobs(fit)
#> [1] 40
```
