# Plot a pooled multiply imputed fit

One row per parameter: each imputation's aligned estimate as an open
circle and the pooled estimate as a filled circle with its interval, so
the between-imputation spread is visible next to the total uncertainty.

## Usage

``` r
# S3 method for class 'latents_pooled'
plot(x, level = c("profile", "group", "measurement"), main = NULL, ...)
```

## Arguments

- x:

  A `latents_pooled` object.

- level:

  Which parameter level to show: `"profile"` (the profile membership
  coefficients or probabilities), `"group"` or `"measurement"`.

- main:

  Plot title; `NULL` for the default.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the side effect of drawing.

## Examples

``` r
set.seed(1)
completed <- lapply(1:2, function(i) {
  within(subset(course_engagement, student <= 40),
         previous_grade <- previous_grade + rnorm(length(previous_grade), sd = 0.1))
})
pooled <- pool_imputations(completed, c("browse", "lectures", "forum_read"),
                           "student", n_profiles = 2, n_group_classes = 1,
                           profile_covariates = "previous_grade",
                           n_starts = 2, seed = 1)
plot(pooled)
```
