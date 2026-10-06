# Plot transition probabilities of a general transition fit

One panel per group class (and occasion, when transitions vary): a tile
per origin and destination, labelled with its probability.

## Usage

``` r
# S3 method for class 'multilpa_lta'
plot(x, main = NULL, subtitle = NULL, ...)
```

## Arguments

- x:

  A `multilpa_lta` fit.

- main, subtitle:

  Optional title and subtitle.

- ...:

  Unused.

## Value

A ggplot object. Raises `latents_missing_package` without ggplot2.

## Examples

``` r
fit <- lta(course_engagement, c("browse", "lectures"), "student",
           n_profiles = 2, time = "sequence",
           transition_covariates = "previous_grade", n_starts = 2, seed = 1)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(fit)
```
