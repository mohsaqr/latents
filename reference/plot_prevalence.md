# Profile prevalence within each group class

The probability of each profile within each group class, the quantity a
two-level model adds to a one-level one. On a cross-level fit it draws
the profile composition of each group class and on pooled imputations
the pooled profile probabilities beside each imputation's.

## Usage

``` r
plot_prevalence(x, labels = TRUE, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A fitted model of any family.

- labels:

  `TRUE` (default) labels each profile directly; `FALSE` uses a legend.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object; a fit with one group class is drawn too. A covariate
model has no view and raises `latents_nothing_to_plot`.

## See also

[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md)
for the class sizes at each level.

## Examples

``` r
fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 2, seed = 1)
plot_prevalence(fit)
```
