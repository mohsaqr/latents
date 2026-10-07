# Plot information criteria across an enumeration grid

Draws information criteria against the number of profiles, one panel per
criterion and one line per covariance model and number of group classes.
It is the view `plot(x)` draws for an enumeration, as a function of its
own so that its arguments are listed (and completed by an editor)
without going through `what =`. A candidate that failed to converge is
marked with a cross on its panel's floor, so a gap in a line reads as a
failure and not as a missing candidate.

## Usage

``` r
plot_enumeration(
  x,
  criterion = c("aic", "bic_groups", "bic_individual", "icl_individual"),
  combine = TRUE,
  labels = TRUE,
  mark_minimum = TRUE,
  main = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A `multilpa_enumeration` result from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  or
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md).

- criterion:

  One or more criterion columns of `as.data.frame(x)`, such as
  `"bic_individual"` or `"sabic_groups"`, or `"all"` for every
  information criterion in the grid. The default draws AIC, BIC under
  both sample-size conventions and ICL counted over individuals. A
  criterion identical at both levels, as in a single-level grid, is
  drawn once.

- combine:

  `TRUE` draws several criteria as panels of one plot; `FALSE` returns
  one plot per criterion.

- labels:

  `TRUE` labels each series at its right end.

- mark_minimum:

  `TRUE` rings each criterion's lowest value among converged candidates.
  This marks an extremum; it does not select a model.

- main, subtitle:

  Title and subtitle. `NULL` uses the view's own.

## Value

A ggplot object, one panel per criterion; with `combine = FALSE` and
several criteria, a `latents_plots` list of ggplot objects, one per
criterion, named by its criterion column.

## Errors

`latents_unknown_criterion` for a criterion that is not an information
criterion of the grid; `latents_nothing_to_plot` when no converged
candidate has a finite value.

## See also

[`plot.multilpa_enumeration()`](https://pak.dynasite.org/latents/reference/plot.multilpa_enumeration.md)
for the profile tree,
[`summary.multilpa_enumeration()`](https://pak.dynasite.org/latents/reference/latents-summary.md)
for the candidate each criterion prefers.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(7)
  example_data <- data.frame(score_a = rnorm(120), score_b = rnorm(120))
  candidates <- enumerate_lpa(example_data, c("score_a", "score_b"),
                              n_profiles = 1:3, n_starts = 2)
  plot_enumeration(candidates)
  plot_enumeration(candidates, criterion = "all")
}
```
