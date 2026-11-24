# Plot every information criterion of a grid, or of a model comparison, against the number of classes

`plot_criteria()` plots every information criterion of an enumeration
grid against the number of classes, or the criteria of the models
compared by
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).

## Usage

``` r
plot_criteria(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - an enumeration grid of any kind, from
    [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
    [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
    [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
    or
    [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md);

  - a comparison from
    [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text.

## Value

A ggplot object.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a fitted model, a diagnostics result,
  pooled imputations, a bootstrap result or a sensitivity analysis. For
  an object of the wrong kind, the message names the plot functions that
  apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

**An enumeration grid.** There is one panel per criterion, each with its
own y-axis; lower values are preferred. The x-axis is the number of
classes (the number of group classes for a grid of group-class families,
or when the grid has one number of profiles and several numbers of group
classes). There is one line per number of group classes, and per
covariance structure when the grid has several. In each panel the lowest
value among the converged candidates is circled, and a candidate that
did not converge is marked with a cross at the bottom of the panel. The
circle marks a minimum; it does not choose a model.

The criteria depend on the kind of grid:

- a profile or class grid: AIC, KIC, BIC, SABIC, CAIC, AWE, ICL and CLC.
  BIC, SABIC, CAIC, AWE and ICL exist in two versions, with the number
  of groups or the number of individuals as the sample size; the two are
  drawn as separate panels where they differ and as one panel where they
  are the same;

- a transition grid or a grid of group-class families: AIC, BIC and BIC
  with the individuals as the sample size;

- a regression grid: AIC, BIC, BIC with the rows as the sample size,
  SABIC and ICL.

On these grids a criterion that takes the same values as one already
drawn is not repeated.

[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
draws the criteria you name, and takes the options of this view
(`labels`, `mark_minimum`, `combine`).

**A comparison.** The x-axis lists the compared models and the y-axis is
the value of the criterion; lower is better. AIC, BIC and SABIC are
drawn as one line each, with their own colour and marker.

## See also

[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
to choose the criteria;
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
and
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
for the values as tables;
[`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_criteria(grid)


# A comparison of two fitted models
one <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 1,
           n_starts = 2, seed = 1)
two <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_criteria(compare_models(one, two))
```
