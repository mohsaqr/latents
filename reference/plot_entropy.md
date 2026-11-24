# Plot how uncertain the classification is, by unit for a fitted model or by candidate for a grid

`plot_entropy()` plots how uncertain the classification is: for a fitted
model, the entropy of each unit, by the class the unit is assigned to;
for an enumeration grid, the relative entropy of every candidate.

## Usage

``` r
plot_entropy(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a fitted model of any family: from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).
    This covers fits with covariates, growth mixtures, and the additive
    and cross-level families of
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);

  - an enumeration grid of any kind, from
    [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
    [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
    [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
    or
    [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md);

  - a
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    result, which is drawn from the fit it checked.

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text. When two plots are returned, both get the same title and
  subtitle.

## Value

A ggplot object. A fit that classifies at both levels returns a
`latents_plots` list of two ggplot objects named `individuals` and
`groups`; printing the list draws both.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a fit with one class or with a noise
  component, a grid in which no converged candidate has two or more
  classes, pooled imputations, a comparison, a bootstrap result or a
  sensitivity analysis. For an object of the wrong kind, the message
  names the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

**A fitted model.** One point per classified unit (a case or a group).
The x-axis is the unit's entropy relative to a flat posterior: 0 means
the unit belongs to one class with certainty, 1 means it is equally
likely to belong to every class (dashed line). The rows are the assigned
(most likely) classes, the largest at the top, each in its own colour;
the points are spread within a row so that they do not overlap. A short
vertical bar marks the mean of each class. The subtitle gives the
relative entropy of the whole classification, which is one minus the
mean of the points: there, higher means sharper.

A fit with two or more group classes besides its profiles is drawn at
both levels, as two plots. A level with one class has no entropy to
draw, so a fit with a single profile or class is refused. A fit with a
noise component is refused too.

**An enumeration grid.** The x-axis is the number of classes (the number
of group classes for a grid of group-class families), the y-axis the
relative entropy of each candidate: higher is sharper. There is one line
per number of group classes, and per covariance structure when the grid
has several. A grid with group classes has one panel for the individuals
and one for the groups. A candidate with one class has no point.

## See also

[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md)
for the relative entropy as a table;
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md),
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_entropy(fit)


# A grid: the relative entropy of every candidate
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_entropy(grid)
```
