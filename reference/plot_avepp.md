# Plot the average posterior probability matrix of a fit, or the lowest average posterior probability of every candidate

`plot_avepp()` plots the average posterior probability (AvePP) matrix of
a fitted model: how much of their probability the units assigned to each
class place on every class. For an enumeration grid it plots the lowest
average posterior probability of every candidate.

## Usage

``` r
plot_avepp(x, main = NULL, subtitle = NULL)
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
  function has no plot of it: a fit with a noise component, a grid in
  which no converged candidate has two or more classes, pooled
  imputations, a comparison, a bootstrap result or a sensitivity
  analysis. For an object of the wrong kind, the message names the plot
  functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

**A fitted model.** A square of tiles. Each row is the class the units
are assigned to (their most likely class), labelled with the number of
units assigned. Each column is a class. A tile holds the mean posterior
probability of the column's class among the units of the row; the number
is printed in the tile and the fill is darker for higher values. The
outlined diagonal is the average posterior probability of each class
among its own units: values near 1 mean that the class is well
separated. A fit with one class draws a single tile.

A fit with two or more group classes besides its profiles is drawn at
both levels, as two plots. A fit with a noise component is refused.

**An enumeration grid.** The x-axis is the number of classes (the number
of group classes for a grid of group-class families), the y-axis the
lowest average posterior probability of each candidate, that is, the
smallest value on its diagonal. There is one line per number of group
classes, and per covariance structure when the grid has several. A
dashed line marks the cutoff 0.8;
`plot_enumeration(x, what = "avepp", cutoff = )` draws the same view
with another cutoff. A grid with group classes has one panel for the
individuals and one for the groups.

## See also

[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md)
for the matrix as a table;
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md),
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_avepp(fit)


# A grid: the lowest average posterior probability of every candidate
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_avepp(grid)
```
