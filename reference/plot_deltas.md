# Plot how much the information criterion falls each time a class or a group class is added

`plot_deltas()` plots how much the information criterion of an
enumeration grid falls each time one class, or one group class, is
added.

## Usage

``` r
plot_deltas(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  or
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).
  A grid of group-class families (`enumerate_classes(family = )`) has no
  such view and is refused.

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text.

## Value

A ggplot object.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a fitted model, a comparison or another
  result, a grid of group-class families, or a grid in which no
  converged candidate has a converged neighbour one step smaller. For an
  object of the wrong kind, the message names the plot functions that
  apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

The plot has up to two panels. In the left panel the x-axis is the
number of classes after the step, and each point is the decrease in the
criterion from the candidate with one class fewer and the same number of
group classes. There is one line per number of group classes, named at
its right end. In the right panel the x-axis is the number of group
classes after the step, and each point is the decrease from the
candidate with one group class fewer and the same number of classes; a
grid with one number of group classes has no such panel. A point above
the zero line means that the larger model has the lower criterion. A
small decrease after a large one is where further classes stop paying.

The criterion is `bic_groups` (BIC with the number of groups as the
sample size) on a profile or class grid and `bic` on a transition or
regression grid. `plot_selection(x, what = "deltas", using = )` draws
the same view with another criterion.

## See also

[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
and
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
for the numbers of a profile or class grid;
[`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md)
for the other selection views.

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_deltas(grid)
```
