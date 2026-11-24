# Plot whether each candidate of a grid converged, lies on a boundary and replicated its best log likelihood

`plot_trust()` plots whether the estimate of each candidate of an
enumeration grid can be relied on: did it converge, is it on a boundary,
and did more than one start reach its best log likelihood.

## Usage

``` r
plot_trust(x, main = NULL, subtitle = NULL)
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
  result, or a grid of group-class families. For an object of the wrong
  kind, the message names the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

The x-axis is the number of classes and the y-axis the number of group
classes, so each tile is one candidate. A grid with several covariance
structures has one panel per structure. The fill and the text of a tile
give the status of the candidate:

- `replicated`: converged, and two or more starts reached the best log
  likelihood;

- `once`: converged, but only one start reached the best log likelihood,
  so a higher maximum may exist;

- `n/a`: converged, with no replication count stored;

- `boundary`: converged with an estimate on the boundary of the
  parameter space;

- `not conv.`: the best start did not converge;

- `not fitted`: the candidate has no log likelihood.

Under the status the tile gives the number of starts that reached the
best log likelihood, for example `2x`. Replication does not prove that
the global maximum was found.

## See also

[`get_replication()`](https://pak.dynasite.org/latents/reference/get_replication.md)
for the counts of starts as a table;
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
for `converged` and `boundary`;
[`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_trust(grid)
```
