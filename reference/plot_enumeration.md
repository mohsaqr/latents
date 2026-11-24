# Plot the information criteria of an enumeration grid against the number of classes

`plot_enumeration()` plots the information criteria of an enumeration
grid against the number of classes: one panel per criterion and one line
per number of group classes (or structure). Lower is preferred; no model
is selected.

## Usage

``` r
plot_enumeration(
  x,
  what = c("aic", "bic_groups", "bic_individual", "icl_individual"),
  combine = TRUE,
  labels = TRUE,
  mark_minimum = TRUE,
  main = NULL,
  subtitle = NULL,
  cutoff = 0.8,
  using = "bic_groups",
  ...
)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  or
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).

- what:

  One or more criteria (columns of `as.data.frame(x)`), or one or more
  views; not both in one call. The criteria come first below; of the
  views, a family grid has only the first three and `"tree"`:

  - criteria of a profile or class grid: `"aic"`, `"kic"`, and `"bic"`,
    `"sabic"`, `"caic"`, `"awe"`, `"icl"`, `"clc"` with the suffix
    `_groups` or `_individual` (sample size: the groups, or the rows).
    Default: `"aic"`, `"bic_groups"`, `"bic_individual"`,
    `"icl_individual"`;

  - criteria of a transition or family grid
    ([`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
    with `time` or `family`): `"aic"`, `"bic"`, `"bic_individual"`, the
    default;

  - criteria of a regression grid: `"aic"`, `"bic"`, `"sabic"`, `"icl"`
    (the default) and `"bic_rows"`;

  - `"all"`: every criterion of the grid, equal ones drawn once;

  - `"entropy"`: each candidate's relative entropy
    ([`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md));

  - `"avepp"`: its lowest average posterior probability
    ([`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md));

  - `"diagnostics"`: `"entropy"` and `"avepp"` together;

  - `"deltas"`: the drop in `using` with each added class
    ([`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md));

  - `"heatmap"`: how far `using` is above its lowest
    ([`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md));

  - `"smallest"`: the smallest class's share, 5% dashed
    ([`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md));

  - `"trust"`: convergence, bounds and replication
    ([`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md));

  - `"selection"`: the last four and `"diagnostics"`
    ([`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md));

  - `"tree"`: how classes split
    ([`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md));
    not with another view.

- combine:

  `TRUE` (default) draws the criteria, or `"entropy"` and `"avepp"`, as
  panels of one plot; `FALSE` returns one plot each, as a list. With
  `"tree"` it goes to
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md),
  whose default is `FALSE`.

- labels:

  `TRUE` (default) labels each line directly; `FALSE` uses a legend
  instead.

- mark_minimum:

  `TRUE` (default) circles each criterion's lowest value among the
  converged candidates; `FALSE` draws no circle.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

- cutoff:

  The AvePP cutoff drawn dashed by `"avepp"`, `"diagnostics"` and
  `"selection"`: one number above 0 and at most 1 (default 0.8).

- using:

  The criterion that `"deltas"`, `"heatmap"` and `"selection"` read; by
  default `"bic_groups"` on a profile or class grid, else `"bic"`.

- ...:

  With `what = "tree"`, the other arguments of
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md).
  Any other argument raises `latents_bad_argument`.

## Value

A ggplot object, or a `latents_plots` list of them when there are
several (`combine = FALSE`, several views, `"selection"`, several
trees). Errors: `latents_unknown_criterion` (`what` is neither a
criterion nor a view), `latents_unsupported_view` (a view this kind of
grid lacks), `latents_nothing_to_plot` (no converged candidate has a
finite value), `latents_bad_argument` (no grid, criteria mixed with
views, a bad value).

## See also

[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)
for any grid or comparison;
[`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md)
lists what a grid offers;
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
for the table.

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_enumeration(grid)

plot_enumeration(grid, what = "bic_individual")

plot_enumeration(grid, what = "entropy")
```
