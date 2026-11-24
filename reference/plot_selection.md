# Plot the evidence for choosing the number of classes from an enumeration grid

`plot_selection()` plots five views of an enumeration grid that bear on
the choice of a solution: how much a criterion improves with each added
class, the criterion over the whole grid, the smallest class, whether
each estimate can be trusted, and how sharply each candidate classifies.
None of them selects a model. Each view is also a function of its own,
named below.

## Usage

``` r
plot_selection(
  x,
  what = c("all", "deltas", "heatmap", "smallest", "trust", "diagnostics"),
  using = "bic_groups",
  cutoff = 0.8,
  main = NULL,
  subtitle = NULL
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
  A family grid (`enumerate_classes(family = )`) has no numbers of
  profiles to lay the views out by and draws `"diagnostics"` only.

- what:

  Which view to draw:

  - `"all"` (default): the five views below, as a list;

  - `"deltas"`: the decrease in the criterion `using` with each added
    class and each added group class
    ([`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md));

  - `"heatmap"`: the criterion `using` over the numbers of classes and
    group classes, relative to its lowest value
    ([`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md));

  - `"smallest"`: the share of the smallest class of every candidate,
    with 5% drawn dashed
    ([`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md));

  - `"trust"`: whether each candidate converged, whether it has an
    estimate at a bound, and how many starts reached its best log
    likelihood
    ([`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md));

  - `"diagnostics"`: the relative entropy and the lowest average
    posterior probability (AvePP) of every candidate, with `cutoff`
    drawn dashed
    ([`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md),
    [`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md)).

- using:

  The information criterion that `"deltas"` and `"heatmap"` read: one
  criterion column of the grid. When not given it is `"bic_groups"` on a
  profile or class grid and `"bic"` on a transition or regression grid.

- cutoff:

  The AvePP cutoff drawn dashed by `"diagnostics"`: one number above 0
  and at most 1 (default 0.8).

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

## Value

A ggplot object, which draws the plot when printed; for `"all"` a
`latents_plots` list of the five ggplot objects, named `deltas`,
`heatmap`, `smallest`, `trust` and `diagnostics`, which draws every plot
when printed. Conditions raised:

- `latents_unsupported_view`: a view the kind of grid does not offer
  (every view but `"diagnostics"` on a family grid, `"all"` included);
  the message gives the reason and what the grid offers;

- `latents_bad_argument`: `x` is not an enumeration grid, or `using` is
  not a criterion of the grid.

## See also

[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
and
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
for the tables;
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
for the information criteria;
[`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md)
lists the views a grid offers.

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_selection(grid, what = "deltas")

plot_selection(grid, what = "trust")
```
