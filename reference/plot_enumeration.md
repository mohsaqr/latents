# Information criteria across a candidate grid

Each information criterion against the number of classes, one panel per
criterion and one line per number of group classes (or covariance
structure), with each criterion's lowest value circled. `what` also
names the other views of a grid, which
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md),
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md),
[`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md),
[`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md),
[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md),
[`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md)
and
[`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md)
draw by name.

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

  Criterion columns of `as.data.frame(x)`, such as `"bic_groups"`, or
  `"all"`. The default draws AIC, both BICs and ICL on a profile or
  class grid; AIC, BIC and BIC over individuals on a transition or
  family grid; AIC, BIC, SABIC and ICL on a regression grid. Otherwise
  one of `"entropy"`, `"avepp"`, `"diagnostics"` (both), `"deltas"`,
  `"heatmap"`, `"smallest"`, `"trust"`, `"selection"` (those four and
  `"diagnostics"`) or `"tree"`; these cannot be mixed with criteria.

- combine:

  `TRUE` (default) draws the criteria as panels of one plot; `FALSE`
  returns one plot per criterion.

- labels:

  `TRUE` (default) labels each line directly; `FALSE` uses a legend.

- mark_minimum:

  `TRUE` (default) circles each criterion's lowest value among converged
  candidates, without selecting a model.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

- cutoff:

  The AvePP cutoff drawn dashed by `"avepp"`, `"diagnostics"` and
  `"selection"` (default 0.8).

- using:

  The criterion `"deltas"`, `"heatmap"` and `"selection"` use: when not
  given, `"bic_groups"` on a profile or class grid and `"bic"` on the
  others.

- ...:

  With `what = "tree"`, the arguments of
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md).

## Value

A ggplot object; a `latents_plots` list for several criteria with
`combine = FALSE`, or for a `what` of several views. A criterion not in
the grid raises `latents_unknown_criterion`; a view the kind of grid
does not offer raises `latents_unsupported_view`.

## See also

[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)
draws every criterion of any grid or comparison;
[`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md)
lists the views a grid offers;
[`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_enumeration(grid)

plot_enumeration(grid, what = "bic_individual")
```
