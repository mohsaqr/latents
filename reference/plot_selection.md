# Model-selection evidence across a candidate grid

Five views that bear on the choice of a solution, none of which selects
one: the decrease in a criterion with each added class (`"deltas"`), the
criterion over profiles and group classes relative to its lowest value
(`"heatmap"`), the smallest class share with 5% dashed (`"smallest"`),
convergence, boundary and replication of every candidate (`"trust"`),
and relative entropy with the lowest AvePP (`"diagnostics"`). Each is
also a function of its own:
[`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md),
[`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md),
[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md),
[`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md),
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
and
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md).

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

- what:

  `"all"` (default) draws the five views as a list, or one of
  `"deltas"`, `"heatmap"`, `"smallest"`, `"trust"` and `"diagnostics"`.

- using:

  The criterion `"deltas"`, `"heatmap"` and `"selection"` use: when not
  given, `"bic_groups"` on a profile or class grid and `"bic"` on the
  others.

- cutoff:

  The AvePP cutoff drawn dashed by `"avepp"`, `"diagnostics"` and
  `"selection"` (default 0.8).

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object; a `latents_plots` list for `"all"`. A view the kind of
grid does not offer raises `latents_unsupported_view`.

## See also

[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
and
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
for the tables;
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
for the information criteria.

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_selection(grid, what = "deltas")
```
