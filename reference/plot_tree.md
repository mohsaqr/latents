# Plot how the classes of each solution split or rearrange as one more class is added

`plot_tree()` plots a clustering tree of an enumeration grid, which
shows whether a solution with one more class splits an existing class or
rearranges them all. A tree has one row per solution and one bar per
class, as wide as the class's share of the units. A band between two
rows shows where the units of a class go in the next solution, as wide
as the share of all units the two classes have in common. Colour follows
a class from row to row, and a new colour marks a class that split off.

## Usage

``` r
plot_tree(
  x,
  main = NULL,
  subtitle = NULL,
  n_group_classes = NULL,
  n_profiles = NULL,
  model = NULL,
  candidates = NULL,
  from = NULL,
  to = NULL,
  level = c("individuals", "groups"),
  combine = FALSE
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
  A family grid (`enumerate_classes(family = )`) has group classes only,
  so its trees are of the groups.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

- n_group_classes, n_profiles, model:

  Which trees and rows to draw; `NULL` (default) keeps every converged
  candidate. At the individual level `n_group_classes` and `model`
  choose the trees and `n_profiles` the rows; at the group level
  `n_profiles` and `model` choose the trees and `n_group_classes` the
  rows. A single count for the rows compares that solution across the
  trees instead, as the rows of one tree:
  `plot_tree(x, n_profiles = 4, n_group_classes = c(2, 5))`.

- candidates:

  Two or more candidate names, drawn as the rows of one tree in the
  order given. Not with `n_group_classes`, `n_profiles` or `model`.

- from, to:

  Two solutions compared as a two-row tree, each named by its class
  counts `c(n_profiles, n_group_classes)` (one count means one group
  class) or by its candidate name:
  `plot_tree(x, from = c(3, 3), to = c(4, 3))`. Given together, not with
  `candidates`.

- level:

  Which units a tree divides:

  - `"individuals"` (default): the observations, into profiles
    (classes). A tree runs down the numbers of profiles at one number of
    group classes and one covariance structure;

  - `"groups"`: the groups, into group classes. A tree runs down the
    numbers of group classes at one number of profiles, so the grid
    needs two or more of them. The default for a family grid.

- combine:

  `FALSE` (default) returns one plot per tree; `TRUE` draws every tree
  as a panel of one plot.

## Value

A ggplot object when there is one tree or `combine = TRUE`; otherwise a
`latents_plots` list of ggplot objects named by tree. Conditions raised:

- `latents_no_tree`: a tree would have fewer than two rows;

- `latents_unknown_candidate`: `candidates`, `from` or `to` names no
  converged candidate of the grid;

- `latents_unsupported_view`: `level = "individuals"` on a family grid;

- `latents_bad_argument`: `x` is not an enumeration grid, or the
  arguments are combined in a way the tree does not take.

## See also

[`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md)
for the flows as a table;
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_tree(grid)

plot_tree(grid, from = 2, to = 3)
```
