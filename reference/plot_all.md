# Plot every view an object has and return the plots as one list

`plot_all()` plots every view that an object of the package has, and
returns the plots as one list. Printing the list shows each plot in
turn.

## Usage

``` r
plot_all(x)
```

## Arguments

- x:

  Any object returned by the package: a fitted model, pooled
  imputations, an enumeration grid, a comparison, a bootstrap result, a
  sensitivity analysis or a diagnostics result.

## Value

A `latents_plots` list of ggplot objects, named by view. Printing it
draws every plot; in an interactive session R waits for a key press
before each new plot.

Errors: `latents_not_a_fit` when `x` is not a fit, grid or other result
of the package, for example a data frame; `latents_nothing_to_plot` when
none of the object's views can be drawn; `latents_missing_package` when
ggplot2 is not installed.

## Details

Which plots the list holds depends on the object. The names of the list
are the names of the views:

- a profile, class or transition fit of
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  or [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) (class
  `multilpa_transitions`): the views its data allow among `profiles`,
  `bars`, `heatmap`, `raincloud`, `parallel`, `pairs`, `responses`,
  `probabilities`, `sequences`, `transitions`, `sizes`, `entropy`,
  `posteriors` and `avepp`. A fit with two or more group classes adds
  `group_sizes`, `group_posteriors`, `group_entropy` and `group_avepp`;

- a fit of another family: the views of its own
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), then the
  four classification views `sizes`, `posteriors`, `entropy` and
  `avepp`. A mixture regression has `fitted` and `coefficients`; one
  that classifies persons also `trajectories`, `individuals` and
  `classification`; a growth mixture has those three, `coefficients` and
  `random`; a general transition fit `transitions`, `profiles` and
  `sequences`; an additive fit `means`, `variances` and `intercepts`; a
  cross-level fit `profiles`, `group_means` and `composition`. When such
  a fit classifies at two levels, each classification view is drawn
  twice, as `sizes_individuals`, `sizes_groups` and so on;

- pooled imputations: the pooled estimates of each level the pooled
  table has, among `profile`, `group` and `measurement`;

- an enumeration grid: `criteria` (every criterion, as
  [`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)),
  `entropy`, `avepp`, `deltas`, `heatmap`, `smallest`, `trust` and
  `tree`. A grid of group-class families has `criteria`, `entropy`,
  `avepp` and `tree`. A grid with several trees has one plot per tree,
  named `tree_` and the tree's name;

- a comparison, a bootstrap test, a bootstrap table of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  or a sensitivity analysis: its one plot, named `plot`;

- a
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  result: the plots of the fit it checked.

Only the views the object can draw are returned. A view that was tried
and refused is left out and named in a message
(`Not drawn for this object: ...`).

## See also

[`overview()`](https://pak.dynasite.org/latents/reference/overview.md)
prints the summary and the checks as well; the functions named after
single views, such as
[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md),
[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md)
and
[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_all(fit)










```
