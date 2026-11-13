# Classification summaries of a fit

The four classification views of a profile, class or transition fit:
effective class sizes, the posterior probability of each case's assigned
class, each case's entropy and the average posterior matrix; and the
same four for the groups of a fit with two or more group classes. Each
view is also a function of its own that takes a fit of any family:
[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md),
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md),
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
and
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md).

## Usage

``` r
plot_diagnostics(
  x,
  what = c("all", "sizes", "posteriors", "entropy", "avepp", "group_sizes",
    "group_posteriors", "group_entropy", "group_avepp"),
  main = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A fitted model of any family.

- what:

  `"all"` (default) draws every view the fit supports, as a list;
  `"sizes"`, `"posteriors"`, `"entropy"` or `"avepp"` one view of the
  cases; `"group_sizes"`, `"group_posteriors"`, `"group_entropy"` or
  `"group_avepp"` one view of the groups.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object; a `latents_plots` list for `"all"`. A group view of a
fit with fewer than two group classes raises `latents_no_group_classes`;
a fit with a noise component (always one group class) draws only
`"sizes"` and refuses the other case views with
`latents_unsupported_noise`.

## See also

[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
and
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md)
for the tables behind the views.

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_diagnostics(fit, what = "avepp")
```
