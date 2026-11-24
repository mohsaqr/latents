# Plot how well a fitted model separates its classes: sizes, posteriors, entropy and average posterior probabilities

`plot_diagnostics()` plots how well a fitted model separates its
classes. A profile, class or transition fit has four views of how its
observations are classified into profiles: the class sizes, the
posterior probability of each observation's assigned profile, each
observation's entropy and the average posterior matrix. A fit with two
or more group classes has the same four views of how its groups are
classified into group classes. Each view is also a function of its own,
which takes models of the other families too:
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

  One of:

  - a profile or class fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    (with or without membership covariates), or an
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
    class `multilpa_transitions`: every view;

  - a growth mixture or a mixture regression: `"posteriors"` only, each
    unit's posterior probability of its most likely class (`"all"` then
    holds that one view).

  Any other model (an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
  class `multilpa_lta`, an additive or cross-level fit, pooled
  imputations) raises `latents_nothing_to_plot`.
  [`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md),
  [`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md),
  [`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
  and
  [`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md)
  draw those fits, pooled imputations excepted.

- what:

  Which view to draw:

  - `"all"` (default): every view the fit supports, as a list;

  - `"sizes"`: the effective number of observations in each profile (the
    sum of their posterior probabilities), as bars labelled with the
    count and its share;

  - `"posteriors"`: each observation's posterior probability of the
    profile it is assigned to, one strip per profile, with chance (one
    over the number of profiles) marked;

  - `"entropy"`: each observation's entropy relative to a flat posterior
    (0 is certain, 1 is flat), one strip per profile;

  - `"avepp"`: the average posterior probability matrix: for the
    observations assigned to each profile (rows), their mean posterior
    probability of every profile (columns). The diagonal is the
    classification accuracy;

  - `"group_sizes"`, `"group_posteriors"`, `"group_entropy"`,
    `"group_avepp"`: the same four views of the groups in their group
    classes.

  `"all"` leaves out what the fit cannot draw: the group views of a fit
  with fewer than two group classes, `"posteriors"` and `"entropy"` of a
  fit with one profile, and every view but `"sizes"` of a fit with a
  noise component.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

## Value

A ggplot object, which draws the plot when printed; for `"all"` a
`latents_plots` list of ggplot objects named by view, which draws every
plot when printed. Conditions raised:

- `latents_no_group_classes`: a group view of a fit with fewer than two
  group classes;

- `latents_unsupported_noise`: a view other than `"sizes"` of a fit with
  a noise component (such a fit always has one group class);

- `latents_nothing_to_plot`: `"posteriors"` or `"entropy"` of a fit with
  one profile, or a model without the view asked for; the message names
  the views and plot functions that apply;

- `latents_not_a_fit`: `x` is not a fitted model.

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

plot_diagnostics(fit, what = "sizes")
```
