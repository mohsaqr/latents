# Plot every pair of continuous indicators against each other, with the cases and an ellipse for each profile

`plot_pairs()` plots every pair of continuous indicators against each
other, with the cases marked by profile and an ellipse for each profile.

## Usage

``` r
plot_pairs(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a profile fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) or
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    with or without covariates, that has two or more continuous
    indicators and was fitted without centering or scaling;

  - a transition fit from
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with
    homogeneous transitions (class `multilpa_transitions`) that has two
    or more continuous indicators;

  - a
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    result, which is drawn from the fit it checked.

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text.

## Value

A ggplot object.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a fit with fewer than two continuous
  indicators, a centred or scaled fit, a general transition fit (class
  `multilpa_lta`), an additive or cross-level fit, a mixture regression,
  a growth mixture or an enumeration grid. For an object of the wrong
  kind, the message names the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

The plot is a matrix of panels with one row and one column per
continuous indicator. A panel off the diagonal is a scatter plot of two
indicators. Each case is a point with the colour and the marker of its
most likely profile, drawn larger the less certain its assignment. Each
profile has an ellipse that holds 95% of its fitted distribution, drawn
around the profile's mean, which is marked. A case missing either
indicator of a panel is not in that panel. A panel on the diagonal shows
the density of that indicator among the cases of each profile that
answered it, scaled by the size of the profile and marked at its peak; a
profile with fewer than three answers there, or all alike, has no
density.

The view needs two or more continuous indicators and draws the stored
data of the fit. A fit with person or grand centering or scaling is
refused, because its stored data are not on the scale of the model. For
a transition fit, a case is one answered occasion of a unit, counted for
its most likely state.

## See also

[`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md)
and
[`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md)
for the other views of the cases;
[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
for the profile means.

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_pairs(fit)
```
