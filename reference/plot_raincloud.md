# Plot the distribution of each continuous indicator in each profile as a density, its quartiles and the points

`plot_raincloud()` plots the observed values of each continuous
indicator, split by profile: a density, a line over the middle half of
the values (first to third quartile) with a point at the median, and
every case as a point.

## Usage

``` r
plot_raincloud(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a profile fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) or
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    with or without covariates, that has at least one continuous
    indicator and was fitted without centering or scaling;

  - a transition fit from
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with
    homogeneous transitions (class `multilpa_transitions`) that has at
    least one continuous indicator;

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
  function has no plot of it: a fit without continuous indicators (a
  latent class fit), a centred or scaled fit, a general transition fit
  (class `multilpa_lta`), an additive or cross-level fit, a mixture
  regression, a growth mixture or an enumeration grid. For an object of
  the wrong kind, the message names the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

There is one panel per continuous indicator, each with its own x-axis in
the indicator's units. The rows of a panel are the profiles, the largest
at the top, each in its own colour. For each profile the panel shows the
density of the indicator as a filled curve, the interquartile range as a
thick line with a white point at the median, and each case as a small
point below. A case counts for its most likely profile.

The view draws the stored data of the fit. A fit with person or grand
centering or scaling is refused, because its stored data are not on the
scale of the model. For a transition fit, a case is one answered
occasion of a unit, counted for its most likely state.

`plot_profiles(x, what = "raincloud", scale = "standardized")` draws the
same view in standard deviations.

## See also

[`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md)
and
[`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md)
for the other views of the cases;
[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
for the profile means.

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_raincloud(fit)
```
