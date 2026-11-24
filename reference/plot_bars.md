# Plot the mean of each continuous indicator in each profile as bars

`plot_bars()` plots the mean of each continuous indicator in each
profile as bars that start at zero.

## Usage

``` r
plot_bars(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a profile fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) or
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    with or without covariates, that has at least one continuous
    indicator;

  - a transition fit from
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) of
    class `multilpa_transitions` with at least one continuous indicator.
    The general transition fits (class `multilpa_lta`) are refused;

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
  latent class fit), an additive or cross-level fit, a mixture
  regression, a growth mixture, a general transition fit or an
  enumeration grid. For an object of the wrong kind, the message names
  the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

The x-axis lists the continuous indicators and the y-axis is the
estimated mean, in each indicator's own units. Within an indicator there
is one bar per profile, side by side. A profile has its own fill colour
and its own marker at the end of its bar, shown in the legend. When the
fit has standard errors, each bar carries its 95% interval; when it has
none, the subtitle says why. Categorical, ordinal and count indicators
are not drawn; a caption names them.

`plot_profiles(x, what = "bars", scale = "standardized")` draws the same
bars in standard deviations, and `intervals = FALSE` there leaves the
intervals out.

## See also

[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
for the same means as lines and for the options of this view;
[`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md),
[`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_bars(fit)
```
