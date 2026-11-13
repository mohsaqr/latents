# Profile means as bars

The profile means of a profile or class fit as grouped bars from zero,
with 95% intervals when the fit has standard errors.

## Usage

``` r
plot_bars(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A profile or class fit with continuous indicators.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. An object with nothing to draw raises
`latents_nothing_to_plot`, naming the plot functions that apply; an
object the package did not return raises `latents_not_a_fit`.

## See also

[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md),
[`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md),
[`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_bars(fit)
```
