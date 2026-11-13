# Profile means, or a criterion, as a heat map

For a profile or class fit, the profile means in observed standard
deviations from the observed mean, or the response probabilities of
categorical indicators. For an enumeration grid, the information
criterion over profiles and group classes, relative to its lowest value
(BIC over groups on a profile or class grid, BIC on the others).

## Usage

``` r
plot_heatmap(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A profile or class fit, or an enumeration grid.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. An object with nothing to draw raises
`latents_nothing_to_plot`, naming the plot functions that apply; an
object the package did not return raises `latents_not_a_fit`.

## See also

[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md),
[`plot_bars()`](https://pak.dynasite.org/latents/reference/plot_bars.md);
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_heatmap(fit)
```
