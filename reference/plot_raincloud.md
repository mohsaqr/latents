# Distribution of each indicator by profile

Each continuous indicator's distribution by assigned profile: density,
box and the observations.

## Usage

``` r
plot_raincloud(x, main = NULL, subtitle = NULL)
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

[`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md),
[`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md),
[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_raincloud(fit)
```
