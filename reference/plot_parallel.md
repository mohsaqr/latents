# Every case as a line across the indicators

Each case's values on the continuous indicators as one line, one panel
per assigned profile.

## Usage

``` r
plot_parallel(x, main = NULL, subtitle = NULL)
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

[`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md),
[`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md),
[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_parallel(fit)
```
