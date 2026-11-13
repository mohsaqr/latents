# Profile means across the indicators

One line per profile (class), joining its mean, median or mode of each
indicator, with 95% intervals when the fit has standard errors. On a
growth mixture or a regression that classifies persons it draws the
class trajectories, on another mixture regression the fitted lines, on
an additive fit the group-class means, on a cross-level fit the profiles
and on pooled imputations the pooled measurement parameters.

## Usage

``` r
plot_profiles(
  x,
  statistic = c("mean", "median", "mode"),
  scale = c("raw", "standardized"),
  intervals = TRUE,
  labels = TRUE,
  data = NULL,
  main = NULL,
  subtitle = NULL,
  what = c("lines", "bars", "heatmap", "raincloud", "parallel", "pairs", "all")
)
```

## Arguments

- x:

  A fitted model of any family.

- statistic:

  `"mean"` (default), `"median"` or `"mode"`; they coincide for
  continuous indicators. For categorical ones the fitted probabilities
  are summarised on the numeric category scale, without intervals.

- scale:

  `"raw"` (default) keeps each indicator in its own units;
  `"standardized"` divides its deviation from the observed mean by the
  observed standard deviation.

- intervals:

  `TRUE` (default) draws the 95% intervals of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  when the fit has standard errors; otherwise the subtitle says why
  there are none.

- labels:

  `TRUE` (default) labels each profile directly; `FALSE` uses a legend.

- data:

  The data the model was fitted to, used only for the intervals. `NULL`
  (default) uses the data stored in the fit.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

- what:

  `"lines"` (default). The other values draw the view of the function
  with that name
  ([`plot_bars()`](https://pak.dynasite.org/latents/reference/plot_bars.md),
  [`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md),
  [`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md),
  [`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md),
  [`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md));
  `"all"` draws every view the fit supports, as a list.

## Value

A ggplot object; a `latents_plots` list for `what = "all"`. A fit the
function has no view for raises `latents_nothing_to_plot`, naming the
plot functions that apply; an argument no view of the fit takes (for
example `scale` on a growth mixture) raises `latents_bad_argument`, and
one that belongs to another view of the same fit is dropped with a
warning.

## See also

[`plot_responses()`](https://pak.dynasite.org/latents/reference/plot_responses.md)
for categorical indicators;
[`plot_bars()`](https://pak.dynasite.org/latents/reference/plot_bars.md),
[`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md),
[`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md),
[`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md),
[`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md);
[`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_profiles(fit)

plot_profiles(fit, scale = "standardized")
```
