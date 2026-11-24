# Plot how the profiles differ on the indicators: each profile's mean, or another summary, on every indicator

`plot_profiles()` plots how the profiles (classes) of a fitted model
differ on its indicators: by default one line per profile through its
mean of each indicator, with 95% intervals. `what` chooses another view
of them.

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

  A fitted model, or pooled imputations. What each draws:

  - a profile or class fit
    ([`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)):
    every view it has the indicators for;

  - an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit:
    `"lines"`, and, when its class is `multilpa_transitions`, `"bars"`,
    `"heatmap"` and the views of its continuous indicators
    (`"raincloud"`, `"parallel"`, `"pairs"`), with each answered
    occasion a case; nothing when it has only count indicators or
    `measurement = "occasion"`;

  - a growth mixture or mixture regression: `"lines"` only, each class's
    trajectory over time or fitted line over the first numeric
    predictor;

  - an additive or cross-level fit: `"lines"` only;

  - pooled imputations
    ([`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)):
    `"lines"` only, the pooled means with their intervals beside each
    imputation's estimates.

- statistic:

  What `"lines"` joins. It matters only for a profile, class or
  `multilpa_lta` fit with only categorical or ordinal indicators:

  - `"mean"` (default): the probability-weighted category score;

  - `"median"`: the first category where cumulative probability reaches
    0.5;

  - `"mode"`: the most probable category.

- scale:

  `"raw"` (default) draws each indicator in its own units;
  `"standardized"` in observed standard deviations from the observed
  mean (continuous indicators only). For `"lines"`, `"bars"` and
  `"raincloud"` of profile, class and `multilpa_transitions` fits.

- intervals:

  `TRUE` (default) draws the 95% intervals of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  on `"lines"` and `"bars"` of continuous indicators (the subtitle says
  why a fit has none); `FALSE` draws none.

- labels:

  `TRUE` (default) writes each profile's name at the end of its line;
  `FALSE` uses a legend instead. For `what = "lines"`.

- data:

  The data frame the model was fitted to, used only to compute the
  intervals. `NULL` (default) uses the data stored in the fit.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

- what:

  Which view to draw. The linked pages say what a view needs:

  - `"lines"` (default): each profile's mean of every indicator, joined;

  - `"bars"`: the same means as bars
    ([`plot_bars()`](https://pak.dynasite.org/latents/reference/plot_bars.md));

  - `"heatmap"`: the means in standard deviations
    ([`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md));

  - `"raincloud"`: each indicator's distribution
    ([`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md));

  - `"parallel"`: every case as a line, by profile
    ([`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md));

  - `"pairs"`: the cases in every pair of indicators
    ([`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md));

  - `"all"`: every view the fit supports, as a list.

## Value

A ggplot object; for `what = "all"` a `latents_plots` list of the views
the fit can draw, named by view (a message names the others). Errors:
`latents_nothing_to_plot` (the fit has no such view),
`latents_no_continuous` (too few continuous indicators),
`latents_bad_argument` (an argument no view of the fit takes, or
`"lines"` for category labels that are not numbers),
`latents_not_a_fit`. An argument of another view gives a
`latents_unused_argument` warning.

## See also

[`plot_responses()`](https://pak.dynasite.org/latents/reference/plot_responses.md)
for categorical indicators;
[`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_profiles(fit)

plot_profiles(fit, scale = "standardized")

plot_profiles(fit, what = "bars")
```
