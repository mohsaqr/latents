# Plot each profile's probability of one response category on every categorical or ordinal indicator

`plot_responses()` plots how likely each profile (class) is to give one
response category. There is one line per profile, joining the
probability of the chosen category on every categorical or ordinal
indicator, on a 0 to 1 axis. The categorical indicators carry 95%
intervals, formed on the logit scale, when the fit has standard errors.

## Usage

``` r
plot_responses(
  x,
  category = "last",
  intervals = TRUE,
  labels = TRUE,
  data = NULL,
  main = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A profile, class or transition fit with categorical or ordinal
  indicators: from
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) or
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (with or without membership covariates), or from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with
  categorical indicators. An
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
  class `multilpa_lta` (ordinal or count indicators, covariates,
  weights, a stayer class, second-order or occasion-specific
  transitions, occasion-specific measurement) is not drawn: read its
  probabilities with `get_results(x, "responses")` or
  `get_results(x, "ordinal")`.

- category:

  Which category's probability is drawn:

  - `"last"` (default): the last category of each indicator, which for a
    no/yes indicator is the probability of "yes";

  - `"first"`: the first category of each indicator;

  - one category label (such as `"yes"`) or one whole-number position
    (such as `2`), which every indicator must have.

- intervals:

  `TRUE` (default) draws the 95% intervals of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md);
  a fit that has no standard errors is drawn without them, and the
  subtitle says why. A probability estimated at 0 or 1 is held there and
  has no interval. `FALSE` draws no intervals and does not compute them.

- labels:

  `TRUE` (default) writes each profile's name at the end of its line;
  `FALSE` uses a legend instead.

- data:

  The data frame the model was fitted to, used only to compute the
  intervals. `NULL` (default) uses the data stored in the fit.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

## Value

A ggplot object, which draws the plot when printed. Conditions raised:

- `latents_no_categorical`: a profile, class or transition fit without
  categorical or ordinal indicators;

- `latents_unknown_category`: `category` is not a category of every
  indicator;

- `latents_nothing_to_plot`: a model of another family, which has no
  response view; the message names the plot functions that apply;

- `latents_not_a_fit`: `x` is not a fitted model.

## See also

[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
for continuous indicators;
[`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md)
for every category at once.

## Examples

``` r
fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
           n_classes = 2, n_starts = 2, seed = 1)
plot_responses(fit)

plot_responses(fit, category = "no")
```
