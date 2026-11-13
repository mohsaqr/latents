# Response probabilities of the categorical indicators

One line per profile (class), joining the probability of one category of
each categorical indicator, with 95% intervals (formed on the logit
scale) when the fit has standard errors.

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

  A fitted model of any family.

- category:

  The category drawn: `"last"` (default), `"first"`, or one label or
  index.

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

## Value

A ggplot object. A fit without categorical indicators raises
`latents_no_categorical`; a `category` no indicator has raises
`latents_unknown_category`.

## See also

[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
for continuous indicators;
[`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md).

## Examples

``` r
fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
           n_classes = 2, n_starts = 2, seed = 1)
plot_responses(fit)
```
