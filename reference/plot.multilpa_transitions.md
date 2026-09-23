# Plot a fitted latent transition model

The measurement model is the one
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
fits, so every measurement and classification view it draws is available
here. `what = "transitions"` is the view this family adds: the estimated
transition matrix, one panel per group class.

## Usage

``` r
# S3 method for class 'multilpa_transitions'
plot(
  x,
  what = c("transitions", "profiles", "bars", "heatmap", "responses", "sequences",
    "sizes", "entropy", "posteriors", "avepp", "all"),
  data = NULL,
  scale = c("raw", "standardized"),
  category = "last",
  labels = TRUE,
  main = NULL,
  subtitle = NULL,
  palette = NULL,
  symbols = NULL,
  linetypes = NULL,
  style = .multilpa_style(),
  cell_labels = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted `multilpa_transitions` model.

- what:

  The view to draw. `"transitions"` draws the estimated transition
  matrix, one panel per group class. `"profiles"`, `"bars"` and
  `"heatmap"` draw the measurement model, `"responses"` the categorical
  response curves, `"sequences"` each group's profile at each occasion,
  and `"sizes"`, `"entropy"`, `"posteriors"` and `"avepp"` the
  classification diagnostics. `"all"` draws every view this fit has the
  ingredients for.

- data:

  Optional. The data the model was fitted to. Accepted for consistency
  with
  [`plot.multilpa()`](https://pak.dynasite.org/latents/reference/plot.multilpa.md);
  `"bars"` draws point estimates without intervals here, because this
  family has no standard errors.

- scale:

  `"raw"` keeps the indicators in their own units; `"standardized"`
  divides by each indicator's observed standard deviation.

- category:

  For `"responses"`, which category to draw.

- labels:

  Whether to label series directly.

- main, subtitle:

  Panel title and secondary line.

- palette, symbols, linetypes, style, cell_labels:

  Visual overrides, as in
  [`plot.multilpa()`](https://pak.dynasite.org/latents/reference/plot.multilpa.md).

- ...:

  Further style overrides.

## Value

The fitted model, invisibly, having drawn the requested view.

## See also

[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) and
[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
to draw the transitions as a network instead, and
[`plot.multilpa()`](https://pak.dynasite.org/latents/reference/plot.multilpa.md)
for the same views on a cross-sectional fit.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
)
example_data$score_a <- stats::rnorm(nrow(example_data))
example_data$score_b <- stats::rnorm(nrow(example_data))
fit <- lta(example_data, c("score_a", "score_b"), "person",
           n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
plot(fit, what = "transitions")

plot(fit, what = "profiles")
```
