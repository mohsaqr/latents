# Certainty of each unit's classification

The posterior probability of each unit's assigned class, by class, at
each level a fit of any family classifies. A level with one class has
nothing to draw.

## Usage

``` r
plot_posteriors(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A fitted model of any family or a
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  result (drawn from its fit).

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object; a `latents_plots` list with elements `individuals` and
`groups` for a fit that classifies at both levels. An object with
nothing to draw raises `latents_nothing_to_plot`, naming the plot
functions that apply; an object the package did not return raises
`latents_not_a_fit`.

## See also

[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md)
for the table;
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md),
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_posteriors(fit)
```
