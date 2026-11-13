# Entropy of the classification

For a fit of any family, each unit's entropy relative to a flat
posterior, by assigned class, at each level the fit classifies. For an
enumeration grid, the relative entropy of every candidate. A level with
one class has no entropy to draw.

## Usage

``` r
plot_entropy(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A fitted model of any family, an enumeration grid or a
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

[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md)
for the table;
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md),
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_entropy(fit)
```
