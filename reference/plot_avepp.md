# Average posterior probabilities

For a fit of any family, the average posterior probability matrix
(assigned class by posterior class), at each level the fit classifies.
For an enumeration grid, the lowest average posterior of every
candidate, with the 0.8 cutoff dashed.

## Usage

``` r
plot_avepp(x, main = NULL, subtitle = NULL)
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

[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md)
for the table;
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md),
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_avepp(fit)
```
