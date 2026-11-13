# Every information criterion

For an enumeration grid, every information criterion against the number
of classes, as panels of one plot with each criterion's lowest value
circled. For a comparison from
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md),
the criteria of each compared model.

## Usage

``` r
plot_criteria(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  An enumeration grid, or a comparison from
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. An object with nothing to draw raises
`latents_nothing_to_plot`, naming the plot functions that apply; an
object the package did not return raises `latents_not_a_fit`.

## See also

[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
to choose the criteria;
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_criteria(grid)
```
