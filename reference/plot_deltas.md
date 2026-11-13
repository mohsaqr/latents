# Change in the criterion with each added class

The decrease in the information criterion of an enumeration grid with
each added profile and each added group class (BIC over groups on a
profile or class grid, BIC on the others).

## Usage

``` r
plot_deltas(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  An enumeration grid.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. An object with nothing to draw raises
`latents_nothing_to_plot`, naming the plot functions that apply; an
object the package did not return raises `latents_not_a_fit`.

## See also

[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
and
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
for the tables;
[`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_deltas(grid)
```
