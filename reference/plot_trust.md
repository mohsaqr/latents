# Convergence, boundary and replication of every candidate

Whether each candidate of an enumeration grid converged, reached a
boundary and replicated its best log likelihood across starts.

## Usage

``` r
plot_trust(x, main = NULL, subtitle = NULL)
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

[`get_replication()`](https://pak.dynasite.org/latents/reference/get_replication.md)
for the table;
[`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_trust(grid)
```
