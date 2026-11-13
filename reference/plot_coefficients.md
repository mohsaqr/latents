# Class coefficients with intervals

Each class's regression coefficients of a mixture regression or growth
mixture with their intervals. The fit's own
`plot(x, what = "coefficients")` takes `level` (default 0.95).

## Usage

``` r
plot_coefficients(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A growth mixture, or a mixture regression from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. An object with nothing to draw raises
`latents_nothing_to_plot`, naming the plot functions that apply; an
object the package did not return raises `latents_not_a_fit`.

## See also

[`plot_trajectories()`](https://pak.dynasite.org/latents/reference/plot_trajectories.md);
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the coefficient table.

## Examples

``` r
fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                          n_classes = 2, id = "student", class_level = "group",
                          random = "wave", n_starts = 2, seed = 1)
plot_coefficients(fit)
```
