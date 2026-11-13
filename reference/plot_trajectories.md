# Class trajectories of the outcome

For a growth mixture or a regression that classifies persons, each
class's mean trajectory over time with its band, the observed
trajectories behind. For another mixture regression, the outcome against
the first numeric predictor, one line per class. The fit's own
`plot(x, what = "trajectories")` also takes `time`, `facet`,
`max_persons` and `spread`, and `plot(x, what = "fitted")` takes
`predictor` and `level`.

## Usage

``` r
plot_trajectories(x, main = NULL, subtitle = NULL)
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

[`plot_individuals()`](https://pak.dynasite.org/latents/reference/plot_individuals.md),
[`plot_coefficients()`](https://pak.dynasite.org/latents/reference/plot_coefficients.md),
[`plot_classification()`](https://pak.dynasite.org/latents/reference/plot_classification.md).

## Examples

``` r
fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                          n_classes = 2, id = "student", class_level = "group",
                          random = "wave", n_starts = 2, seed = 1)
plot_trajectories(fit)
```
