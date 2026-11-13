# Class probabilities of every person

One stacked bar per person of a growth mixture or a regression that
classifies persons: the probability of each class, persons ordered by
assigned class and certainty.

## Usage

``` r
plot_classification(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A growth mixture, or a mixture regression that classifies persons
  (`class_level = "group"`).

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. An object with nothing to draw raises
`latents_nothing_to_plot`, naming the plot functions that apply; an
object the package did not return raises `latents_not_a_fit`.

## See also

[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md)
for the certainty of the assigned class alone;
[`plot_individuals()`](https://pak.dynasite.org/latents/reference/plot_individuals.md).

## Examples

``` r
fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                          n_classes = 2, id = "student", class_level = "group",
                          random = "wave", n_starts = 2, seed = 1)
plot_classification(fit)
```
