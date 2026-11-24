# Plot each person's posterior probability of every class in a growth model as a stacked bar

`plot_classification()` plots one stacked bar per person of a growth
model, split by the person's posterior probability of each class.

## Usage

``` r
plot_classification(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a growth mixture: a fit from
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
    with `random`;

  - a fit from
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
    with `class_level = "group"`, which classifies persons;

  - a
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    result, which is drawn from the fit it checked.

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text.

## Value

A ggplot object.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a mixture regression that classifies
  observations (including one with `n_group_classes` of two or more), a
  fit of any other family or an enumeration grid. For an object of the
  wrong kind, the message names the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

The x-axis holds the persons and the y-axis the posterior probability,
from 0 to 1. Each bar is one person and is divided into one segment per
class, filled in the class's colour, so that a bar of a single colour is
a person classified with certainty. There is one panel per assigned
(most likely) class, the largest first; within a panel the persons are
ordered from the most to the least certain. On the bar of the least
certain person of each panel the segments are labelled with their class,
where they are wide enough. The subtitle gives, for each class, the
average probability of the assigned class.

[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md)
shows only the probability of the assigned class, and takes fits of
every family.

## See also

[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md)
for the certainty of the assigned class alone;
[`plot_individuals()`](https://pak.dynasite.org/latents/reference/plot_individuals.md);
[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md)
for the table.

## Examples

``` r
fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                          n_classes = 2, id = "student", class_level = "group",
                          random = "wave", n_starts = 2, seed = 1)
plot_classification(fit)
```
