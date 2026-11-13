# Trajectories of individual persons

The observed and predicted trajectory of each of several persons of a
growth mixture or a regression that classifies persons, coloured by
class. Each panel is titled with the person's identifier, named as in
the assignment table, and the person's class and posterior probability.

## Usage

``` r
plot_individuals(x, persons = 12L, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A growth mixture, or a mixture regression that classifies persons
  (`class_level = "group"`).

- persons:

  Which persons to draw. A single number is how many (default 12), taken
  from each class's most to least certain members. Two or more values,
  or text, are the persons' identifiers as they appear in the data; a
  single numeric identifier is given as text, `persons = "1001"`.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object, one panel per person. An object with nothing to draw
raises `latents_nothing_to_plot`, naming the plot functions that apply;
an object the package did not return raises `latents_not_a_fit`. A
`persons` that is not a positive number or a vector of identifiers, or
names none of the fitted persons, raises `latents_bad_argument`.

## See also

[`plot_trajectories()`](https://pak.dynasite.org/latents/reference/plot_trajectories.md),
[`plot_classification()`](https://pak.dynasite.org/latents/reference/plot_classification.md).

## Examples

``` r
fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                          n_classes = 2, id = "student", class_level = "group",
                          random = "wave", n_starts = 2, seed = 1)
plot_individuals(fit)

plot_individuals(fit, persons = c(3, 17, 42))
```
