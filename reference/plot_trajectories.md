# Plot the curve of the outcome in each class of a mixture regression model

`plot_trajectories()` plots the curve of the outcome in each class of a
model fitted with
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md):
over time when the classes are classes of persons, over a predictor
otherwise.

## Usage

``` r
plot_trajectories(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  of any kind, growth mixtures included, or a
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  result, which is drawn from the fit it checked.

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text.

## Value

A ggplot object.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a profile, class, transition, additive or
  cross-level fit, or an enumeration grid. For an object of the wrong
  kind, the message names the plot functions that apply to it;

- `latents_bad_argument` when a mixture regression that classifies
  observations has no numeric predictor in its formula;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

**Classes of persons.** This is a growth mixture (a fit with `random`)
or a fit with `class_level = "group"`. The x-axis is the time variable:
the first numeric variable among the random-effect terms and then the
predictors. The y-axis is the outcome. For each class the plot shows:

- a thick line, the mean trajectory of the class, with the class named
  at its right end;

- a dark band, the 95% confidence band of that mean;

- for a growth mixture, a light band in which 80% of the trajectories of
  the class's persons lie;

- points with a range, the observed mean of the class at each time, with
  its interval (drawn when there are at most 30 distinct times);

- thin lines, the observed trajectories of up to 120 persons, coloured
  by their most likely class.

`plot(x, what = "trajectories")` draws the same view and takes `time`
(another time variable), `facet = TRUE` (one panel per class),
`max_persons`, `level` and, for a growth mixture, `spread`.

**Classes of observations.** For any other mixture regression the x-axis
is the first numeric predictor of the formula and the y-axis the
outcome. Each row of the data is a point with the colour and marker of
its most likely class, and each class has one fitted line, drawn with
the other predictors held at their means (factors at their first level).
`plot(x, what = "fitted", predictor = )` draws the same view over
another numeric predictor.

## See also

[`plot_individuals()`](https://pak.dynasite.org/latents/reference/plot_individuals.md)
for single persons,
[`plot_coefficients()`](https://pak.dynasite.org/latents/reference/plot_coefficients.md)
for the coefficients behind the curves,
[`plot_classification()`](https://pak.dynasite.org/latents/reference/plot_classification.md).

## Examples

``` r
fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                          n_classes = 2, id = "student", class_level = "group",
                          random = "wave", n_starts = 2, seed = 1)
plot_trajectories(fit)


# Classes of observations: the outcome against the first numeric predictor
regression <- mixture_regression(score ~ hours, data = study_hours,
                                 n_classes = 2, n_starts = 2, seed = 1)
plot_trajectories(regression)
```
