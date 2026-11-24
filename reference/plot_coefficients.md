# Plot the regression coefficients of each class of a mixture regression, with 95% confidence intervals

`plot_coefficients()` plots the regression coefficients of each class of
a model fitted with
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
each with its 95% confidence interval.

## Usage

``` r
plot_coefficients(x, main = NULL, subtitle = NULL)
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

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

**Classes of persons** (a growth mixture, or a fit with
`class_level = "group"`). There is one panel per term of the model, each
with its own x-axis. The rows of a panel are the classes, and a row
`Shared by all classes` holds a term that is common to them. A point is
the estimate and the horizontal bar its interval; each class has its own
colour and marker. The estimate, the interval and the Wald p-value are
printed beside each row, and a dashed line marks zero.

**Classes of observations** (any other mixture regression). There is one
panel. The rows are the terms of the model; within a term the classes
are drawn one under the other, each with its own colour and marker. A
point is the estimate and the horizontal line its interval, and a dotted
line marks zero.

`plot(x, what = "coefficients", level = )` draws the same view with
another confidence level.

## See also

[`plot_trajectories()`](https://pak.dynasite.org/latents/reference/plot_trajectories.md)
for the curves the coefficients give;
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the coefficients as a table.

## Examples

``` r
fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                          n_classes = 2, id = "student", class_level = "group",
                          random = "wave", n_starts = 2, seed = 1)
plot_coefficients(fit)
```
