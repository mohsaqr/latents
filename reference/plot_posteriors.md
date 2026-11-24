# Plot the posterior probability of each unit's assigned class, by class

`plot_posteriors()` plots how certain the classification of each unit
is: the posterior probability of the class the unit is assigned to, by
class.

## Usage

``` r
plot_posteriors(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a fitted model of any family: from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).
    This covers fits with covariates, growth mixtures, and the additive
    and cross-level families of
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);

  - a
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    result, which is drawn from the fit it checked.

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text. When two plots are returned, both get the same title and
  subtitle.

## Value

A ggplot object. A fit that classifies at both levels returns a
`latents_plots` list of two ggplot objects named `individuals` and
`groups`; printing the list draws both.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a fit with one class or with a noise
  component, an enumeration grid, pooled imputations, a comparison, a
  bootstrap result or a sensitivity analysis. For an object of the wrong
  kind, the message names the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

One point per classified unit (a case or a group). The x-axis is the
posterior probability of the unit's assigned (most likely) class. It
runs from chance, one divided by the number of classes (dashed line), to
1, which is certainty. The rows are the assigned classes, the largest at
the top, each in its own colour; the points are spread within a row so
that they do not overlap. A short vertical bar marks the mean of each
class, which is the class's average posterior probability.

A fit with two or more group classes besides its profiles is drawn at
both levels, as two plots. A level with one class has nothing to draw,
so a fit with a single profile or class is refused. A fit with a noise
component is refused too. Enumeration grids have no such view; use
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md)
or
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
for a grid.

## See also

[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md)
for the average posterior probability of each class as a table;
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md),
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md);
[`plot_classification()`](https://pak.dynasite.org/latents/reference/plot_classification.md)
for every class probability of each person of a growth model.

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_posteriors(fit)
```
