# Plot how many units each class holds, or the size of the smallest class of every candidate in a grid

`plot_sizes()` plots how many units belong to each class of a fitted
model, as bars, or how large the smallest class of every candidate of an
enumeration grid is, as a line.

## Usage

``` r
plot_sizes(x, main = NULL, subtitle = NULL)
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

  - an enumeration grid from
    [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
    [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
    [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
    or
    [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).
    A grid of group-class families (`enumerate_classes(family = )`) has
    no such view and is refused;

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
  function has no plot of it: a grid of group-class families, a grid in
  which no converged candidate has two or more classes, pooled
  imputations, a comparison, a bootstrap result or a sensitivity
  analysis. For an object of the wrong kind, the message names the plot
  functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

**A fitted model.** One horizontal bar per class, the largest at the
top. The length of a bar is the effective number of units in the class:
the sum of the units' posterior probabilities of that class. Each bar is
labelled with that number and its share of the total, and each class has
its own colour.

A fit classifies at one or two levels. The individual level sorts the
cases (rows) into profiles or classes. The group level sorts the groups
(the persons, for example) into group classes. A fit with two or more
group classes besides its profiles is drawn at both levels, as two
plots. An additive fit and a mixture regression with
`class_level = "group"` classify the groups only and give one plot. A
fit with a noise component draws the sizes of its profiles; the noise
cases are left out of the bars and a caption says so.

**An enumeration grid.** The x-axis is the number of classes, the y-axis
the estimated share of the smallest class of each candidate. There is
one line per number of group classes, and per covariance structure when
the grid has several. A dashed line marks 5%. A grid with group classes
has a second panel for the smallest group class.

## See also

[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md)
for the sizes as a table;
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md),
[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md),
[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_sizes(fit)


# Two levels: profiles of the rows and group classes of the students
two_level <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                      "student", n_profiles = 2, n_group_classes = 2,
                      n_starts = 2, seed = 1)
plot_sizes(two_level)



# A grid: the smallest class of every candidate
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_sizes(grid)
```
