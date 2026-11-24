# Plot how common each profile is within each group class of a two-level model

`plot_prevalence()` plots how common each profile is in each group
class: the quantity a two-level model adds to a one-level one. There is
one line per group class, labelled with the class's share of the groups,
joining the probability of each profile within that class.

## Usage

``` r
plot_prevalence(x, labels = TRUE, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a profile or class fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    without membership covariates: the profile probabilities within each
    group class. A fit with one group class (or none) is drawn as a
    single line, the overall profile shares;

  - a cross-level fit: the share of each profile among the members of
    each group class, as stacked bars;

  - pooled imputations
    ([`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)):
    the pooled profile-level parameters with their intervals, beside
    each imputation's. They are the profile probabilities within each
    group class, or the membership coefficients of a model with profile
    covariates.

- labels:

  `TRUE` (default) writes each group class's name at the end of its
  line; `FALSE` uses a legend instead. For a profile or class fit; a
  cross-level fit or pooled imputations raise `latents_bad_argument`
  when it is given.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

## Value

A ggplot object, which draws the plot when printed. Conditions raised:

- `latents_nothing_to_plot`: a model with no such view. A fit with
  membership covariates has none because its prevalence varies with each
  unit's covariates; a transition fit from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), a growth
  mixture, a mixture regression and an additive fit have none either.
  The message names the plot functions that apply;

- `latents_not_a_fit`: `x` is not a fitted model.

## See also

[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md)
for the class sizes at each level;
`get_results(x, "profile_probabilities")` for the table.

## Examples

``` r
fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 2, seed = 1)
plot_prevalence(fit)
```
