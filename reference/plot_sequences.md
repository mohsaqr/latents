# Plot each unit's sequence of assigned profiles across the occasions as a row of coloured cells

`plot_sequences()` plots each group's (person's) sequence of profiles as
a row of coloured cells: one row per group, one cell per occasion,
coloured by the profile the observation at that occasion is assigned to
(its most likely one). Within a block the rows are sorted by each
group's most frequent profile, so groups with similar sequences sit
together. An occasion at which nothing was observed is shown grey, and
one absent from the data is left blank.

## Usage

``` r
plot_sequences(
  x,
  labels = TRUE,
  cell_labels = TRUE,
  main = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A fit with ordered occasions:

  - a transition fit from
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), of
    either class (`multilpa_transitions` or `multilpa_lta`), except one
    with `measurement = "occasion"`, whose profiles differ by occasion;

  - a profile or class fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    (with or without membership covariates) made with `time`.

- labels:

  `TRUE` (default) puts the rows in blocks under labelled strips: one
  block per group class when the fit has two or more group classes,
  otherwise one per most frequent profile. `FALSE` draws the rows in the
  same order as one block without strips.

- cell_labels:

  `TRUE` (default) prints the profile number in each cell when the grid
  has at most 60 rows and 40 occasions; a larger grid shows the profiles
  by colour, with a dash still marking an occasion that observed nothing
  when the grid also has noise or absent occasions. `FALSE` never prints
  the numbers.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

## Value

A ggplot object, which draws the plot when printed. Conditions raised:

- `latents_no_time`: a profile or class fit made without `time`;

- `latents_nothing_to_plot`: a model with no sequence view (an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit with
  `measurement = "occasion"`, a growth mixture, a mixture regression, an
  additive or cross-level fit, pooled imputations);

- `latents_not_a_fit`: `x` is not a fitted model.

## See also

[`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
for the transition matrix. The sequences as a table:
`get_results(x, "sequences")` for a profile or class model and a
transition fit of class `multilpa_transitions`;
`get_results(x, "assignments")` (one row per observation, with its
profile) for a general transition fit of class `multilpa_lta`.

## Examples

``` r
fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
           id = "student", n_profiles = 2, time = "sequence",
           n_starts = 2, seed = 1)
plot_sequences(fit)
```
