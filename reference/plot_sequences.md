# Assigned profiles over the occasions

One row per group (person) and one cell per occasion, coloured by the
profile assigned at that occasion, for a fit with `time`.

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

  A fitted model of any family.

- labels:

  `TRUE` (default) groups the rows by group class under labelled strips;
  `FALSE` does not.

- cell_labels:

  `TRUE` (default) prints the profile in each cell when the cells are
  large enough for the text.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. A fit without occasions raises `latents_no_time`.

## See also

[`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
for the transition matrix.

## Examples

``` r
fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
           id = "student", n_profiles = 2, time = "sequence",
           n_starts = 2, seed = 1)
plot_sequences(fit)
```
