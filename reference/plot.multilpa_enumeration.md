# Plot a class-enumeration grid

Draws information criteria against the number of profiles as line plots,
one panel per criterion, with one line per number of group classes and
covariance model. Candidates that failed to converge are marked on the
baseline, so a gap in a line reads as a failure and not as a missing
candidate.

## Usage

``` r
# S3 method for class 'multilpa_enumeration'
plot(
  x,
  criterion = c("aic", "bic_groups", "bic_individual", "icl_individual"),
  combine = TRUE,
  labels = TRUE,
  mark_minimum = TRUE,
  main = NULL,
  subtitle = NULL,
  palette = NULL,
  symbols = NULL,
  linetypes = NULL,
  style = .multilpa_style(),
  ...
)
```

## Arguments

- x:

  An `multilpa_enumeration` result from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).

- criterion:

  One or more criterion columns of `as.data.frame(x)`, such as
  `"bic_individual"` or `"sabic_groups"`. Several names draw one panel
  each; the default draws AIC, BIC under both sample-size conventions,
  and ICL counted over individuals.

- combine:

  `TRUE`, the default, draws several criteria as panels of one figure.
  `FALSE` draws each criterion as its own full-size figure, one after
  another, so a report shows them as separate images.

- labels:

  `TRUE` prints a direct label at the right end of each series.

- mark_minimum:

  `TRUE` rings the lowest value among converged candidates. This marks
  an extremum, it does not select a model.

- main, subtitle:

  Panel title and secondary line.

- palette, symbols, linetypes:

  Series aesthetics, recycled over combinations of group-class count and
  covariance model.

- style:

  A list of visual constants, as built by `.multilpa_style()`.

- ...:

  Further named visual constants, merged into `style`.

## Value

The enumeration result, invisibly. Called for its drawing side effect.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
candidates <- enumerate_classes(
  example_data, c("score_a", "score_b"), "school", n_profiles = 1:3,
  n_group_classes = 1, n_starts = 2, seed = 1
)
plot(candidates)
```
