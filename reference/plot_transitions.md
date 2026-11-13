# Transition probabilities between profiles

The estimated transition matrix of a latent transition model: from the
profile at one occasion (rows) to the profile at the next (columns), one
panel per group class.

## Usage

``` r
plot_transitions(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A fitted model of any family.

- main, subtitle:

  Title and subtitle; `NULL` (default) keeps the view's own.

## Value

A ggplot object. A fit without transitions raises
`latents_nothing_to_plot`.

## See also

[`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md);
[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md)
hands the transitions to the tna package.

## Examples

``` r
fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
           id = "student", n_profiles = 2, time = "sequence",
           n_starts = 2, seed = 1)
plot_transitions(fit)
```
