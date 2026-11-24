# Plot the estimated probability of moving from each profile to each profile between occasions

`plot_transitions()` plots the estimated transition matrix of a latent
transition model as a grid of tiles: the profile at one occasion in the
rows, the profile at the next occasion in the columns, and in each tile
the probability of that move. Each row sums to one. There is one panel
per group class.

## Usage

``` r
plot_transitions(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  A transition fit from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), of
  either class:

  - `multilpa_transitions` (the basic transition model: continuous or
    categorical indicators and one first-order transition matrix for all
    occasions): the estimated matrix of each group class, the panel
    labelled with the class's share of the groups. A row with no
    expected moves out of its profile is not estimated from the data,
    and its label is put in parentheses;

  - `multilpa_lta` (any extension of the basic model: ordinal or count
    indicators, covariates, weights, a stayer class, second-order or
    occasion-specific transitions, occasion-specific measurement): the
    probabilities the model implies, averaged over the observed groups,
    with one panel per group class and, when the transitions differ by
    occasion, per occasion.

- main, subtitle:

  Title and subtitle of the plot, each one string; `NULL` (default)
  keeps the view's own.

## Value

A ggplot object, which draws the plot when printed. Conditions raised:

- `latents_nothing_to_plot`: a model without transitions, which is every
  model not fitted by
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md),
  including a profile or class fit made with `time`; the message names
  the plot functions that apply;

- `latents_not_a_fit`: `x` is not a fitted model.

## See also

[`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md)
for the assigned profiles over time; `get_results(x, "transitions")` for
the probabilities as a table;
[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md)
hands the transitions to the tna package.

## Examples

``` r
fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
           id = "student", n_profiles = 2, time = "sequence",
           n_starts = 2, seed = 1)
plot_transitions(fit)
```
