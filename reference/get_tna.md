# Get the transitions between a fitted model's profiles as one transition network for the whole sample

`get_tna()` returns the transitions between the profiles (classes) of a
fitted model as one transition network for the whole sample, an object
of the tna package. The profiles are the nodes, named `profile_1`,
`profile_2`, and so on. The weight of the edge from one profile to
another is the probability of moving to it at the next occasion, and the
weight of a loop is the probability of staying.

Where the probabilities come from depends on the model:

- for a transition fit from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), they are
  the transition probabilities the model estimated, pooled over its
  group classes, and the network starts from the estimated initial
  probabilities;

- for a profile or class fit made with `time`, they are counted by
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html)
  from each group's sequence of assigned (most likely) profiles.

## Usage

``` r
get_tna(x, ...)
```

## Arguments

- x:

  A transition fit from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), or a
  profile or class fit made with `time`. See *Which fits are supported*.

- ...:

  Further arguments passed to
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html).

## Value

An object of class `tna`, as
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html)
returns it: a list with `weights` (the square matrix of transition
probabilities, rows the current profile and columns the next, each row
summing to one), `inits` (the probability of starting in each profile),
`labels` (the profile names) and `data` (the sequences, for a profile or
class fit; `NULL` for a transition fit). The functions of the tna
package apply to it, among them `centralities()`, `communities()`,
`cliques()` and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

Conditions raised: `latents_no_time`, `latents_unsupported_tna` and
`latents_not_a_fit` as described under *Which fits are supported*;
`latents_missing_package` when the tna package is not installed; and the
warning `latents_empty_transition_row`.

## Details

For a transition fit, the pooled matrix is built from the expected
transition counts, summed over the group classes and divided by their
row totals. Each group class therefore counts in proportion to the
transitions it is expected to make, not to its class probability. If a
profile has no expected moves out of it in any class, its row is the
class rows averaged by class probability, and a warning of class
`latents_empty_transition_row` says that the row is not estimated from
transitions; a missing transition is never counted as staying.

## Which fits are supported

- An [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
  class `multilpa_transitions` (the basic transition model: continuous
  or categorical indicators and one first-order transition matrix for
  all occasions, without covariates, weights or a stayer class): always.

- An [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
  class `multilpa_lta` (ordinal or count indicators, weights, covariates
  of the initial profiles, occasion-specific measurement): when its
  transitions are first order, the same at every occasion and without
  transition covariates. A fit whose transitions vary by occasion, by
  covariate or by the previous profile (`order = 2`) has no single
  matrix and raises `latents_unsupported_tna`; read its transitions with
  `get_results(x, "transitions")`.

- A mover-stayer fit (`mover_stayer = TRUE`): not by `get_tna()`, which
  raises `latents_unsupported_tna` because the stayer class records no
  expected moves to pool.
  [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  returns one network per mover class and one for the stayers.

- A profile or class fit from
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  (with or without membership covariates): when it was made with `time`;
  without it the fit raises `latents_no_time`.

- Any other model (a mixture regression, a growth mixture, an additive
  or cross-level fit, pooled imputations) raises `latents_not_a_fit`.

## Interpretation

A network counted from assigned profiles treats each assignment as known
and does not carry the uncertainty of the classification. The network of
a transition fit uses probabilities estimated from the posterior
probabilities, so an observation that may belong to two profiles
contributes to both. A network pooled over the group classes can hide
differences between them;
[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
returns their separate networks.

## See also

[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
for one network per latent group class;
[`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
for the transition matrix as a plot; `get_results(x, "transitions")` for
it as a table.

## Examples

``` r
if (requireNamespace("tna", quietly = TRUE)) {
  activity <- c("browse", "lectures", "forum_read")
  moves <- lta(course_engagement, vars = activity, id = "student",
               time = "sequence", n_profiles = 2, n_group_classes = 2,
               n_starts = 2, seed = 1)
  get_tna(moves)
}
#> State Labels : 
#> 
#>    profile_1, profile_2 
#> 
#> Transition Probability Matrix :
#> 
#>           profile_1 profile_2
#> profile_1 0.8626339 0.1373661
#> profile_2 0.1848036 0.8151964
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.6324903 0.3675097 
```
