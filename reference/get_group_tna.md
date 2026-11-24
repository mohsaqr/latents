# Get one transition network per latent group class, to compare how their members move between profiles

`get_group_tna()` returns one transition network per latent group class,
as a grouped object of the tna package, so that the group classes can be
compared on how their members stay in and move between profiles. Every
network has the same nodes, the profiles (classes) `profile_1`,
`profile_2`, and so on, and its edge weights are the probabilities of
moving from one profile to another at the next occasion.

Where the probabilities come from depends on the model:

- for a transition fit from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), each
  network is the transition matrix the model estimated for that group
  class, and starts from the class's estimated initial probabilities;

- for a profile or class fit made with `time`, the groups are split by
  their assigned (most likely) group class, and
  [`tna::group_tna()`](https://sonsoles.me/tna/reference/group_model.html)
  counts each network from the sequences of assigned profiles of the
  groups in that class.

## Usage

``` r
get_group_tna(x, ...)

# S3 method for class 'multilpa_transitions'
get_group_tna(x, label = "Group class", ...)
```

## Arguments

- x:

  A transition fit from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), or a
  profile or class fit made with `time`. See *Which fits are supported*.

- ...:

  Further arguments passed to
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html) for
  each class of a transition fit, or to
  [`tna::group_tna()`](https://sonsoles.me/tna/reference/group_model.html)
  for a profile or class fit.

- label:

  One string: what the classes are called in the names of the networks
  and in tna's output (default `"Group class"`, which gives
  `"Group class 1"`, `"Group class 2"`, and so on).

## Value

An object of class `group_tna`: a named list of `tna` networks, one per
latent group class, each as
[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md)
describes. The grouped functions of the tna package apply to it:
`centralities()` returns one table with a `group` column, `compare()`
contrasts two classes, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws every
network.

Conditions raised: `latents_no_time`, `latents_unsupported_tna` and
`latents_not_a_fit` as described under *Which fits are supported*;
`latents_missing_package` when the tna package is not installed;
`latents_bad_argument` when `label` is not one string; and the warning
`latents_empty_transition_row`.

## Details

A transition-matrix row with no expected moves out of its profile keeps
its fitted probabilities, and a warning of class
`latents_empty_transition_row` says that the row is not estimated from
transitions.

## Which fits are supported

The same fits as
[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md),
and mover-stayer fits as well:

- An [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
  class `multilpa_transitions`: always.

- An [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
  class `multilpa_lta`: when its transitions are first order, the same
  at every occasion and without transition covariates; otherwise it
  raises `latents_unsupported_tna`.

- A mover-stayer fit (`mover_stayer = TRUE`): one network per mover
  class and one named `"Stayers"`, whose matrix is the identity.

- A profile or class fit from
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  (with or without membership covariates): when it was made with `time`;
  without it the fit raises `latents_no_time`.

- Any other model raises `latents_not_a_fit`.

## Interpretation

The networks describe latent group classes, not the individual observed
groups. Networks counted from assigned profiles treat both the profile
and the group-class assignments as known and carry the uncertainty of
neither.

## See also

[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) for
one network over the whole sample;
[`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
for the transition matrices as a plot.

## Examples

``` r
if (requireNamespace("tna", quietly = TRUE)) {
  activity <- c("browse", "lectures", "forum_read")
  moves <- lta(course_engagement, vars = activity, id = "student",
               time = "sequence", n_profiles = 2, n_group_classes = 2,
               n_starts = 2, seed = 1)
  get_group_tna(moves)
}
#> Group class 1 :
#> State Labels : 
#> 
#>    profile_1, profile_2 
#> 
#> Transition Probability Matrix :
#> 
#>            profile_1 profile_2
#> profile_1 0.70397865 0.2960213
#> profile_2 0.05768424 0.9423158
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#>  0.245967  0.754033 
#> 
#> Group class 2 :
#> State Labels : 
#> 
#>    profile_1, profile_2 
#> 
#> Transition Probability Matrix :
#> 
#>           profile_1 profile_2
#> profile_1 0.8805107 0.1194893
#> profile_2 0.3951415 0.6048585
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.8082994 0.1917006 
#> 
```
