# Construct transition networks for latent group classes

`get_group_tna()` constructs separate transition networks for latent
group classes, allowing patterns of persistence and movement between
states to be compared. Each network represents the same latent profiles
or classes, making differences in transition probabilities directly
comparable.

For a supported latent transition model, each network uses the estimated
transition matrix and initial state distribution of its group class. For
a profile or class model fitted with `time`, it summarises the modal
state sequences of observed groups assigned to each latent group class.
The returned object supports grouped network analysis and visualisation
with the tna package.

## Usage

``` r
get_group_tna(x, ...)

# S3 method for class 'multilpa_transitions'
get_group_tna(x, label = "Group class", ...)
```

## Arguments

- x:

  A supported time-homogeneous model from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), or a
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  fit made with `time`, whose modal profile sequences are then used.

- ...:

  Passed to
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html) for
  each class.

- label:

  Name given to the classes in tna's output (default `"Group class"`).

## Value

An object of class `group_tna`, a list of `tna` models, one per latent
group class. `centralities()` returns one table with a `group` column,
`compare()` contrasts two classes, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws all
networks.

## Details

For a transition-matrix row with no expected outgoing transitions, the
fitted probabilities are retained and a `latents_empty_transition_row`
warning identifies the row as unestimated. General
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) models are
currently unsupported; see
[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) for
the supported model scope.

## Interpretation

Networks represent latent group classes, rather than individual observed
groups. Networks constructed from modal sequences treat both state and
group-class assignments as known and do not account for uncertainty in
either classification.

## See also

[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) for
one network over the whole sample.

## Examples

``` r
if (requireNamespace("tna", quietly = TRUE)) {
  activity <- c("browse", "lectures", "forum_read")
  moves <- lta(course_engagement, vars = activity, id = "student",
                           time = "sequence", n_profiles = 2,
                           n_group_classes = 2, n_starts = 2, seed = 1)
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
