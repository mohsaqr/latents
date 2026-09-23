# A transition network for each latent group class

Builds one
[`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) model
per latent group class of a fitted transition model, collected into the
`group_tna` object tna's grouped verbs expect. A class row with no
expected outgoing moves keeps the fit's unestimated transition
probabilities and raises `latents_empty_transition_row`.

## Usage

``` r
get_group_tna(x, ...)

# S3 method for class 'multilpa_transitions'
get_group_tna(x, label = "Group class", ...)
```

## Arguments

- x:

  A fitted model from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md).

- ...:

  Passed to
  [`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) for
  each class.

- label:

  What the classes are called in tna's output.

## Value

An object of class `group_tna`, one `tna` model per latent class:
`centralities()` returns one tidy table with a `group` column,
`compare()` contrasts two classes, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws them
together.

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
#>           profile_1  profile_2
#> profile_1 0.9421132 0.05788681
#> profile_2 0.2952643 0.70473570
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.7534105 0.2465895 
#> 
#> Group class 2 :
#> State Labels : 
#> 
#>    profile_1, profile_2 
#> 
#> Transition Probability Matrix :
#> 
#>           profile_1 profile_2
#> profile_1 0.6043167 0.3956833
#> profile_2 0.1194480 0.8805520
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#>  0.191336  0.808664 
#> 
```
