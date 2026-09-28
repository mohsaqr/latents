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

# S3 method for class 'multilpa'
get_group_tna(x, label = "Group class", ...)

# S3 method for class 'multilpa_covariates'
get_group_tna(x, label = "Group class", ...)
```

## Arguments

- x:

  A fitted model from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), or a
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  /
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  fit made with `time =`, whose modal profile sequences are then used.

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
#>           profile_1 profile_2
#> profile_1 0.6043454 0.3956546
#> profile_2 0.1194501 0.8805499
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.1913555 0.8086445 
#> 
#> Group class 2 :
#> State Labels : 
#> 
#>    profile_1, profile_2 
#> 
#> Transition Probability Matrix :
#> 
#>           profile_1  profile_2
#> profile_1 0.9421240 0.05787602
#> profile_2 0.2953049 0.70469511
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.7534435 0.2465565 
#> 
```
