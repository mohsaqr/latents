# Build one transition network per latent group class

Builds one
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html) model
per latent group class of a fitted transition model and returns them as
a `group_tna` object, the input of tna's grouped functions. A class row
with no expected outgoing transitions keeps the fit's unestimated
transition probabilities and raises the warning
`latents_empty_transition_row`.

## Usage

``` r
get_group_tna(x, ...)

# S3 method for class 'multilpa_transitions'
get_group_tna(x, label = "Group class", ...)

# S3 method for class 'multilpa_lta'
get_group_tna(x, ...)

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
#> profile_1 0.8805151 0.1194849
#> profile_2 0.3951996 0.6048004
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.8083383 0.1916617 
#> 
#> Group class 2 :
#> State Labels : 
#> 
#>    profile_1, profile_2 
#> 
#> Transition Probability Matrix :
#> 
#>            profile_1 profile_2
#> profile_1 0.70405918 0.2959408
#> profile_2 0.05770587 0.9422941
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.2460338 0.7539662 
#> 
```
