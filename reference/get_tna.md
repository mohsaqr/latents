# Build one transition network for the whole sample

Builds a single
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html) model
from a fitted transition model, aggregated over all latent group
classes, and returns it.

## Usage

``` r
get_tna(x, ...)

# S3 method for class 'multilpa_transitions'
get_tna(x, ...)

# S3 method for class 'multilpa_lta'
get_tna(x, ...)

# S3 method for class 'multilpa'
get_tna(x, ...)

# S3 method for class 'multilpa_covariates'
get_tna(x, ...)
```

## Arguments

- x:

  A fitted model from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md).

- ...:

  Passed to
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html).

## Value

An object of class `tna`, as returned by
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html), to
which the functions of the tna package apply, including
`centralities()`, `communities()`, `cliques()` and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

The network's transition probabilities are the expected transition
counts summed across group classes and normalised by row, so each class
contributes in proportion to its share of the expected transitions, not
its class probability. If a profile has no expected outgoing transitions
in any class, its row is the class rows averaged by class probability,
and a warning states that the row is unestimated; absent transitions are
not treated as self-transitions.

A general [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
fit (occasion-varying or covariate-dependent transitions, second order,
movers and stayers, or ordinal, count or weighted data) has no single
transition matrix and is refused with `latents_unsupported_tna`; its
transitions are returned by `get_results(x, "transitions")`.

A [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
or
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
fit made with `time =` has no estimated transition matrix: the network
is built by
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html) from
each group's sequence of modal profiles in time order.
[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
splits those sequences by each group's modal group class. A fit made
without `time` is refused.

## See also

[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
for one network per latent class.

## Examples

``` r
if (requireNamespace("tna", quietly = TRUE)) {
  activity <- c("browse", "lectures", "forum_read")
  moves <- lta(course_engagement, vars = activity, id = "student",
                           time = "sequence", n_profiles = 2,
                           n_group_classes = 2, n_starts = 2, seed = 1)
  get_tna(moves)
}
#> State Labels : 
#> 
#>    profile_1, profile_2 
#> 
#> Transition Probability Matrix :
#> 
#>           profile_1 profile_2
#> profile_1 0.8626333 0.1373667
#> profile_2 0.1848046 0.8151954
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.6324896 0.3675104 
```
