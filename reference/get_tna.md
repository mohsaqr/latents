# Represent transitions between latent states as a network

`get_tna()` constructs a network summarising transitions between latent
profiles or classes across ordered observations. The profiles or classes
form the nodes, and directed edge weights represent transition
probabilities, including the probability of remaining in the same state.
The returned network supports visualisation and analysis with the tna
package.

For a supported latent transition model, the network uses
model-estimated transitions, pooling expected transitions across latent
group classes. For a profile or class model fitted with `time`, it
instead summarises transitions between the most probable class
assignments in each group's ordered sequence.

## Usage

``` r
get_tna(x, ...)
```

## Arguments

- x:

  A time-homogeneous fitted model from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), or a
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  model fitted with `time`, for which the network is based on modal
  profile sequences.

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

For a supported latent transition model, transition probabilities are
the expected transition counts summed across group classes and
normalised by row, so each class contributes in proportion to its share
of the expected transitions, not its class probability. If a profile has
no expected outgoing transitions in any class, its row is the class rows
averaged by class probability, and a warning states that the row is
unestimated; absent transitions are not treated as self-transitions.

For a
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
or
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
fit, [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html)
builds the network from the modal sequences. A fit made without `time`
is refused.

A general [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
fit (occasion-varying or covariate-dependent transitions, second order,
movers and stayers, or ordinal, count or weighted data) is currently
unsupported and raises `latents_unsupported_tna`; its transitions are
returned by `get_results(x, "transitions")`.

## Interpretation

A network built from modal assignments is a descriptive summary and does
not account for classification uncertainty. A pooled network may conceal
differences between latent group classes; use
[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
to examine their separate networks.

## See also

[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
for one network per latent group class.

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
#> profile_1 0.8626339 0.1373661
#> profile_2 0.1848036 0.8151964
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.6324903 0.3675097 
```
