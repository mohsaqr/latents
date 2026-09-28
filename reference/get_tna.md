# A transition network for the whole sample

Builds one
[`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) model
from a fitted transition model, aggregated over every latent group
class.

## Usage

``` r
get_tna(x, ...)

# S3 method for class 'multilpa_transitions'
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
  [`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md).

## Value

An object of class `tna`, as
[`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) returns:
every verb of that package applies to it, including `centralities()`,
`communities()`, `cliques()` and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

The aggregate is formed from the expected transition counts, summed
across classes and then normalised, rather than by averaging the class
transition matrices. The two differ: averaging weights each class by how
probable it is, and summing counts weights it by how much transition
mass it actually contributes, which is what a marginal transition
probability means. A class holding a tenth of the students but a fifth
of the observed moves counts for the latter. If a profile has no
expected outgoing moves in any class, its transition row has no
count-based estimate; the network uses the fitted class rows averaged by
class probability and warns that the row is unestimated. It never
interprets absent moves as a certain self-transition.

For a
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
or
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
fit made with `time =`, there is no estimated transition matrix: the
network is built by
[`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) from
each group's sequence of modal profiles in time order, the observed
transitions between the profiles the observations were assigned to.
[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
splits those sequences by each group's modal group class. A fit made
without `time` has no order and is refused.

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
#> profile_1 0.8151872 0.1848128
#> profile_2 0.1373717 0.8626283
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#> 0.3675166 0.6324834 
```
