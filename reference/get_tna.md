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
#> profile_1 0.8151867 0.1848133
#> profile_2 0.1373721 0.8626279
#> 
#> Initial Probabilities : 
#> 
#> profile_1 profile_2 
#>  0.367517  0.632483 
```
