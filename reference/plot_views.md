# The plots this package can draw

Returns the catalogue of `what =` values accepted by the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods, with
the group each belongs to, so the available views can be listed rather
than recalled from a help page.

## Usage

``` r
plot_views()
```

## Value

A base `data.frame` with one row per plot type and the columns `type`
(the value to pass as `what`), `group` (`"measurement"`, `"structure"`,
`"diagnostics"`, `"selection"` or `"every"`), and `description`.

## Examples

``` r
plot_views()
#>             type       group
#> 1       profiles measurement
#> 2           bars measurement
#> 3        heatmap measurement
#> 4      raincloud measurement
#> 5       parallel measurement
#> 6          pairs measurement
#> 7      responses measurement
#> 8  probabilities   structure
#> 9      sequences   structure
#> 10   transitions   structure
#> 11         sizes diagnostics
#> 12       entropy diagnostics
#> 13    posteriors diagnostics
#> 14         avepp diagnostics
#> 15   enumeration   selection
#> 16          tree   selection
#> 17           all       every
#>                                                                      description
#> 1                 Profile means across indicators, one labelled line per profile
#> 2                    Profile means as grouped bars from zero, with 95% intervals
#> 3           Profile means in observed standard deviations from the observed mean
#> 4  Each indicator's distribution by assigned profile: density, box, observations
#> 5                  Every case as a line across indicators, one panel per profile
#> 6                 Scatter-plot matrix with each profile's 95% covariance ellipse
#> 7                       Categorical response probabilities, one line per profile
#> 8             Profile prevalence within each group class, the two-level quantity
#> 9                       Each group's profile at each occasion, one row per group
#> 10     Estimated transition matrix, one panel per group class (a transition fit)
#> 11                     Effective number of cases in each profile, with its share
#> 12                  Each case's entropy relative to a flat posterior, by profile
#> 13             Posterior probability of each case's assigned profile, by profile
#> 14          Average posterior probability: assigned profile by posterior profile
#> 15            Information criteria across a candidate grid (plot an enumeration)
#> 16                    How profiles split as more are added (plot an enumeration)
#> 17           Every view above that this fit has the ingredients for, in one call
```
