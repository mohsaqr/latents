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
`"diagnostics"` or `"selection"`), and `description`.

## Examples

``` r
plot_views()
#>             type       group
#> 1       profiles measurement
#> 2           bars measurement
#> 3        heatmap measurement
#> 4      raincloud measurement
#> 5      responses measurement
#> 6  probabilities   structure
#> 7      sequences   structure
#> 8    transitions   structure
#> 9          sizes   structure
#> 10       entropy diagnostics
#> 11    posteriors diagnostics
#> 12         avepp diagnostics
#> 13   enumeration   selection
#> 14           all       every
#>                                                                      description
#> 1         Profile means across indicators, point size showing profile prevalence
#> 2         Profile means as grouped bars, with 95% intervals when `data` is given
#> 3          Profile means as standard deviations from each indicator's grand mean
#> 4  Each indicator's distribution by assigned profile: density, box, observations
#> 5                       Categorical response probabilities, one line per profile
#> 6             Profile prevalence within each group class, the two-level quantity
#> 7                       Each group's profile at each occasion, one row per group
#> 8      Estimated transition matrix, one panel per group class (a transition fit)
#> 9                      Effective number of cases in each profile, with its share
#> 10                  Per-case entropy contribution within each profile, as ridges
#> 11                      Posterior probability of the assigned profile, as ridges
#> 12          Average posterior probability: assigned profile by posterior profile
#> 13            Information criteria across a candidate grid (plot an enumeration)
#> 14           Every view above that this fit has the ingredients for, in one call
```
