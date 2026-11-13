# Every plot of an object

Every view the object can draw, as a list that draws each plot when
printed: for a fit, its views and the classification views at each
level; for an enumeration grid, the criteria, entropy, AvePP, smallest
class, criterion changes, trust and tree; for a comparison, bootstrap
test or sensitivity analysis, its one plot; for a diagnostics result,
the views of the fit it checked. A view the object lacks the ingredients
for is left out and named in a message.

## Usage

``` r
plot_all(x)
```

## Arguments

- x:

  Any object returned by the package.

## Value

A `latents_plots` list of ggplot objects named by view. An object the
package did not return raises `latents_not_a_fit`.

## See also

[`overview()`](https://pak.dynasite.org/latents/reference/overview.md)
prints the summary and checks as well.

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_all(fit)










```
