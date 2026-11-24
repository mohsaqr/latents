# List the views a plot call can draw for a fitted model or a grid, with what each one shows

`plot_views()` lists the values that `what` takes in a plot call, with
one line on what each plots. Without an argument it lists the views of
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for a profile,
class or transition fit, and the two views of
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for an
enumeration grid. Given a grid, it lists every `what` value of
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
for that grid and says whether the grid offers it.

## Usage

``` r
plot_views(x = NULL)
```

## Arguments

- x:

  What to list the views of:

  - `NULL` (default): the views of a fit's
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
    (`"profiles"`, `"bars"`, `"heatmap"`, `"raincloud"`, `"parallel"`,
    `"pairs"`, `"responses"`, `"probabilities"`, `"sequences"`,
    `"transitions"`, `"sizes"`, `"entropy"`, `"posteriors"`, `"avepp"`),
    a grid's `"enumeration"` and `"tree"`, and `"all"`;

  - an enumeration grid from
    [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
    [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
    [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
    or
    [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md):
    that grid's criteria, `"all"`, and the views `"entropy"`, `"avepp"`,
    `"diagnostics"`, `"deltas"`, `"heatmap"`, `"smallest"`, `"trust"`,
    `"selection"` and `"tree"`.

  A fitted model is not taken: the views without an argument are the
  same for every profile, class and homogeneous transition fit, and a
  fit draws the ones it has the indicators, occasions or group classes
  for. A general transition fit (class `multilpa_lta`) has its own
  three: `"transitions"`, `"profiles"` and `"sequences"`.

## Value

A data frame with one row per view, in the order given above.

Without `x`, three character columns:

- `type`: the value to pass as `what`;

- `group`: `"measurement"`, `"structure"`, `"diagnostics"`,
  `"selection"` or `"every"`;

- `description`: what the view draws.

With a grid, five columns:

- `type`: the value to pass as `what`;

- `group`: `"criteria"`, `"quality"`, `"selection"` or `"tree"`;

- `supported`: `TRUE` when this kind of grid offers the view;

- `description`: what the view draws;

- `reason`: why a view is not offered; `NA` when it is.

An `x` that is neither `NULL` nor an enumeration grid raises
`latents_bad_argument`.

## See also

[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md),
[`plot_responses()`](https://pak.dynasite.org/latents/reference/plot_responses.md),
[`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md),
[`plot_prevalence()`](https://pak.dynasite.org/latents/reference/plot_prevalence.md),
[`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
and
[`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/plot_diagnostics.md)
draw a fit's views by name;
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md),
[`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md)
and
[`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md)
a grid's.

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
#>                                                                                                                               description
#> 1                                                                 Profile means, medians or modes across indicators, one line per profile
#> 2                                                                             Profile means as grouped bars from zero, with 95% intervals
#> 3                                                                    Profile means in observed standard deviations from the observed mean
#> 4                                                           Each indicator's distribution by assigned profile: density, box, observations
#> 5                                                                           Every case as a line across indicators, one panel per profile
#> 6                                                                          Scatter-plot matrix with each profile's 95% covariance ellipse
#> 7                                                                                Categorical response probabilities, one line per profile
#> 8                                                                      Profile prevalence within each group class, the two-level quantity
#> 9                                                                                Each group's profile at each occasion, one row per group
#> 10                                                              Estimated transition matrix, one panel per group class (a transition fit)
#> 11                                                                              Effective number of cases in each profile, with its share
#> 12                                                                           Each case's entropy relative to a flat posterior, by profile
#> 13                                                                      Posterior probability of each case's assigned profile, by profile
#> 14                                                                   Average posterior probability: assigned profile by posterior profile
#> 15 Information criteria across a candidate grid (plot any enumeration; or name criteria in `what`; plot_views(grid) lists a grid's views)
#> 16                                           How classes split as more are added, one plot per tree (a profile, class or transition grid)
#> 17                                                                    Every view above that this fit has the ingredients for, in one call

grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
plot_views(grid)
#>                type     group supported
#> 1               aic  criteria      TRUE
#> 2               kic  criteria      TRUE
#> 3        bic_groups  criteria      TRUE
#> 4    bic_individual  criteria      TRUE
#> 5      sabic_groups  criteria      TRUE
#> 6  sabic_individual  criteria      TRUE
#> 7       caic_groups  criteria      TRUE
#> 8   caic_individual  criteria      TRUE
#> 9        awe_groups  criteria      TRUE
#> 10   awe_individual  criteria      TRUE
#> 11       icl_groups  criteria      TRUE
#> 12   icl_individual  criteria      TRUE
#> 13       clc_groups  criteria      TRUE
#> 14   clc_individual  criteria      TRUE
#> 15              all  criteria      TRUE
#> 16          entropy   quality      TRUE
#> 17            avepp   quality      TRUE
#> 18      diagnostics   quality      TRUE
#> 19           deltas selection      TRUE
#> 20          heatmap selection      TRUE
#> 21         smallest selection      TRUE
#> 22            trust selection      TRUE
#> 23        selection selection      TRUE
#> 24             tree      tree      TRUE
#>                                                                                  description
#> 1                                              AIC across the candidates; lower is preferred
#> 2                                              KIC across the candidates; lower is preferred
#> 3                                     BIC (groups) across the candidates; lower is preferred
#> 4                                BIC (individuals) across the candidates; lower is preferred
#> 5                                   SABIC (groups) across the candidates; lower is preferred
#> 6                              SABIC (individuals) across the candidates; lower is preferred
#> 7                                    CAIC (groups) across the candidates; lower is preferred
#> 8                               CAIC (individuals) across the candidates; lower is preferred
#> 9                                     AWE (groups) across the candidates; lower is preferred
#> 10                               AWE (individuals) across the candidates; lower is preferred
#> 11                                    ICL (groups) across the candidates; lower is preferred
#> 12                               ICL (individuals) across the candidates; lower is preferred
#> 13                                    CLC (groups) across the candidates; lower is preferred
#> 14                               CLC (individuals) across the candidates; lower is preferred
#> 15                                        Every criterion of the grid, as panels of one plot
#> 16               Relative entropy of every candidate, at each level with two or more classes
#> 17       The lowest average posterior probability of every candidate, with the cutoff dashed
#> 18                                                         Entropy and lowest AvePP together
#> 19 The decrease in the criterion in `using` with each added class and each added group class
#> 20           The criterion in `using` over classes and group classes, relative to the lowest
#> 21                        The share of the smallest class of every candidate, with 5% dashed
#> 22                                  Convergence, boundary and replication of every candidate
#> 23                          Every selection view the grid offers, and diagnostics, as a list
#> 24                  How classes split as more are added, one plot per tree (see plot_tree())
#>    reason
#> 1    <NA>
#> 2    <NA>
#> 3    <NA>
#> 4    <NA>
#> 5    <NA>
#> 6    <NA>
#> 7    <NA>
#> 8    <NA>
#> 9    <NA>
#> 10   <NA>
#> 11   <NA>
#> 12   <NA>
#> 13   <NA>
#> 14   <NA>
#> 15   <NA>
#> 16   <NA>
#> 17   <NA>
#> 18   <NA>
#> 19   <NA>
#> 20   <NA>
#> 21   <NA>
#> 22   <NA>
#> 23   <NA>
#> 24   <NA>
```
