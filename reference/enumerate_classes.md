# Enumerate numbers of individual profiles and group classes

Fits every requested combination and retains errors and warnings
alongside successful fits. Information criteria are descriptive: neither
the smallest BIC nor high entropy guarantees the correct number of
classes. Unconverged, boundary, and unreplicated fits are reported, not
silently selected.

## Usage

``` r
enumerate_classes(
  data,
  vars,
  id,
  n_profiles = 1:4,
  n_group_classes = 1:3,
  model = "basic",
  seed = NULL,
  family = "profiles",
  between_variance = c("varying", "equal"),
  ...
)
```

## Arguments

- data:

  Data frame.

- vars:

  Continuous indicator names.

- id:

  The column that identifies the groups the observations are nested in
  (students in schools, reports in students). `NULL` enumerates
  single-level models, with a message;
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  and
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  fit those by name.

- n_profiles:

  Positive integer profile counts to try.

- n_group_classes:

  Positive integer group-class counts to try. With `id = NULL` it is 1
  and may be left out.

- model:

  The covariance structures to cross with the class counts. `"basic"`,
  the default, fits the four structures that combine variances equal or
  varying across profiles with covariances absent or present: `"EEI"`
  (equal variances, no covariances), `"VVI"` (varying variances, no
  covariances), `"EEE"` (one covariance matrix shared by all profiles)
  and `"VVV"` (a covariance matrix for each profile). They are the
  choices that matter most in practice, since omitting covariances
  between correlated indicators makes the grid favour additional
  profiles. `"all"` fits all 14 structures of Celeux and Govaert (1995).
  Otherwise name the structures by their three-letter codes, as
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)'s
  `model` takes them: the letters give the volume, shape and orientation
  of each profile's covariance matrix, each equal (`E`) or varying (`V`)
  across profiles, with `I` for no covariances. `NULL` crosses no
  structures and fits whatever the other arguments ask for. With one
  continuous indicator the structures reduce to equal or varying
  variance, and codes that coincide are fitted once; with only
  categorical indicators there is no structure to cross. A structure set
  through `variance_model`, `covariance_model`, `volume`, `shape` or
  `orientation` replaces the default.

- seed:

  Optional reproducible seed for each fit.

- family:

  `"profiles"` (the default) enumerates the profile model over
  `n_profiles`, `n_group_classes` and `model`. One or more of the
  group-class families `"additive"`, `"dispersion"` and
  `"additive_dispersion"` enumerates those instead, over `family`,
  `n_group_classes` and `between_variance`; `n_profiles` and `model` are
  then refused, and the result is a `latents_family_enumeration` read
  with
  [`get_results.latents_family_enumeration()`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  and
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md).

- between_variance:

  For group-class families: the between-variance restrictions to cross,
  `"varying"`, `"equal"` or both (the default). The dispersion family is
  always fitted with `"equal"`.

- ...:

  Further arguments to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);
  for group-class families only `n_starts`, `max_iter`, `tol` and
  `min_variance`.

## Value

An object of class `multilpa_enumeration`. Read it with the verbs that
describe it rather than by reaching into it:
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per candidate model with every criterion and diagnostic,
[`summary()`](https://rdrr.io/r/base/summary.html) gives one row per
information criterion naming the candidate that minimises it,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one
criterion across the grid, and
[`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
returns the fitted model for one cell of the grid.

## Conditions

`latents_bad_argument` when `model` names an unknown structure, when it
is given together with another way of setting the structure, or when it
is given although every indicator is categorical.

## References

Celeux, G., & Govaert, G. (1995). Gaussian parsimonious clustering
models. *Pattern Recognition*, 28(5), 781–793.

## See also

[`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
to take one fitted model out of the grid,
[`summary.multilpa_enumeration()`](https://pak.dynasite.org/latents/reference/summary.multilpa_enumeration.md)
for the criterion-by-criterion comparison.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
                              n_group_classes = 1, n_starts = 2, seed = 1)
as.data.frame(candidates)
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic
#> 1          1               1   EEI      -130.6550            2 265.3100
#> 2          2               1   EEI      -130.5791            4 269.1581
#> 3          1               1   VVI      -130.6550            2 265.3100
#> 4          2               1   VVI      -130.5920            5 271.1839
#>        kic bic_groups bic_individual sabic_groups sabic_individual caic_groups
#> 1 270.3100   265.9152       270.5204     259.9237         264.2039    267.9152
#> 2 276.1581   270.3684       279.5788     258.3855         266.9458    274.3684
#> 3 270.3100   265.9152       270.5204     259.9237         264.2039    267.9152
#> 4 279.1839   272.6968       284.2098     257.7182         268.4185    277.6968
#>   caic_individual awe_groups awe_individual icl_groups icl_individual
#> 1        272.5204   276.5204       285.7307   265.9152       270.5204
#> 2        283.5788   291.5788       386.0639   270.3684       355.6432
#> 3        272.5204   276.5204       285.7307   265.9152       270.5204
#> 4        289.2098   299.2098       449.6156   272.6968       411.5897
#>   clc_groups clc_individual profile_entropy group_entropy converged boundary
#> 1   261.3100       261.3100              NA            NA      TRUE    FALSE
#> 2   261.1581       337.2225      0.45131114            NA      TRUE    FALSE
#> 3   261.3100       261.3100              NA            NA      TRUE    FALSE
#> 4   261.1839       388.5639      0.08114782            NA      TRUE    FALSE
#>   n_best_replicated warnings error
#> 1                 2           <NA>
#> 2                 2           <NA>
#> 3                 2           <NA>
#> 4                 2           <NA>
summary(candidates)
#> Class enumeration: 4 candidates, 4 converged, 0 failed to fit
#> 3 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 11 of its 26 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>  n_profiles model log_likelihood n_parameters   aic bic_groups bic_individual
#>           1   EEI         -130.7            2 265.3      265.9          270.5
#>           2   EEI         -130.6            4 269.2      270.4          279.6
#>           1   VVI         -130.7            2 265.3      265.9          270.5
#>           2   VVI         -130.6            5 271.2      272.7          284.2
#>  icl_individual profile_entropy converged boundary
#>           270.5              NA      TRUE    FALSE
#>           355.6         0.45131      TRUE    FALSE
#>           270.5              NA      TRUE    FALSE
#>           411.6         0.08115      TRUE    FALSE
#> 
#> -- criteria --------------------------------------------------------
#>  criterion  convention n_profiles n_group_classes model value
#>        aic        <NA>          1               1   EEI 265.3
#>        kic        <NA>          1               1   EEI 270.3
#>        bic      groups          1               1   EEI 265.9
#>        bic individuals          1               1   EEI 270.5
#>      sabic      groups          2               1   VVI 257.7
#>      sabic individuals          1               1   EEI 264.2
#>       caic      groups          1               1   EEI 267.9
#>       caic individuals          1               1   EEI 272.5
#>        awe      groups          1               1   EEI 276.5
#>        awe individuals          1               1   EEI 285.7
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> 2 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.

# Single level, crossing class counts with two named structures:
single <- enumerate_lpa(iris, c("Sepal.Length", "Sepal.Width",
                                "Petal.Length", "Petal.Width"),
                        n_profiles = 1:3, model = c("VVV", "EEE"),
                        n_starts = 2, seed = 1)
summary(single)
#> Class enumeration: 6 candidates, 6 converged, 0 failed to fit
#> 2 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 10 of its 26 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>  n_profiles model log_likelihood n_parameters   aic   bic icl_individual
#>           1   VVV         -379.9           14 787.8 830.0          830.0
#>           2   VVV         -214.4           29 486.7 574.0          574.0
#>           3   VVV         -180.2           44 448.4 580.8          590.6
#>           1   EEE         -379.9           14 787.8 830.0          830.0
#>           2   EEE         -296.4           19 630.9 688.1          688.1
#>           3   EEE         -256.4           24 560.7 633.0          645.6
#>  profile_entropy converged boundary
#>               NA      TRUE    FALSE
#>           0.9999      TRUE    FALSE
#>           0.9704      TRUE    FALSE
#>               NA      TRUE    FALSE
#>           0.9999      TRUE    FALSE
#>           0.9616      TRUE    FALSE
#> 
#> -- criteria --------------------------------------------------------
#>  criterion  convention n_profiles n_group_classes model value
#>        aic        <NA>          3               1   VVV 448.4
#>        kic        <NA>          3               1   VVV 495.4
#>        bic      groups          2               1   VVV 574.0
#>        bic individuals          2               1   VVV 574.0
#>      sabic      groups          3               1   VVV 441.6
#>      sabic individuals          3               1   VVV 441.6
#>       caic      groups          2               1   VVV 603.0
#>       caic individuals          2               1   VVV 603.0
#>        awe      groups          2               1   VVV 806.3
#>        awe individuals          2               1   VVV 806.3
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> 2 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
