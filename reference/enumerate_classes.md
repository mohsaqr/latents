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
  time = NULL,
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

  Optional reproducible seed for each fit: any whole number
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts. Anything
  else raises `latents_bad_argument` before any candidate is fitted.

- family:

  `"profiles"` (the default) enumerates the profile model over
  `n_profiles`, `n_group_classes` and `model`. One or more of the
  group-class families `"additive"`, `"dispersion"` and
  `"additive_dispersion"` enumerates those instead, over `family`,
  `n_group_classes` and `between_variance`; `n_profiles` and `model` are
  then refused, and the result is a `latents_family_enumeration` read
  with
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  and
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

- between_variance:

  For group-class families: the between-variance restrictions to cross,
  `"varying"`, `"equal"` or both (the default). The dispersion family is
  always fitted with `"equal"`.

- time:

  The occasion column: enumerate latent transition models
  ([`lta()`](https://pak.dynasite.org/latents/reference/lta.md)) over
  `n_profiles` (default `2:4`) and `n_group_classes` (default `1`),
  holding every other
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) argument
  given in `...` fixed. `model`, when given, crosses those covariance
  structures. The result is a `latents_transition_enumeration`, read
  with
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  and
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

- ...:

  Further arguments to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (or to [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  with `time`), applied to every candidate alike, such as `n_starts`
  (random starts per candidate; more starts make a local maximum less
  likely, and every candidate gets the same number), `max_iter`,
  `missing`, `categorical`, `centering` and `scaling`. For example,
  `enumerate_classes(data, vars, id, n_profiles = 2:6, n_starts = 20, seed = 1)`.
  Group-class families take only `n_starts`, `max_iter`, `tol` and
  `min_variance`.

## Value

An object of class `multilpa_enumeration`. Read it with the verbs that
describe it rather than by reaching into it:
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per candidate model with every criterion and diagnostic (including
both entropies and, at each level, the smallest, largest and mean
average posterior probability and the number of classes below a cutoff),
`get_results(x, "best")` names the candidate each criterion's minimum
falls on, `get_results(x, "classification")` gives every class's average
posterior probability and relative entropy for every candidate,
`get_results(x, "deltas")` and `"knees"` give the model-selection
evidence, [`summary()`](https://rdrr.io/r/base/summary.html) prints each
criterion's minimum and every table,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
information criteria (and, with `what`, classification quality and
model-selection views),
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
returns the fitted model of one candidate, and
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
puts candidates side by side. Each candidate has a name, the first
column of the table (`"profiles3_groups2"`), which these verbs take.

## Conditions

`latents_bad_argument` when `model` names an unknown structure, when it
is given together with another way of setting the structure, when it is
given although every indicator is categorical, or when `seed` is not a
seed [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts.

## References

Celeux, G., & Govaert, G. (1995). Gaussian parsimonious clustering
models. *Pattern Recognition*, 28(5), 781–793.

## See also

[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
to take one fitted model out of the grid,
[`summary()`](https://rdrr.io/r/base/summary.html) for the
criterion-by-criterion comparison.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
                              n_group_classes = 1, n_starts = 2, seed = 1)
as.data.frame(candidates)
#>               candidate n_profiles n_group_classes model log_likelihood
#> 1 profiles1_groups1_EEI          1               1   EEI      -130.6550
#> 2 profiles2_groups1_EEI          2               1   EEI      -130.5791
#> 3 profiles1_groups1_VVI          1               1   VVI      -130.6550
#> 4 profiles2_groups1_VVI          2               1   VVI      -130.5920
#>   n_parameters      aic      kic bic_groups bic_individual sabic_groups
#> 1            2 265.3100 270.3100   265.9152       270.5204     259.9237
#> 2            4 269.1581 276.1581   270.3684       279.5788     258.3855
#> 3            2 265.3100 270.3100   265.9152       270.5204     259.9237
#> 4            5 271.1839 279.1839   272.6968       284.2098     257.7182
#>   sabic_individual caic_groups caic_individual awe_groups awe_individual
#> 1         264.2039    267.9152        272.5204   276.5204       285.7307
#> 2         266.9458    274.3684        283.5788   291.5788       386.0639
#> 3         264.2039    267.9152        272.5204   276.5204       285.7307
#> 4         268.4185    277.6968        289.2098   299.2098       449.6156
#>   icl_groups icl_individual clc_groups clc_individual profile_entropy
#> 1   265.9152       270.5204   261.3100       261.3100              NA
#> 2   270.3684       355.6432   261.1581       337.2225      0.45131114
#> 3   265.9152       270.5204   261.3100       261.3100              NA
#> 4   272.6968       411.5897   261.1839       388.5639      0.08114782
#>   group_entropy avepp_min avepp_max avepp_mean avepp_below group_avepp_min
#> 1            NA        NA        NA         NA          NA              NA
#> 2            NA 0.6413011 0.8527078  0.7470045           1              NA
#> 3            NA        NA        NA         NA          NA              NA
#> 4            NA 0.6380900 0.6405463  0.6393181           2              NA
#>   group_avepp_max group_avepp_mean group_avepp_below converged boundary
#> 1              NA               NA                NA      TRUE    FALSE
#> 2              NA               NA                NA      TRUE    FALSE
#> 3              NA               NA                NA      TRUE    FALSE
#> 4              NA               NA                NA      TRUE    FALSE
#>   n_best_replicated warnings error
#> 1                 2           <NA>
#> 2                 2           <NA>
#> 3                 2           <NA>
#> 4                 2           <NA>
summary(candidates)
#> Class enumeration: 4 candidates, 4 converged, 0 failed to fit
#> 3 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 15 of its 35 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>              candidate n_profiles model log_likelihood n_parameters   aic
#>  profiles1_groups1_EEI          1   EEI         -130.7            2 265.3
#>  profiles2_groups1_EEI          2   EEI         -130.6            4 269.2
#>  profiles1_groups1_VVI          1   VVI         -130.7            2 265.3
#>  profiles2_groups1_VVI          2   VVI         -130.6            5 271.2
#>  bic_groups bic_individual icl_individual profile_entropy avepp_min avepp_mean
#>       265.9          270.5          270.5              NA        NA         NA
#>       270.4          279.6          355.6         0.45131    0.6413     0.7470
#>       265.9          270.5          270.5              NA        NA         NA
#>       272.7          284.2          411.6         0.08115    0.6381     0.6393
#>  avepp_below converged boundary
#>           NA      TRUE    FALSE
#>            1      TRUE    FALSE
#>           NA      TRUE    FALSE
#>            2      TRUE    FALSE
#> 
#> -- best ------------------------------------------------------------
#>  criterion  convention             candidate n_profiles n_group_classes model
#>        aic        <NA> profiles1_groups1_EEI          1               1   EEI
#>        kic        <NA> profiles1_groups1_EEI          1               1   EEI
#>        bic      groups profiles1_groups1_EEI          1               1   EEI
#>        bic individuals profiles1_groups1_EEI          1               1   EEI
#>      sabic      groups profiles2_groups1_VVI          2               1   VVI
#>      sabic individuals profiles1_groups1_EEI          1               1   EEI
#>       caic      groups profiles1_groups1_EEI          1               1   EEI
#>       caic individuals profiles1_groups1_EEI          1               1   EEI
#>        awe      groups profiles1_groups1_EEI          1               1   EEI
#>        awe individuals profiles1_groups1_EEI          1               1   EEI
#>  value
#>  265.3
#>  270.3
#>  265.9
#>  270.5
#>  257.7
#>  264.2
#>  267.9
#>  272.5
#>  276.5
#>  285.7
#>    ... 4 more rows.  get_results(x, what = "best")
#> 
#> -- classification --------------------------------------------------
#>              candidate n_profiles n_group_classes model       level class
#>  profiles2_groups1_EEI          2               1   EEI individuals     1
#>  profiles2_groups1_EEI          2               1   EEI individuals     2
#>  profiles2_groups1_VVI          2               1   VVI individuals     1
#>  profiles2_groups1_VVI          2               1   VVI individuals     2
#>  n_modal proportion_modal estimated_n estimated_proportion average_posterior
#>       94             0.94       82.31               0.8231            0.8527
#>        6             0.06       17.69               0.1769            0.6413
#>       31             0.31       44.58               0.4458            0.6381
#>       69             0.69       55.42               0.5542            0.6405
#>  odds_correct_classification relative_entropy
#>                        1.244          0.45131
#>                        8.317          0.45131
#>                        2.192          0.08115
#>                        1.434          0.08115
#> 
#> -- deltas ----------------------------------------------------------
#>              candidate n_profiles n_group_classes model  criterion value
#>  profiles1_groups1_EEI          1               1   EEI bic_groups 265.9
#>  profiles2_groups1_EEI          2               1   EEI bic_groups 270.4
#>  profiles1_groups1_VVI          1               1   VVI bic_groups 265.9
#>  profiles2_groups1_VVI          2               1   VVI bic_groups 272.7
#>  delta_profile          from_profile delta_group_class from_group_class
#>             NA                  <NA>                NA             <NA>
#>         -4.453 profiles1_groups1_EEI                NA             <NA>
#>             NA                  <NA>                NA             <NA>
#>         -6.782 profiles1_groups1_VVI                NA             <NA>
#> 
#> -- knees -----------------------------------------------------------
#>   direction                line  criterion                  from
#>    profiles  1 group class, EEI bic_groups profiles1_groups1_EEI
#>    profiles  1 group class, VVI bic_groups profiles1_groups1_VVI
#>  parameters all candidates, EEI bic_groups profiles1_groups1_EEI
#>  parameters all candidates, VVI bic_groups profiles1_groups1_VVI
#>                     to n_points knee_candidate distance knee_without_largest
#>  profiles2_groups1_EEI        2           <NA>       NA                 <NA>
#>  profiles2_groups1_VVI        2           <NA>       NA                 <NA>
#>  profiles2_groups1_EEI        2           <NA>       NA                 <NA>
#>  profiles2_groups1_VVI        2           <NA>       NA                 <NA>
#>  stable                        note
#>      NA fewer than three candidates
#>      NA fewer than three candidates
#>      NA fewer than three candidates
#>      NA fewer than three candidates
#> 
#> 5 tables above, truncated to fit. get_results(x, what = ) returns any
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
#> The candidates table shows 14 of its 35 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>      candidate n_profiles model log_likelihood n_parameters   aic   bic
#>  profiles1_VVV          1   VVV         -379.9           14 787.8 830.0
#>  profiles2_VVV          2   VVV         -214.4           29 486.7 574.0
#>  profiles3_VVV          3   VVV         -180.2           44 448.4 580.8
#>  profiles1_EEE          1   EEE         -379.9           14 787.8 830.0
#>  profiles2_EEE          2   EEE         -296.4           19 630.9 688.1
#>  profiles3_EEE          3   EEE         -256.4           24 560.7 633.0
#>  icl_individual profile_entropy avepp_min avepp_mean avepp_below converged
#>           830.0              NA        NA         NA          NA      TRUE
#>           574.0          0.9999    1.0000     1.0000           0      TRUE
#>           590.6          0.9704    0.9823     0.9900           0      TRUE
#>           830.0              NA        NA         NA          NA      TRUE
#>           688.1          0.9999    1.0000     1.0000           0      TRUE
#>           645.6          0.9616    0.9742     0.9855           0      TRUE
#>  boundary
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#> 
#> -- best ------------------------------------------------------------
#>  criterion  convention     candidate n_profiles n_group_classes model value
#>        aic        <NA> profiles3_VVV          3               1   VVV 448.4
#>        kic        <NA> profiles3_VVV          3               1   VVV 495.4
#>        bic      groups profiles2_VVV          2               1   VVV 574.0
#>        bic individuals profiles2_VVV          2               1   VVV 574.0
#>      sabic      groups profiles3_VVV          3               1   VVV 441.6
#>      sabic individuals profiles3_VVV          3               1   VVV 441.6
#>       caic      groups profiles2_VVV          2               1   VVV 603.0
#>       caic individuals profiles2_VVV          2               1   VVV 603.0
#>        awe      groups profiles2_VVV          2               1   VVV 806.3
#>        awe individuals profiles2_VVV          2               1   VVV 806.3
#>    ... 4 more rows.  get_results(x, what = "best")
#> 
#> -- classification --------------------------------------------------
#>      candidate n_profiles n_group_classes model       level class n_modal
#>  profiles2_VVV          2               1   VVV individuals     1      50
#>  profiles2_VVV          2               1   VVV individuals     2     100
#>  profiles3_VVV          3               1   VVV individuals     1      50
#>  profiles3_VVV          3               1   VVV individuals     2      45
#>  profiles3_VVV          3               1   VVV individuals     3      55
#>  profiles2_EEE          2               1   EEE individuals     1      50
#>  profiles2_EEE          2               1   EEE individuals     2     100
#>  profiles3_EEE          3               1   EEE individuals     1      50
#>  profiles3_EEE          3               1   EEE individuals     2      49
#>  profiles3_EEE          3               1   EEE individuals     3      51
#>  proportion_modal estimated_n estimated_proportion average_posterior
#>            0.3333       50.00               0.3333            1.0000
#>            0.6667      100.00               0.6667            1.0000
#>            0.3333       50.00               0.3333            1.0000
#>            0.3000       44.88               0.2992            0.9823
#>            0.3667       55.12               0.3675            0.9877
#>            0.3333       50.00               0.3333            1.0000
#>            0.6667      100.00               0.6667            1.0000
#>            0.3333       50.00               0.3333            1.0000
#>            0.3267       49.44               0.3296            0.9822
#>            0.3400       50.56               0.3371            0.9742
#>  odds_correct_classification relative_entropy
#>                    1.580e+05           0.9999
#>                           NA           0.9999
#>                    1.801e+16           0.9704
#>                    1.300e+02           0.9704
#>                    1.384e+02           0.9704
#>                    3.297e+05           0.9999
#>                    1.477e+05           0.9999
#>                    6.775e+12           0.9616
#>                    1.122e+02           0.9616
#>                    7.441e+01           0.9616
#> 
#> -- deltas ----------------------------------------------------------
#>      candidate n_profiles n_group_classes model  criterion value delta_profile
#>  profiles1_VVV          1               1   VVV bic_groups 830.0            NA
#>  profiles2_VVV          2               1   VVV bic_groups 574.0       255.960
#>  profiles3_VVV          3               1   VVV bic_groups 580.8        -6.821
#>  profiles1_EEE          1               1   EEE bic_groups 830.0            NA
#>  profiles2_EEE          2               1   EEE bic_groups 688.1       141.881
#>  profiles3_EEE          3               1   EEE bic_groups 633.0        55.134
#>   from_profile delta_group_class from_group_class
#>           <NA>                NA             <NA>
#>  profiles1_VVV                NA             <NA>
#>  profiles2_VVV                NA             <NA>
#>           <NA>                NA             <NA>
#>  profiles1_EEE                NA             <NA>
#>  profiles2_EEE                NA             <NA>
#> 
#> -- knees -----------------------------------------------------------
#>   direction                line  criterion          from            to n_points
#>    profiles  1 group class, EEE bic_groups profiles1_EEE profiles3_EEE        3
#>    profiles  1 group class, VVV bic_groups profiles1_VVV profiles3_VVV        3
#>  parameters all candidates, EEE bic_groups profiles1_EEE profiles3_EEE        3
#>  parameters all candidates, VVV bic_groups profiles1_VVV profiles3_VVV        3
#>  knee_candidate distance knee_without_largest stable
#>   profiles2_EEE   0.1557                 <NA>     NA
#>   profiles2_VVV   0.3678                 <NA>     NA
#>   profiles2_EEE   0.1557                 <NA>     NA
#>   profiles2_VVV   0.3678                 <NA>     NA
#>                                      note
#>  three candidates: stability not assessed
#>  three candidates: stability not assessed
#>  three candidates: stability not assessed
#>  three candidates: stability not assessed
#> 
#> 5 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
