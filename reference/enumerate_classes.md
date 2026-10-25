# Enumerate numbers of individual profiles and group classes

Fits every requested combination of profile count, group-class count and
covariance structure and returns the candidates with their information
criteria and classification diagnostics. Errors and warnings of each
candidate are kept in the table alongside the successful fits.
Information criteria are descriptive: neither the smallest BIC nor high
entropy guarantees the correct number of classes. Unconverged, boundary
and unreplicated fits are reported, not silently selected.

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
  ...,
  cores = 1L,
  type = NULL
)
```

## Arguments

- data:

  Data frame.

- vars:

  Indicator column names.

- id:

  The column identifying the groups the observations are nested in
  (students in schools, reports in students). Required; `NULL`
  enumerates single-level models, with a message.
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  and
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  fit single-level models directly.

- n_profiles:

  Positive integer profile counts to try. Default `1:4`.

- n_group_classes:

  Positive integer group-class counts to try. Default `1:3`. With
  `id = NULL` it is 1 and may be left out.

- model:

  The covariance structures to cross with the class counts. `"basic"`
  (the default) fits four structures: `"EEI"` (equal variances, no
  covariances), `"VVI"` (varying variances, no covariances), `"EEE"`
  (one covariance matrix shared by all profiles) and `"VVV"` (a
  covariance matrix for each profile). Omitting covariances between
  correlated indicators tends to favour additional profiles. `"all"`
  fits all 14 structures of Celeux and Govaert (1995). Otherwise, name
  structures by their three-letter codes, as
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)'s
  `model` takes them: the letters give the volume, shape and orientation
  of each profile's covariance matrix, each equal (`E`) or varying (`V`)
  across profiles, with `I` for no covariances. `NULL` crosses no
  structures. With one continuous indicator, codes that coincide (equal
  or varying variance) are fitted once; with only categorical indicators
  no structure is crossed. A structure set through `variance_model`,
  `covariance_model`, `volume`, `shape` or `orientation` in `...`
  replaces the default.

- seed:

  Optional seed for each fit: a whole number
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts. Anything
  else raises `latents_bad_argument` before any candidate is fitted.

- family:

  `"profiles"` (the default) enumerates the profile model over
  `n_profiles`, `n_group_classes` and `model`. One or more of the
  group-class families `"additive"`, `"dispersion"` and
  `"additive_dispersion"` enumerates those instead, over `family`,
  `n_group_classes` and `between_variance`; these need `id`,
  `n_profiles` and `model` are then refused, and the result is a
  `latents_family_enumeration`, read with
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  and
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

- between_variance:

  For group-class families: the between-variance restrictions to cross,
  `"varying"`, `"equal"` or both (the default). The dispersion family is
  always fitted with `"equal"`. Refused for the profile family.

- time:

  The occasion column. Enumerates latent transition models
  ([`lta()`](https://pak.dynasite.org/latents/reference/lta.md)) over
  `n_profiles` (default `2:4`) and `n_group_classes` (default `1`), with
  every other
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) argument
  given in `...` held fixed. `model`, when given, crosses those
  covariance structures. Only the profile family is available. The
  result is a `latents_transition_enumeration`, read with
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  and
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

- ...:

  Further arguments to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (or to [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  with `time`), applied to every candidate, such as `n_starts` (random
  starts per candidate), `max_iter`, `missing`, `categorical`,
  `centering`, `scaling`, `drop_empty` and `boundary` (`"exclude"` sets
  aside starts with a variance at `min_variance`, degenerate solutions
  that become more likely as `n_starts` grows). For example,
  `enumerate_classes(data, vars, id, n_profiles = 2:6, n_starts = 20, seed = 1)`.
  Group-class families take only `n_starts`, `max_iter`, `tol` and
  `min_variance`.

- cores:

  Number of processor cores to fit the candidates on (default 1). Each
  candidate sets and restores its own random seed, so with `seed` given
  the result is identical to `cores = 1`. Candidates are fitted by
  forked processes on macOS and Linux, and by a socket cluster on
  Windows and on macOS when R uses Apple's Accelerate BLAS (socket
  workers load the installed package). Warnings and messages from the
  workers are passed on in grid order once all have finished. Without
  `seed`, the random starts differ between runs on any number of cores.
  Lost worker results are retried once with half as many workers
  (rounded up). One `latents_worker_lost` warning reports how many tasks
  were lost and recovered, and suggests fewer cores. A task still lost
  becomes a failed start (for flat start grids) or a failed candidate;
  the other results are retained. On sockets, bounded chunks preserve
  earlier results when a connection breaks; completed tasks within the
  failed chunk may also be retried. Each core holds its own copy of the
  fits it computes;
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  gives the number available. Transition grids (`time`) and profile
  grids (including categorical-only latent class models) distribute
  individual starts across one pool, longest estimated first, with at
  most `cores` workers and no nested workers. Both emEM screening and
  finalist runs use the pool; `search` keeps the same defaults as
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
  Setup, ordered starting values, start selection and finalisation stay
  in the current process. Profile models with membership covariates or
  joined categorical `local_dependence` pairs instead occupy one worker
  per candidate in the same pool, with serial starts. Group-class
  families also keep candidate-level parallelism.

- type:

  Optional model type: `"lca"`, `"lpa"` or `"lta"`. `NULL` (the default)
  keeps the existing routing and measurement settings. `"lca"` treats
  every indicator as categorical, as
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  does: indicators named in `ordinal` or `count` keep those
  distributions, and `categorical` adds to the set rather than replacing
  it. Numeric response codes are categorical with `type = "lca"`.
  `"lpa"` requires continuous indicators and refuses non-empty
  `categorical`, `ordinal` or `count`. Both require `time = NULL`.
  `"lta"` requires `time` and takes its measurement settings from
  `categorical`, as
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) does; it
  does not make numeric indicators categorical. An explicit type
  requires `family = "profiles"`; leave `type = NULL` for group-class
  families.

## Value

An object of class `multilpa_enumeration` (or the transition or family
enumeration described above). Its `type` records the resolved model:
`"lca"`, `"lpa"` or `"lta"`, or `"mixed"`/`"family"` for those grids
under `type = NULL`. The print header shows this type; table columns are
unchanged.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per candidate, starting with its name (`"profiles3_groups2"`), with
every criterion and diagnostic (including both entropies and, at each
level, the smallest, largest and mean average posterior probability and
the number of classes below a cutoff) and the candidate's warnings and
error. `get_results(x, "best")` names the candidate at each criterion's
minimum; `get_results(x, "classification")` gives every class's average
posterior probability and relative entropy for every candidate;
`get_results(x, "deltas")` and `"knees"` give model-selection evidence.
[`summary()`](https://rdrr.io/r/base/summary.html) prints each
criterion's minimum and every table,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
information criteria (and, with `what`, classification quality and
model-selection views),
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
returns the fitted model of one candidate, and
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
puts candidates side by side; these verbs take candidates by name.

## Conditions

`latents_bad_argument` when `model` names an unknown structure, when it
is given together with another way of setting the structure, or when it
is given although no indicator is continuous; when `seed` is not a seed
[`set.seed()`](https://rdrr.io/r/base/Random.html) accepts; when `...`
contains an argument
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
does not take; when `n_group_classes` other than 1 is given with
`id = NULL`; when `id` is missing; when `between_variance` is given for
the profile family; when `n_profiles` or `model` is given for a
group-class family, or `id` is `NULL`; when `time` is combined with a
group-class family or with `between_variance`; and when `cores` is not a
whole number of at least 1; and when `type` is invalid or conflicts with
`time`, `family` or the requested measurement distributions.

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
#> 2         266.9458    274.3684        283.5788   291.5788       385.9162
#> 3         264.2039    267.9152        272.5204   276.5204       285.7307
#> 4         268.4185    277.6968        289.2098   299.2098       449.6156
#>   icl_groups icl_individual clc_groups clc_individual profile_entropy
#> 1   265.9152       270.5204   261.3100       261.3100              NA
#> 2   270.3684       355.4955   261.1581       337.0748      0.45237653
#> 3   265.9152       270.5204   261.3100       261.3100              NA
#> 4   272.6968       411.5897   261.1839       388.5639      0.08114782
#>   group_entropy avepp_min avepp_max avepp_mean avepp_below group_avepp_min
#> 1            NA        NA        NA         NA          NA              NA
#> 2            NA 0.6402627 0.8532616  0.7467621           1              NA
#> 3            NA        NA        NA         NA          NA              NA
#> 4            NA 0.6380900 0.6405463  0.6393181           2              NA
#>   group_avepp_max group_avepp_mean group_avepp_below converged boundary
#> 1              NA               NA                NA      TRUE    FALSE
#> 2              NA               NA                NA      TRUE    FALSE
#> 3              NA               NA                NA      TRUE    FALSE
#> 4              NA               NA                NA      TRUE    FALSE
#>   min_standardized_eigenvalue min_gaussian_mass n_best_replicated warnings
#> 1                   1.0000000         100.00000                 2         
#> 2                   0.8403157          17.63499                 2         
#> 3                   1.0000000         100.00000                 2         
#> 4                   0.8201227          44.58310                 2         
#>   error
#> 1  <NA>
#> 2  <NA>
#> 3  <NA>
#> 4  <NA>
summary(candidates)
#> Class enumeration: 4 candidates, 4 converged, 0 failed to fit
#> 3 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 15 of its 37 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>              candidate n_profiles model log_likelihood n_parameters   aic
#>  profiles1_groups1_EEI          1   EEI         -130.7            2 265.3
#>  profiles2_groups1_EEI          2   EEI         -130.6            4 269.2
#>  profiles1_groups1_VVI          1   VVI         -130.7            2 265.3
#>  profiles2_groups1_VVI          2   VVI         -130.6            5 271.2
#>  bic_groups bic_individual icl_individual profile_entropy avepp_min avepp_mean
#>       265.9          270.5          270.5              NA        NA         NA
#>       270.4          279.6          355.5         0.45238    0.6403     0.7468
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
#>        6             0.06       17.63               0.1763            0.6403
#>       94             0.94       82.37               0.8237            0.8533
#>       31             0.31       44.58               0.4458            0.6381
#>       69             0.69       55.42               0.5542            0.6405
#>  odds_correct_classification relative_entropy
#>                        8.313          0.45238
#>                        1.245          0.45238
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
#> -- tree ------------------------------------------------------------
#>  tree       level        from_candidate from_class from_share
#>   EEI individuals profiles1_groups1_EEI          1          1
#>   EEI individuals profiles1_groups1_EEI          1          1
#>   VVI individuals profiles1_groups1_VVI          1          1
#>   VVI individuals profiles1_groups1_VVI          1          1
#>           to_candidate to_class to_share   mass share_of_from share_of_to
#>  profiles2_groups1_EEI        1   0.1763 0.1763        0.1763           1
#>  profiles2_groups1_EEI        2   0.8237 0.8237        0.8237           1
#>  profiles2_groups1_VVI        1   0.4458 0.4458        0.4458           1
#>  profiles2_groups1_VVI        2   0.5542 0.5542        0.5542           1
#> 
#> 6 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.

# Numeric response codes, explicitly categorical:
codes <- data.frame(g = rep(seq_len(10), each = 10), item = rep(0:1, 50))
classes <- enumerate_classes(codes, "item", "g", type = "lca",
                             n_profiles = 1, n_group_classes = 1,
                             n_starts = 1, seed = 1)

# Single level, crossing class counts with two named structures:
single <- enumerate_lpa(iris, c("Sepal.Length", "Sepal.Width",
                                "Petal.Length", "Petal.Width"),
                        n_profiles = 1:3, model = c("VVV", "EEE"),
                        n_starts = 2, seed = 1)
summary(single)
#> Class enumeration: 6 candidates, 6 converged, 0 failed to fit
#> 2 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 14 of its 37 columns; get_results(x, what = "candidates") returns all of them.
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
#> -- tree ------------------------------------------------------------
#>  tree       level from_candidate from_class from_share  to_candidate to_class
#>   VVV individuals  profiles1_VVV          1     1.0000 profiles2_VVV        1
#>   VVV individuals  profiles1_VVV          1     1.0000 profiles2_VVV        2
#>   VVV individuals  profiles2_VVV          1     0.3333 profiles3_VVV        1
#>   VVV individuals  profiles2_VVV          2     0.6667 profiles3_VVV        1
#>   VVV individuals  profiles2_VVV          1     0.3333 profiles3_VVV        2
#>   VVV individuals  profiles2_VVV          2     0.6667 profiles3_VVV        2
#>   VVV individuals  profiles2_VVV          1     0.3333 profiles3_VVV        3
#>   VVV individuals  profiles2_VVV          2     0.6667 profiles3_VVV        3
#>   EEE individuals  profiles1_EEE          1     1.0000 profiles2_EEE        1
#>   EEE individuals  profiles1_EEE          1     1.0000 profiles2_EEE        2
#>  to_share      mass share_of_from share_of_to
#>    0.3333 3.333e-01     3.333e-01   1.000e+00
#>    0.6667 6.667e-01     6.667e-01   1.000e+00
#>    0.3333 3.333e-01     1.000e+00   1.000e+00
#>    0.3333 4.220e-06     6.330e-06   1.266e-05
#>    0.2992 2.274e-17     6.822e-17   7.601e-17
#>    0.2992 2.992e-01     4.488e-01   1.000e+00
#>    0.3675 2.016e-21     6.047e-21   5.485e-21
#>    0.3675 3.675e-01     5.512e-01   1.000e+00
#>    0.3333 3.333e-01     3.333e-01   1.000e+00
#>    0.6667 6.667e-01     6.667e-01   1.000e+00
#>    ... 6 more rows.  get_results(x, what = "tree")
#> 
#> 6 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
