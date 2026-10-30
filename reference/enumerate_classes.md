# Class enumeration and comparison of mixture models

Class enumeration evaluates alternative numbers of latent classes to
inform model selection. `enumerate_classes()` fits a specified set of
candidate models and summarizes their fit, classification quality, and
estimation diagnostics. It supports latent profile, latent class, and
latent transition models, including multilevel specifications.

For multilevel profile and class models, candidates can differ in the
number of observation-level profiles or classes, the number of
group-level classes, and, for continuous indicators, the within-profile
covariance structure. Comparing these specifications allows the analyst
to assess both the latent class structure and the assumptions used to
represent variation within classes.

The function retains the fitted candidates and provides tables for
comparing information criteria, posterior classification probabilities,
class sizes, and convergence diagnostics. Individual models can be
extracted with
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
for further examination and analysis.

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
  type = NULL,
  verbose = cores > 1L
)
```

## Arguments

- data:

  A data frame with one row per observation and columns containing the
  indicators, group identifier and, when applicable, occasion order and
  covariates. Use the same observations and measurement scale throughout
  a candidate comparison.

- vars:

  Character vector naming the indicator columns in `data`. Indicator
  distributions are determined by `type` and the measurement arguments
  supplied through `...`.

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

  Optional integer random seed accepted by
  [`set.seed()`](https://rdrr.io/r/base/Random.html). Used to reproduce
  the starting values for the same data and fitting controls. Invalid
  values raise `latents_bad_argument` before fitting begins.

- family:

  `"profiles"` (the default) enumerates the profile model over
  `n_profiles`, `n_group_classes` and `model`. One or more of the
  group-class families `"additive"`, `"dispersion"` and
  `"additive_dispersion"` enumerates those instead, over `family`,
  `n_group_classes` and `between_variance`; these need `id`,
  `n_profiles` and `model` are then rejected, and the result is a
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

  Additional fitting arguments passed to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  or to [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  when `time` is supplied. The same settings are used for each
  candidate. Common arguments include `n_starts`, `max_iter`, `missing`,
  `categorical`, `centering`, `scaling`, `drop_empty` and `boundary`.
  The group-class families accept only `n_starts`, `max_iter`, `tol` and
  `min_variance` through `...`.

- cores:

  Positive integer maximum number of parallel workers; the default, `1`,
  fits models serially. With a fixed `seed`, parallel and serial
  execution produce the same results. Parallel fitting requires
  additional memory. See *Parallel estimation* for platform requirements
  and the handling of interrupted workers.

- type:

  Optional model type: `"lca"`, `"lpa"` or `"lta"`. `NULL` (the default)
  determines the model from the other arguments. `"lca"` treats every
  indicator as categorical, as
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  does: indicators named in `ordinal` or `count` keep those
  distributions, and `categorical` adds to the set rather than replacing
  it. Numeric response codes are categorical with `type = "lca"`.
  `"lpa"` requires continuous indicators and rejects non-empty
  `categorical`, `ordinal` or `count`. Both require `time = NULL`.
  `"lta"` requires `time` and takes its measurement settings from
  `categorical`, as
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) does; it
  does not make numeric indicators categorical. An explicit type
  requires `family = "profiles"`; leave `type = NULL` for group-class
  families.

- verbose:

  Whether to report progress: one message each when the candidates are
  prepared, when their starts are fitted, and when the grid is done,
  with the elapsed time. The default reports when `cores > 1`, where a
  grid is usually long. Group-class families (`family`) are not
  reported.

## Value

An object of class `multilpa_enumeration` (or the transition or family
enumeration described above). Its `type` records the resolved model:
`"lca"`, `"lpa"` or `"lta"`, or `"mixed"`/`"family"` for those grids
under `type = NULL`. The print header reports this type.
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
extracts a fitted candidate, and
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
compares candidates selected by name.

## Model selection

Information criteria summarize relative support for the candidate
specifications, whereas entropy and posterior membership probabilities
describe classification uncertainty within each fitted model. These
measures address different aspects of the analysis. Model selection
should also examine whether the selected likelihood is replicated across
starts, whether estimates approach parameter bounds, and whether
additional classes yield interpretable differences in response patterns.
The function reports this evidence without prescribing a final number of
classes.

Comparisons of class counts depend on the measurement specification. For
example, additional profiles can represent residual associations that a
diagonal covariance model cannot accommodate within profiles. Evaluate
the number of classes alongside plausible covariance structures and
indicator distributions. The candidate range and fitting controls are
supplied by the analyst; the function does not expand the range
automatically.

## Parallel estimation

Socket workers are used on Windows and on macOS with Apple's Accelerate
BLAS; these workers load the installed package, which must match the
version used in the main R session. Other supported Unix configurations
use forked processes. Each worker requires memory for its fitting tasks.

Depending on the model family, parallel tasks consist of individual
starts or complete candidate models. The total number of workers does
not exceed `cores`. Warnings and messages are reported in candidate
order. Without a fixed `seed`, starting values can differ between runs.

Lost worker results are retried once with fewer workers, and a
`latents_worker_lost` warning reports the interruption. Unrecovered
tasks are recorded as failed starts or candidates; other results are
retained.

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
group-class family or with `between_variance`; and when `cores` is not
an integer of at least 1; and when `type` is invalid or conflicts with
`time`, `family` or the requested measurement distributions.

## References

Celeux, G., & Govaert, G. (1995). Gaussian parsimonious clustering
models. *Pattern Recognition*, 28(5), 781–793.

## See also

[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
to extract a fitted candidate,
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
#> 2         266.9458    274.3684        283.5788   367.5789       385.9996
#> 3         264.2039    267.9152        272.5204   276.5204       285.7307
#> 4         268.4185    277.6968        289.2098   426.5893       449.6151
#>   icl_groups icl_individual clc_groups clc_individual profile_entropy
#> 1   265.9152       270.5204   261.3100       261.3100              NA
#> 2   346.3686       355.5789   337.1583       337.1583      0.45177478
#> 3   265.9152       270.5204   261.3100       261.3100              NA
#> 4   400.0763       411.5893   388.5634       388.5634      0.08115115
#>   group_entropy avepp_min avepp_max avepp_mean avepp_below group_avepp_min
#> 1            NA        NA        NA         NA          NA              NA
#> 2            NA 0.6411851 0.8529074  0.7470462           1              NA
#> 3            NA        NA        NA         NA          NA              NA
#> 4            NA 0.6380736 0.6405599  0.6393168           2              NA
#>   group_avepp_max group_avepp_mean group_avepp_below converged boundary
#> 1              NA               NA                NA      TRUE    FALSE
#> 2              NA               NA                NA      TRUE    FALSE
#> 3              NA               NA                NA      TRUE    FALSE
#> 4              NA               NA                NA      TRUE    FALSE
#>   min_standardized_eigenvalue min_gaussian_mass n_best_replicated warnings
#> 1                   1.0000000         100.00000                 2         
#> 2                   0.8400023          17.67382                 2         
#> 3                   1.0000000         100.00000                 2         
#> 4                   0.8201279          44.58165                 2         
#>   error
#> 1  <NA>
#> 2  <NA>
#> 3  <NA>
#> 4  <NA>
summary(candidates)
#> Class enumeration: 4 candidates, 4 converged, 0 failed to fit
#> 2 distinct candidate(s) are minimal under some criterion.
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
#>       270.4          279.6          355.6         0.45177    0.6412     0.7470
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
#> -- criteria --------------------------------------------------------
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
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> -- classification --------------------------------------------------
#>              candidate n_profiles n_group_classes model       level class
#>  profiles2_groups1_EEI          2               1   EEI individuals     1
#>  profiles2_groups1_EEI          2               1   EEI individuals     2
#>  profiles2_groups1_VVI          2               1   VVI individuals     1
#>  profiles2_groups1_VVI          2               1   VVI individuals     2
#>  n_modal proportion_modal estimated_n estimated_proportion average_posterior
#>        6             0.06       17.67               0.1767            0.6412
#>       94             0.94       82.33               0.8233            0.8529
#>       31             0.31       44.58               0.4458            0.6381
#>       69             0.69       55.42               0.5542            0.6406
#>  odds_correct_classification relative_entropy
#>                        8.324          0.45177
#>                        1.245          0.45177
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
#>  profiles2_groups1_EEI        1   0.1767 0.1767        0.1767           1
#>  profiles2_groups1_EEI        2   0.8233 0.8233        0.8233           1
#>  profiles2_groups1_VVI        1   0.4458 0.4458        0.4458           1
#>  profiles2_groups1_VVI        2   0.5542 0.5542        0.5542           1
#> 
#> 7 tables above, truncated to fit. get_results(x, what = ) returns any
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
#> -- criteria --------------------------------------------------------
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
#>    ... 4 more rows.  get_results(x, what = "criteria")
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
#>                    157977.78           0.9999
#>                          Inf           0.9999
#>                          Inf           0.9704
#>                       129.97           0.9704
#>                       138.39           0.9704
#>                    329661.62           0.9999
#>                    147652.77           0.9999
#>                          Inf           0.9616
#>                       112.24           0.9616
#>                        74.41           0.9616
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
#>    0.2992 2.275e-17     6.824e-17   7.603e-17
#>    0.2992 2.992e-01     4.488e-01   1.000e+00
#>    0.3675 2.015e-21     6.044e-21   5.482e-21
#>    0.3675 3.675e-01     5.512e-01   1.000e+00
#>    0.3333 3.333e-01     3.333e-01   1.000e+00
#>    0.6667 6.667e-01     6.667e-01   1.000e+00
#>    ... 6 more rows.  get_results(x, what = "tree")
#> 
#> 7 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
