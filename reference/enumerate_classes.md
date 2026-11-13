# Class enumeration and comparison of mixture models

`enumerate_classes()` fits every candidate in a grid of class counts
(and, for continuous indicators, covariance structures) and collects
their information criteria, classification quality and estimation
diagnostics in one object. It covers latent profile, latent class and
latent transition models, single-level or two-level, and the additive,
dispersion and additive-dispersion families of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
The fitted candidates stay in the result (see *Value*).

Information criteria summarise relative support; entropy and average
posteriors describe classification certainty within each candidate.
Neither picks a solution: also check that the best likelihood was
replicated across starts, that no estimate sits at a bound, and that an
added class brings interpretable differences in the response patterns.
More profiles can stand in for covariances a diagonal structure cannot
hold, so compare class counts alongside plausible structures. The grid
is the analyst's: when every criterion is lowest at its edge, extend it.

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

  A data frame with one row per observation: the indicators, the group
  identifier and, when used, the occasion column and covariates. Use the
  same rows and scale for every candidate.

- vars:

  Names of the indicator columns.

- id:

  The column identifying the groups the rows are nested in. `NULL`
  enumerates single-level models, with a message;
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  and
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  do that directly.

- n_profiles:

  Profile (class) counts to try (default `1:4`; `2:4` for transition
  models).

- n_group_classes:

  Group-class counts to try (default `1:3`; default `1` with `time`, and
  only `1` with `id = NULL`).

- model:

  Covariance structures to cross with the class counts. `"basic"` (the
  default for profile grids) fits EEI (equal variances, no covariances),
  VVI (varying variances), EEE (one covariance matrix shared by all
  profiles) and VVV (a covariance matrix per profile). `"all"` fits all
  fourteen structures. Otherwise mclust's three-letter codes (volume,
  shape and orientation, each equal `E` or varying `V`, `I` for no
  covariances), or `NULL` for none, which is the default for a
  transition grid. With one continuous indicator, codes that coincide
  are fitted once; with only categorical indicators no structure is
  crossed. Giving `variance_model`, `covariance_model`, `volume`,
  `shape` or `orientation` in `...` sets one structure instead.

- seed:

  Optional integer seed, reproducing the starting values for the same
  data and controls.

- family:

  `"profiles"` (default) enumerates the profile model over `n_profiles`,
  `n_group_classes` and `model`. One or more of `"additive"`,
  `"dispersion"` and `"additive_dispersion"` enumerates those
  group-class families over `family`, `n_group_classes` and
  `between_variance` instead; they need `id`, refuse `n_profiles` and
  `model`, and return a `latents_family_enumeration`.

- between_variance:

  For the group-class families, the restrictions to cross: `"varying"`,
  `"equal"` or both (default). The dispersion family is fitted with
  `"equal"` only.

- time:

  The occasion column. Enumerates latent transition models
  ([`lta()`](https://pak.dynasite.org/latents/reference/lta.md)) over
  `n_profiles` and `n_group_classes`, every other
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) argument
  in `...` held fixed; `model`, when given, crosses structures. Profile
  family only; returns a `latents_transition_enumeration`.

- ...:

  Further arguments of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  or of [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  with `time`, used for every candidate: typically `n_starts`,
  `max_iter`, `missing`, `categorical`, `centering`, `scaling`,
  `drop_empty` and `boundary`. The group-class families take only
  `n_starts`, `max_iter`, `tol` and `min_variance`.

- cores:

  Maximum number of parallel workers (default 1, serial). With `seed`,
  parallel and serial runs give the same result. On Windows, and on
  macOS with Apple's Accelerate BLAS, each worker is a separate R
  process, so memory grows with `cores`.

- type:

  `"lca"`, `"lpa"` or `"lta"`; `NULL` (default) reads the type from the
  other arguments. `"lca"` treats every indicator as categorical,
  numeric codes included, as
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  does (indicators in `ordinal` or `count` keep those distributions, and
  `categorical` adds to the set). `"lpa"` needs continuous indicators
  and no `categorical`, `ordinal` or `count`. Both need `time = NULL`.
  `"lta"` needs `time` and reads the measurement from `categorical`, as
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) does,
  without making numeric indicators categorical. An explicit type needs
  `family = "profiles"`.

- verbose:

  Whether to report progress: one message each when the candidates are
  prepared, when starts are screened, when they are fitted and when the
  grid is done. Default on when `cores > 1`. Not reported for the
  group-class families.

## Value

A `multilpa_enumeration` object, or the transition or family enumeration
named under `time` and `family`; the print header names the resolved
type (LPA, LCA, LTA, mixed or family). On every grid:
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
gives one row per candidate (named like `"profiles3_groups2"`) with its
class counts, log likelihood, criteria, entropy, convergence, boundary
and replication, leaving out columns empty in every row;
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
the best candidate per criterion;
[`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md)
and
[`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md)
how classes split;
[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md),
[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md)
and
[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md)
the classification tables of every candidate;
[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)
the criteria;
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
one fitted candidate;
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
candidates by name (a family grid is compared in its own candidates
table). A profile or class grid also has the AvePP summaries at each
level in its candidates table,
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
and
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md),
and a [`summary()`](https://rdrr.io/r/base/summary.html) that prints one
screen of the evidence; a transition or family grid's
[`summary()`](https://rdrr.io/r/base/summary.html) is its candidates
table. Invalid arguments are refused with the classed conditions listed
on
[latents-conditions](https://pak.dynasite.org/latents/reference/latents-conditions.md).

## See also

[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md),
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
to test one class count against the next,
[`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
for mixture regressions.

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
#> 
#> ! 2 candidates have a class with average posterior below 0.8
#>     (profiles2_groups1_EEI, profiles2_groups1_VVI).
#> 
#> -- lowest value of each criterion --------------------------------
#>  criterion  convention             candidate value
#>        aic        <NA> profiles1_groups1_EEI 265.3
#>        kic        <NA> profiles1_groups1_EEI 270.3
#>        bic      groups profiles1_groups1_EEI 265.9
#>        bic individuals profiles1_groups1_EEI 270.5
#>      sabic      groups profiles2_groups1_VVI 257.7
#>      sabic individuals profiles1_groups1_EEI 264.2
#>       caic      groups profiles1_groups1_EEI 267.9
#>       caic individuals profiles1_groups1_EEI 272.5
#>        awe      groups profiles1_groups1_EEI 276.5
#>        awe individuals profiles1_groups1_EEI 285.7
#>        icl      groups profiles1_groups1_EEI 265.9
#>        icl individuals profiles1_groups1_EEI 270.5
#>        clc      groups profiles1_groups1_EEI 261.3
#>        clc individuals profiles1_groups1_EEI 261.3
#> 
#> -- knees of bic_groups ----------------------------------------
#>   direction                line knee_candidate knee_without_largest stable
#>    profiles  1 group class, EEI           <NA>                 <NA>     NA
#>    profiles  1 group class, VVI           <NA>                 <NA>     NA
#>  parameters all candidates, EEI           <NA>                 <NA>     NA
#>  parameters all candidates, VVI           <NA>                 <NA>     NA
#> 
#> Tables: get_candidates(x), get_best(x), get_deltas(x), get_knees(x),
#>         get_tree(x), get_classification(x), get_entropy(x), get_checks(x).
#> Plots:  plot_criteria(x), plot_deltas(x), plot_avepp(x), plot_tree(x), plot_all(x).

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
#> 
#> ! 1 of 6 candidates reached their best log likelihood from one start only
#>     (profiles3_VVV); refit those you consider with more starts.
#> 
#> -- lowest value of each criterion --------------------------------
#>  criterion  convention     candidate value
#>        aic        <NA> profiles3_VVV 448.4
#>        kic        <NA> profiles3_VVV 495.4
#>        bic      groups profiles2_VVV 574.0
#>        bic individuals profiles2_VVV 574.0
#>      sabic      groups profiles3_VVV 441.6
#>      sabic individuals profiles3_VVV 441.6
#>       caic      groups profiles2_VVV 603.0
#>       caic individuals profiles2_VVV 603.0
#>        awe      groups profiles2_VVV 806.3
#>        awe individuals profiles2_VVV 806.3
#>        icl      groups profiles2_VVV 574.0
#>        icl individuals profiles2_VVV 574.0
#>        clc      groups profiles3_VVV 370.1
#>        clc individuals profiles3_VVV 370.1
#> 
#> -- knees of bic_groups ----------------------------------------
#>   direction                line knee_candidate knee_without_largest stable
#>    profiles  1 group class, EEE  profiles2_EEE                 <NA>     NA
#>    profiles  1 group class, VVV  profiles2_VVV                 <NA>     NA
#>  parameters all candidates, EEE  profiles2_EEE                 <NA>     NA
#>  parameters all candidates, VVV  profiles2_VVV                 <NA>     NA
#> 
#> Tables: get_candidates(x), get_best(x), get_deltas(x), get_knees(x),
#>         get_tree(x), get_classification(x), get_entropy(x), get_checks(x).
#> Plots:  plot_criteria(x), plot_deltas(x), plot_avepp(x), plot_tree(x), plot_all(x).
```
