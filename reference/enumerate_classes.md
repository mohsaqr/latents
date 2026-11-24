# Fit models with different numbers of classes and collect their criteria and diagnostics for comparison

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
  models, which need at least two profiles).

- n_group_classes:

  Group-class counts to try (default `1:3`; default `1` with `time`, and
  only `1` with `id = NULL`). Every count is crossed with every value of
  `n_profiles`. A candidate with one profile and more than one group
  class is not identified: it is recorded as failed, with the reason in
  the `error` column.

- model:

  The covariance structures of the continuous indicators to fit at every
  combination of class counts:

  - `"basic"` (default for profile models): four structures. `"EEI"`:
    the profiles share their variances and the indicators are
    uncorrelated within a profile. `"VVI"`: each profile has its own
    variances, still uncorrelated. `"EEE"`: the profiles share one
    covariance matrix. `"VVV"`: each profile has its own covariance
    matrix.

  - `"all"`: all fourteen structures of mclust.

  - One or more of the codes `"EII"`, `"VII"`, `"EEI"`, `"VEI"`,
    `"EVI"`, `"VVI"`, `"EEE"`, `"VEE"`, `"EVE"`, `"VVE"`, `"EEV"`,
    `"VEV"`, `"EVV"`, `"VVV"`: those structures only. The three letters
    are the volume, shape and orientation, each equal (`E`) or varying
    (`V`) across profiles, with `I` for a spherical shape or no
    covariances;
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    explains them.

  - `NULL`: no structure is crossed, and every candidate has the
    structure
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    fits by default (`"VVI"`). This is the default with `time`.

  With `time`, only codes and `NULL` are accepted (`"basic"` and `"all"`
  are not). With one continuous indicator, codes that coincide are
  fitted once. With no continuous indicator there is no structure to
  cross, and an explicit `model` is an error. Giving `variance_model`,
  `covariance_model`, `volume`, `shape` or `orientation` in `...` sets
  one structure for every candidate; `model` must then be left out.

- seed:

  Optional integer seed, reproducing the starting values for the same
  data and controls.

- family:

  Which model family to enumerate:

  - `"profiles"` (default): the profile model, over `n_profiles`,
    `n_group_classes` and `model`;

  - one or more of `"additive"`, `"dispersion"` and
    `"additive_dispersion"`: those group-class families of
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    over `family`, `n_group_classes` and `between_variance`. They need
    `id`, refuse `n_profiles`, `model` and `time`, and return a
    `latents_family_enumeration`.

  The cross-level families of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  are not enumerated.

- between_variance:

  For the group-class families only, the restrictions on the
  between-group variances to cross:

  - `"varying"`: each group class has its own between-group variances;

  - `"equal"`: all group classes share them.

  The default, `c("varying", "equal")`, fits both. The dispersion family
  is fitted with `"equal"` only. Giving the argument with
  `family = "profiles"` is an error.

- time:

  The occasion column. Enumerates latent transition models
  ([`lta()`](https://pak.dynasite.org/latents/reference/lta.md)) over
  `n_profiles` and `n_group_classes`, every other
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) argument
  in `...` held fixed; `model`, when given, crosses structures. Profile
  family only; returns a `latents_transition_enumeration`. `NULL`
  (default) enumerates profile or class models.

- ...:

  Further arguments of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  or of [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  with `time`, used for every candidate: typically `n_starts`,
  `max_iter`, `missing` (`"error"` or `"fiml"`), `categorical`,
  `centering`, `scaling`, `drop_empty` and `boundary` (`"allow"` or
  `"exclude"`). The group-class families take only `n_starts`,
  `max_iter`, `tol` and `min_variance`. A name that is not an argument
  of the fitting function is refused.

- cores:

  Maximum number of parallel workers (default 1, serial). With `seed`,
  parallel and serial runs give the same result. On Windows, and on
  macOS with Apple's Accelerate BLAS, each worker is a separate R
  process, so memory grows with `cores`.

- type:

  The kind of model, stated explicitly:

  - `NULL` (default): read from the other arguments. `time` gives
    transition models; otherwise the indicators named in `categorical`,
    `ordinal` and `count` decide between profile, class and mixed
    models;

  - `"lpa"`: latent profile models. Every indicator is continuous;
    `categorical`, `ordinal` and `count` are refused;

  - `"lca"`: latent class models. Every indicator is categorical,
    numeric codes included, as
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    does (indicators in `ordinal` or `count` keep those distributions,
    and `categorical` adds to the set);

  - `"lta"`: latent transition models. Needs `time`, and reads the
    measurement from `categorical`, as
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) does,
    without making numeric indicators categorical.

  `"lpa"` and `"lca"` need `time = NULL`, and an explicit type needs
  `family = "profiles"`.

- verbose:

  `TRUE` reports progress with one message each when the candidates are
  prepared, when starts are screened, when they are fitted and when the
  grid is done; `FALSE` is silent. The default is `TRUE` when
  `cores > 1` and `FALSE` otherwise. Not reported for the group-class
  families.

## Value

An enumeration grid that keeps every fitted candidate. Its class depends
on what was enumerated:

- `multilpa_enumeration`: profile, class or mixed models (the default);

- `latents_transition_enumeration`: transition models (with `time`);

- `latents_family_enumeration`: group-class families (with `family`).

The print header names the resolved type (LPA, LCA, LTA, mixed or
family).

On every grid,
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
(also [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html))
gives one row per candidate, named like `"profiles3_groups2"`: with the
structure appended, as in `"profiles3_groups2_VVI"`, when several are
crossed; with `classes` for `profiles` when every indicator is
categorical; without the `groups` part on a single-level grid
(`"profiles3"`); and by family on a family grid (`"additive_groups2"`).
Its columns:

- `n_profiles`, `n_group_classes`, `model`, `log_likelihood` and
  `n_parameters`;

- the information criteria `aic`, `kic`, `bic`, `sabic`, `caic`, `awe`,
  `icl` and `clc`. On a profile or class grid each except `aic` and
  `kic` has two columns: `_groups` uses the number of groups as the
  sample size and `_individual` the number of rows. Lower is better. A
  transition or family grid has `aic`, `bic` (groups) and
  `bic_individual`;

- `profile_entropy` and `group_entropy` (`entropy` on a transition or
  family grid): relative entropy between 0 and 1, higher meaning clearer
  classification;

- on a profile or class grid, the average posterior probability of the
  assigned class (`avepp_min`, `avepp_max`, `avepp_mean`, and
  `avepp_below`, the number of classes under 0.8), with the same for
  group classes as `group_avepp_*`;

- `converged`, `boundary` (an estimate at its bound) and
  `n_best_replicated` (the number of starts that reached the best log
  likelihood);

- `error`: why a candidate could not be fitted. Such a candidate keeps
  its row, with `NA` values. A transition or family grid also has
  `warnings`, the warnings of each candidate.

Columns empty in every row are left out. The other accessors:
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
gives the best candidate per criterion;
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
table). A profile or class grid also has
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
and
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md),
and a [`summary()`](https://rdrr.io/r/base/summary.html) that prints one
screen of the evidence; a transition or family grid's
[`summary()`](https://rdrr.io/r/base/summary.html) is its candidates
table. No candidate is selected automatically.

## Conditions

`latents_bad_argument` for an argument value or combination that is not
accepted: a `type` that conflicts with `time`, `family` or the
measurement arguments; `model` given with another structure argument,
with no continuous indicator, or naming an unknown structure;
`n_profiles` or `model` with a group-class family; `between_variance`
with the profile family; more than one group class with `id = NULL`; a
name in `...` that the fitting function does not have.
`latents_bad_transition` for a transition grid with fewer than two
profiles. A candidate that fails for a reason of its own stays in the
grid with its `error`. When every candidate fails for the same reason in
the data or the arguments (a column that does not exist, for example),
that error is raised instead: `latents_bad_data`,
`latents_bad_categorical`, `latents_missing_data`, `latents_bad_weights`
or `latents_bad_argument`. The candidates' warnings are raised once the
grid is fitted. When a converged candidate fitted from two or more
starts reached its best log likelihood from one of them only, a
`latents_single_start` message names it once per grid.
[latents-conditions](https://pak.dynasite.org/latents/reference/latents-conditions.md)
lists every condition.

## See also

[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md),
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
to test one class count against the next,
[`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
for mixture regressions.

## Examples

``` r
# Course sessions nested in students: numbers of session profiles and of
# student classes.
activity <- c("browse", "lectures")
candidates <- enumerate_classes(course_engagement, activity, id = "student",
                                n_profiles = 2:3, n_group_classes = 1:2,
                                model = "VVI", n_starts = 2, seed = 1)
#> 2 of 4 candidates reached their best log likelihood from one start only (profiles3_groups1, profiles3_groups2); refit those you consider with more starts.
candidates
#> Class enumeration: 4 candidate models (LPA)
#>          candidate n_profiles n_group_classes model log_likelihood n_parameters
#>  profiles2_groups1          2               1   VVI       -3894.21            9
#>  profiles3_groups1          3               1   VVI       -3889.10           14
#>  profiles2_groups2          2               2   VVI       -3799.97           11
#>  profiles3_groups2          3               2   VVI       -3795.95           17
#>      aic bic_groups bic_individual icl_individual profile_entropy group_entropy
#>  7806.43    7830.40        7853.77        8701.82           0.570            NA
#>  7806.20    7843.49        7879.84        8846.75           0.691            NA
#>  7621.94    7651.24        7679.80        8374.77           0.642         0.860
#>  7625.89    7671.17        7715.31        8759.53           0.663         0.863
#>  avepp_min avepp_max avepp_mean avepp_below group_avepp_min group_avepp_max
#>      0.848     0.886      0.867           0              NA              NA
#>      0.507     0.895      0.737           1              NA              NA
#>      0.882     0.905      0.893           0           0.942           0.968
#>      0.624     0.904      0.779           1           0.949           0.965
#>  group_avepp_mean group_avepp_below converged boundary n_best_replicated
#>                NA                NA      TRUE    FALSE                 2
#>                NA                NA      TRUE    FALSE                 1
#>             0.955                 0      TRUE    FALSE                 2
#>             0.957                 0      TRUE    FALSE                 1
#> Every column, unrounded: as.data.frame(x).
#> No candidate is selected automatically. Compare one convention consistently.
get_best(candidates)
#>    criterion  convention         candidate n_profiles n_group_classes model
#> 1        aic        <NA> profiles2_groups2          2               2   VVI
#> 2        kic        <NA> profiles2_groups2          2               2   VVI
#> 3        bic      groups profiles2_groups2          2               2   VVI
#> 4        bic individuals profiles2_groups2          2               2   VVI
#> 5      sabic      groups profiles2_groups2          2               2   VVI
#> 6      sabic individuals profiles2_groups2          2               2   VVI
#> 7       caic      groups profiles2_groups2          2               2   VVI
#> 8       caic individuals profiles2_groups2          2               2   VVI
#> 9        awe      groups profiles2_groups2          2               2   VVI
#> 10       awe individuals profiles2_groups2          2               2   VVI
#> 11       icl      groups profiles2_groups2          2               2   VVI
#> 12       icl individuals profiles2_groups2          2               2   VVI
#> 13       clc      groups profiles2_groups2          2               2   VVI
#> 14       clc individuals profiles2_groups2          2               2   VVI
#>       value
#> 1  7621.938
#> 2  7635.938
#> 3  7651.236
#> 4  7679.796
#> 5  7616.483
#> 6  7644.853
#> 7  7662.236
#> 8  7690.796
#> 9  8430.504
#> 10 8487.624
#> 11 8346.206
#> 12 8374.766
#> 13 8294.908
#> 14 8294.908
get_candidate(candidates, "profiles2_groups2")
#> Two-level latent profile analysis: 2 profiles, 2 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -3799.969233 | AIC: 7621.938 | BIC (groups): 7651.236
#> Converged: TRUE | iterations: 24 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.5058 (profiles, profile_2); Gaussian mass: 580.3; weight effective n: 712.4.
#> 
#> Class sizes   profile_1 59.2% (842)   profile_2 40.8% (580)
#> Group classes group_class_1 32.4% (34)   group_class_2 67.6% (72)
#> 
#> Means
#>  indicator profile_1 profile_2
#>     browse     0.549    -0.797
#>   lectures     0.405    -0.587
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").

# Latent transition models over each student's sequence of sessions:
moves <- enumerate_classes(course_engagement,
                           c(activity, "forum_read", "forum_post",
                             "attendance"),
                           id = "student", time = "sequence",
                           n_profiles = 2:3, n_starts = 2, seed = 1)
get_candidates(moves)
#>           candidate n_profiles n_group_classes log_likelihood n_parameters
#> 1 profiles2_groups1          2               1      -8335.287           23
#> 2 profiles3_groups1          3               1      -8261.976           38
#>        aic      bic bic_individual   entropy converged boundary
#> 1 16716.57 16777.83       16837.55 0.9611453      TRUE    FALSE
#> 2 16599.95 16701.16       16799.83 0.8005931      TRUE    FALSE
#>   n_best_replicated
#> 1                 2
#> 2                 2

# Latent class models, stated explicitly with `type`:
items <- c("time_with_friends", "on_social_media", "tv_video_games",
           "sports")
classes <- enumerate_classes(subset(student_esm, day <= 1), items,
                             id = "student", type = "lca",
                             n_profiles = 1:2, n_group_classes = 1,
                             n_starts = 2, seed = 1)
#> 1 of 2 candidates reached their best log likelihood from one start only (classes2_groups1); refit those you consider with more starts.
classes
#> Class enumeration: 2 candidate models (LCA)
#>         candidate n_profiles log_likelihood n_parameters     aic bic_groups
#>  classes1_groups1          1        -555.14            4 1118.28    1128.41
#>  classes2_groups1          2        -547.83            9 1113.66    1136.45
#>  bic_individual icl_individual profile_entropy avepp_min avepp_max avepp_mean
#>         1133.48        1133.48              NA        NA        NA         NA
#>         1147.85        1154.46           0.986     0.998         1      0.999
#>  avepp_below converged boundary n_best_replicated
#>           NA      TRUE    FALSE                 2
#>            0      TRUE    FALSE                 1
#> Every column, unrounded: as.data.frame(x).
#> No candidate is selected automatically. Compare one convention consistently.
```
