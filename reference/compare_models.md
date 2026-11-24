# Compare the fit, classification and convergence of several fitted models in one table

`compare_models()` returns one table of the fit, complexity,
classification and convergence of several models fitted to the same
data: separately fitted models, or candidates of an enumeration grid.
Each row holds one model's number of parameters, log likelihood,
information criteria, entropy, average posterior probabilities and
smallest class share. BIC differences give each model's distance from
the smallest BIC among those compared, and the Schwarz weights
`exp(-delta_bic / 2)`, normalised over them, give relative support. The
function does not select a model: a lower criterion favours a model
among those compared and says nothing about absolute fit.

## Usage

``` r
compare_models(
  ...,
  n_profiles = NULL,
  n_group_classes = NULL,
  cutoff = 0.8,
  test = c("none", "lrt")
)
```

## Arguments

- ...:

  The models to compare, in one of three forms:

  - an enumeration grid from
    [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
    [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
    or
    [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
    followed by the candidates to compare (none for all); see *Choosing
    candidates*;

  - two or more fits from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    (with or without membership covariates),
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
    [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md);

  - two or more fits from
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
    with or without `random`.

  The two kinds of fit cannot be mixed in one call. Named fits are
  labelled by their names, unnamed ones by their expression, candidates
  by their candidate name; labels must be distinct. Additive, dispersion
  and cross-level fits of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  are refused (`latents_bad_argument`).

- n_profiles, n_group_classes:

  With a grid, the class counts of the candidates to compare; `NULL`
  (default) takes every count. Candidates whose fit failed are left out,
  with a message. Without a grid they raise `latents_bad_argument`.

- cutoff:

  For profile and class models, the average posterior probability below
  which a class is counted in `avepp_below` and `group_avepp_below`
  (default 0.8; a number above 0 and at most 1). Not used for mixture
  regressions.

- test:

  Whether the models are also tested against each other:

  - `"none"` (default): the comparison table alone;

  - `"lrt"`: adds a likelihood-ratio test of each model against the one
    in the row before it (`lr_statistic`, `lr_df`, `lr_p_value`), for
    nested profile, class and transition models. See *Likelihood-ratio
    tests*.

## Value

A data frame of class `latents_comparison` (also a `latents_table`,
which prints as a formatted table), one row per model, in the order
given. [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
its information criteria, and `get_candidate(comparison, "<model>")`
returns the fitted model of a row.

For profile, class and transition models the columns are `model` (the
label), `n_profiles` (classes, for a class model), `n_group_classes`,
`structure` (the covariance structure, `NA` without continuous
indicators), `n_parameters` (see *Held measurement*), `log_likelihood`;
the information criteria `aic`, `kic`, and `bic`, `sabic`, `caic`,
`awe`, `icl` and `clc`, each of these six in two columns, `_groups` (the
number of groups as the sample size) and `_individual` (the number of
rows); `profile_entropy`, `group_entropy`; `avepp_min`, `avepp_max`,
`avepp_mean` (the smallest, largest and mean average posterior
probability of the profiles) and `avepp_below` (the number of profiles
below `cutoff`), and the same four prefixed `group_` for the group
classes; `smallest_share` (the smallest estimated profile share),
`converged`, `n_best_replicated` (the starts that reached the best log
likelihood), `delta_bic` and `bic_weight` (on `bic_groups`, the BIC that
counts groups as the independent units; the ordinary BIC for
single-level models).

For mixture regressions: `model`, `n_classes`, `random` (the
random-effects specification, `"none"` without one), `n_parameters`,
`log_likelihood`, `aic`, `bic`, `sabic`, `icl`, `entropy`,
`smallest_share`, `delta_bic`, `bic_weight` and `converged`.

With `test = "lrt"`, also `lr_statistic`, `lr_df` and `lr_p_value`, `NA`
in the first row.

Errors: `latents_bad_argument` for fewer than two models, models of
mixed or unsupported kinds, duplicated labels, a group-class family fit
or grid, a mixture-regression grid, or `test = "lrt"` on mixture
regressions (the family and regression grids have their own tables,
`as.data.frame(grid)`); `latents_incomparable_models` for fits to
different data; and, under `test = "lrt"`, `latents_bad_nesting` for
models that are not nested and `latents_singular_restriction` for a
singular restriction.
[latents-conditions](https://pak.dynasite.org/latents/reference/latents-conditions.md)
describes each class. The criteria compare the models with each other; a
final choice also needs the estimation to be stable and the classes to
mean something.

## Held measurement

Fits that hold the same measurement
([`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md),
`fixed`) are compared on the parameters they estimate. In any other
comparison, such as a staged fit against a jointly estimated one, every
criterion and `n_parameters` count the measurement parameters too, and
the table's title says so. Such a comparison is descriptive, since
staged estimates do not maximise the joint likelihood.

## Choosing candidates

Candidates of a grid are named as its `candidate` column shows them
(`compare_models(grid, "profiles3_groups4", "profiles3_groups5")`), by
class counts (`compare_models(grid, c(3, 4), c(3, 5))`, a single number
meaning one group class), or through `n_profiles` and `n_group_classes`
(`compare_models(grid, n_profiles = 3)` compares every candidate with
three profiles). Latent class grids work the same way with `classes` in
the names.

## Likelihood-ratio tests

With `test = "lrt"`, each model is tested against the model in the row
before it, as [`stats::anova()`](https://rdrr.io/r/stats/anova.html)
does: `lr_statistic` is twice the difference in log likelihoods, `lr_df`
the difference in free parameters and `lr_p_value` the upper chi-square
tail. List the models from the most restricted to the most general. Two
models are nested when they share indicators, coding, estimation
settings and numbers of profiles and group classes, and the first only
restricts the second: a covariance structure restricted letter by letter
(EEI within VVI or EEE, both within VVV), `local_dependence` pairs
within more pairs or a full covariance, fewer direct effects or
membership covariates, profile slopes shared across group classes, count
dispersions shared across profiles, and for
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits also
invariant within occasion-specific measurement, homogeneous within
occasion-specific transitions, first within second order, fewer
transition or initial covariates. A covariate both models use must enter
both with the same values (or an affine recoding).

Nesting is not enough: the restriction must also be regular. A spherical
structure (EII, VII) within EVE, VVE, EEV or VEV is singular, because
the orientation those structures add has no value when every profile is
spherical; such pairs are refused, so compare them on the criteria. A
different number of profiles or group classes, or a stayer class, puts
the null on the boundary of the parameter space, where the chi-square
reference does not hold either: those comparisons are refused and belong
to
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
(or, for class counts, `get_results(grid, "lmr")`). Fits with sampling
weights, a prior or a noise component, and mixture regressions, are
refused under `test = "lrt"`.

## See also

[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md);
[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)
draws the criteria of a comparison.

## Examples

``` r
# Compare trajectory and random-intercept growth models for 60 students.
few_students <- subset(growth_scores, student <= 60)
trajectories <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                                   id = "student", class_level = "group",
                                   seed = 1)
growth <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                             id = "student", class_level = "group",
                             random = "intercept", n_starts = 3, seed = 1)
compare_models(trajectories = trajectories, growth = growth)
#> Model comparison
#> 
#> Model         Classes  Random     k    LogLik      BIC    ΔBIC  Weight  Entropy
#> ------------  -------  ---------  -  --------  -------  ------  ------  -------
#> trajectories        2  none       7  -1277.94  2584.53  196.11   0.000     0.94
#> growth              2  intercept  9  -1175.79  2388.43    0.00   1.000     0.97

# Candidates of an enumeration grid, by name, by class counts, or selected
grid <- enumerate_classes(subset(course_engagement, student <= 30),
                          c("browse", "lectures"), "student",
                          n_profiles = 2:3, n_group_classes = 1:2,
                          model = "VVI", n_starts = 2, seed = 1)
#> 1 of 4 candidates reached their best log likelihood from one start only (profiles3_groups1); refit those you consider with more starts.
compare_models(grid, "profiles2_groups2", "profiles3_groups2")
#> Model comparison
#> 
#> Model              Profiles  Group classes   k    LogLik      BIC   ΔBIC  Weight      ICL  Entropy  Min AvePP  < 0.8  Group entropy  Group min AvePP
#> -----------------  --------  -------------  --  --------  -------  -----  ------  -------  -------  ---------  -----  -------------  ---------------
#> profiles2_groups2         2              2  11  -1097.18  2231.78   0.00   0.996  2450.46     0.67       0.90      0           0.97             0.99
#> profiles3_groups2         3              2  17  -1092.40  2242.62  10.84   0.004  2578.50     0.68       0.70      1           0.96             0.99
compare_models(grid, c(2, 2), c(3, 2))
#> Model comparison
#> 
#> Model              Profiles  Group classes   k    LogLik      BIC   ΔBIC  Weight      ICL  Entropy  Min AvePP  < 0.8  Group entropy  Group min AvePP
#> -----------------  --------  -------------  --  --------  -------  -----  ------  -------  -------  ---------  -----  -------------  ---------------
#> profiles2_groups2         2              2  11  -1097.18  2231.78   0.00   0.996  2450.46     0.67       0.90      0           0.97             0.99
#> profiles3_groups2         3              2  17  -1092.40  2242.62  10.84   0.004  2578.50     0.68       0.70      1           0.96             0.99
compare_models(grid, n_group_classes = 2)
#> Model comparison
#> 
#> Model              Profiles  Group classes   k    LogLik      BIC   ΔBIC  Weight      ICL  Entropy  Min AvePP  < 0.8  Group entropy  Group min AvePP
#> -----------------  --------  -------------  --  --------  -------  -----  ------  -------  -------  ---------  -----  -------------  ---------------
#> profiles2_groups2         2              2  11  -1097.18  2231.78   0.00   0.996  2450.46     0.67       0.90      0           0.97             0.99
#> profiles3_groups2         3              2  17  -1092.40  2242.62  10.84   0.004  2578.50     0.68       0.70      1           0.96             0.99

# Is a covariance structure that varies by profile needed? A
# likelihood-ratio test of each structure against the one before it.
shared <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
              model = "EEI", n_starts = 2, seed = 1)
varying <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
               model = "VVI", n_starts = 2, seed = 1)
compare_models(shared = shared, varying = varying, test = "lrt")
#> Model comparison
#> 
#> Model    Profiles  k    LogLik      BIC  ΔBIC  Weight      ICL  Entropy  Min AvePP  < 0.8     LR  df      p
#> -------  --------  -  --------  -------  ----  ------  -------  -------  ---------  -----  -----  --  -----
#> shared          2  7  -1157.68  2354.60  8.54   0.014  2357.27     0.99       1.00      0                  
#> varying         2  9  -1147.81  2346.06  0.00   0.986  2346.52     1.00       1.00      0  19.75   2  <.001
```
