# Compare the fit and classification of mixture models

`compare_models()` compares alternative mixture models fitted to the
same observed data. It supports assessment of the number of latent
classes and other model specifications by presenting evidence about
model fit, complexity, classification, and convergence.

Information criteria summarise the trade-off between fit and the number
of estimated parameters. Classification measures provide complementary
information about the separation and relative size of the classes, while
convergence summaries help identify estimation problems that may
compromise the comparison. These measures address different aspects of a
solution and should be interpreted together.

The function accepts separately fitted models or candidates from an
enumeration grid and returns their comparison statistics in a single
table. BIC differences express each model's distance from the smallest
BIC among the candidates; normalised BIC weights summarise their
relative support within that set. The function does not select a final
model.

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

  Either an enumeration grid from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  followed by the candidates to compare (none for all of them, or a
  selection through `n_profiles` and `n_group_classes`); or two or more
  fits from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md);
  or two or more fits from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  (with or without `random`). Named fits are labelled by their names,
  unnamed fits by their expression, and candidates by their candidate
  name. Labels must be distinct.

- n_profiles, n_group_classes:

  With a grid, the class counts of the candidates to compare; `NULL`
  (the default) takes every count. Candidates whose fit failed are left
  out, with a message.

- cutoff:

  For profile and class models, the average posterior probability below
  which a class is counted in `avepp_below` and `group_avepp_below`; a
  single number in (0, 1\], default 0.8. Not used for mixture
  regressions.

- test:

  `"none"` (the default) for the comparison table alone, or `"lrt"` to
  add a likelihood-ratio test of each profile, class or transition model
  against the one before it (see *Details*).

## Value

A base `data.frame` of class `latents_comparison`, one row per model.

For mixture regressions the columns are `model`, `n_classes`, `random`,
`n_parameters`, `log_likelihood`, `aic`, `bic`, `sabic`, `icl`,
`entropy`, `smallest_share`, `delta_bic`, `bic_weight` and `converged`.

For profile and class models they are `model`, `n_profiles`,
`n_group_classes`, `structure` (the covariance structure, `NA` without
continuous indicators), `n_parameters`, `log_likelihood`, every
information criterion of
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)'s
candidates table, `profile_entropy`, `group_entropy`, `avepp_min`,
`avepp_max`, `avepp_mean`, `avepp_below` and the same four prefixed
`group_`, `smallest_share` (the smallest estimated profile share,
relative to all observation mass including any noise component),
`converged`, `n_best_replicated`, and `delta_bic` and `bic_weight`,
computed on `bic_groups` (the BIC that counts groups as the independent
units; on a single-level model, the ordinary BIC). The values are those
of the grid's candidates table.

With `test = "lrt"`, three more columns: `lr_statistic`, `lr_df` and
`lr_p_value`, `NA` in the first row.

Raises `latents_bad_argument` for fits that cannot be compared, for
fewer than two models, for duplicate labels, for a group-class family
grid (its own table, `as.data.frame(grid)`, compares its candidates),
for candidates named together with `n_profiles` or `n_group_classes`,
and for `n_profiles` or `n_group_classes` without a grid;
`latents_unknown_candidate` for a name, counts or selection not in the
grid; `latents_failed_candidate` for a named candidate that could not be
fitted; `latents_incomparable_models` for fits to different data. A
`cutoff` outside (0, 1\] fails the argument checks.

With `test = "lrt"`: `latents_bad_nesting` for consecutive models that
are not nested as described in *Details*, or that differ in their
numbers of profiles or group classes; `latents_reversed_likelihood` when
the more general model has the lower likelihood beyond the fits'
tolerance; `latents_no_converge` for a fit that did not converge;
`latents_unsupported_weights`, `latents_unsupported_prior` and
`latents_unsupported_noise` for fits with sampling weights, a prior or a
noise component; `latents_bad_argument` for mixture regressions. A fit
with a variance at `min_variance` warns `latents_boundary`.

## Details

`delta_bic` is the difference from the smallest BIC among the supplied
candidates. `bic_weight` is the Schwarz weight `exp(-delta_bic / 2)`,
normalised across those candidates (Wagenmakers and Farrell, 2004).
Incompatible data produce an error.

Candidates of an enumeration grid are selected by name, as the grid's
`candidate` column shows them
(`compare_models(grid, "profiles3_groups4", "profiles3_groups5")`), or
by class counts (`compare_models(grid, c(3, 4), c(3, 5))`, where a
single number means one group class). Alternatively, `n_profiles` and
`n_group_classes` select candidates:
`compare_models(grid, n_profiles = 3)` compares every candidate with 3
profiles, and `compare_models(grid)` compares them all. Latent class
grids work the same way, with `classes` in the names, and the
`n_profiles` column counting classes. Fits from
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
[`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
can also be passed directly.

With `test = "lrt"`, each model is tested against the model in the row
before it, as [`stats::anova()`](https://rdrr.io/r/stats/anova.html)
does: `lr_statistic` is twice the difference in their log likelihoods,
`lr_df` the difference in their numbers of free parameters, and
`lr_p_value` its upper chi-square tail (Wilks, 1938). List the models
from the most restricted to the most general. A model is nested in the
next when the two have the same indicators, coding, estimation settings
and numbers of profiles and group classes, and the first only restricts
the second: a covariance structure restricted letter by letter (EEI
within VVI or EEE, and both within VVV), covariances of chosen pairs
(`local_dependence`) within more pairs or a full covariance, fewer joint
item pairs, direct effects or membership covariates, profile slopes
shared across group classes, count dispersions shared across profiles;
and for [`lta()`](https://pak.dynasite.org/latents/reference/lta.md),
measurement invariant within measurement by occasion, homogeneous within
occasion-specific transitions, first within second order, and fewer
transition or initial covariates. These restrictions hold parameters at
interior values, where the chi-square reference applies. A different
number of profiles or group classes, or a stayer class, puts the null on
the boundary of the parameter space, where it does not: those
comparisons are refused and belong to
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
or `get_results(grid, "lmr")`.

## Interpretation

Models must use the same outcomes or indicators and observational units,
and their likelihoods must be comparable. Lower information criteria
favour a model relative to the alternatives considered, but do not
establish adequate fit or substantive validity. Classification
separation should inform interpretation rather than serve as the sole
basis for choosing the number of classes. Final selection also requires
examination of estimation stability, model assumptions, and the meaning
of the resulting classes.

## References

Wilks, S. S. (1938). The large-sample distribution of the likelihood
ratio for testing composite hypotheses. *The Annals of Mathematical
Statistics*, 9(1), 60–62.

Wagenmakers, E.-J., & Farrell, S. (2004). AIC model selection using
Akaike weights. *Psychonomic Bulletin & Review*, 11, 192–196.

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
