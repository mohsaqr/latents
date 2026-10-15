# Compare fitted models on the same data

One row per model with its log likelihood, number of parameters,
information criteria, classification quality and convergence, so
alternative specifications (for example trajectory classes alone against
a growth mixture model) are compared without assembling tables by hand.
`delta_bic` is each model's BIC above the smallest, and `bic_weight` the
Schwarz weight, exp(-delta/2) normalized over the models compared
(Wagenmakers and Farrell 2004): the approximate posterior probability of
each model if one of them is true, with equal prior probabilities.

## Usage

``` r
compare_models(..., n_profiles = NULL, n_group_classes = NULL, cutoff = 0.8)
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
  (with or without `random`). Name fits to label the rows; unnamed fits
  are labelled by their expression, candidates by their name.

- n_profiles, n_group_classes:

  With a grid, the class counts of the candidates to compare, instead of
  naming them; `NULL` takes every count. Candidates whose fit failed are
  left out, with a message.

- cutoff:

  For profile and class models, the average posterior probability below
  which a class is counted in `avepp_below` and `group_avepp_below`. Not
  used for mixture regressions.

## Value

A base `data.frame` of class `latents_comparison` with one row per
model. For mixture regressions: `model`, `n_classes`, `random`,
`n_parameters`, `log_likelihood`, `aic`, `bic`, `sabic`, `icl`,
`entropy`, `smallest_share`, `delta_bic`, `bic_weight`, `converged`. For
profile and class models: `model`, `n_profiles`, `n_group_classes`,
`structure` (the covariance structure, `NA` without continuous
indicators), `n_parameters`, `log_likelihood`, every information
criterion of
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)'s
candidates table, `profile_entropy`, `group_entropy`, `avepp_min`,
`avepp_max`, `avepp_mean`, `avepp_below` and the same four prefixed
`group_`, `smallest_share` (the smallest estimated profile share,
relative to all observation mass including any noise component),
`converged`, `n_best_replicated`, and `delta_bic` and `bic_weight` on
`bic_groups` (the BIC that counts groups as the independent units; on a
single-level model it is the ordinary BIC). The numbers are the ones the
grid's candidates table reports. The printout marks the row with the
lowest BIC; that is a fact about the criterion, not a choice. Raises
`latents_bad_argument` for fits it cannot compare, for a group-class
family grid (compare those in their own table), for candidates named
together with `n_profiles` or `n_group_classes`, and for `n_profiles` or
`n_group_classes` without a grid; `latents_unknown_candidate` for a name
or counts not in the grid; `latents_failed_candidate` for a named
candidate that could not be fitted; and `latents_incomparable_models`
for fits to different data. A `cutoff` outside (0, 1\] fails the
argument checks.

## Details

Information criteria are comparable only on the same data: every model
must be fitted to the same outcome, persons (or groups) and
observations, or the comparison is refused.

Candidates of an enumeration grid are compared by naming them, as the
grid's `candidate` column shows them:
`compare_models(grid, "classes3_groups4", "classes3_groups5")`, or by
their class counts, `compare_models(grid, c(3, 4), c(3, 5))` (one number
is one group class). `n_profiles` and `n_group_classes` select
candidates instead: `compare_models(grid, n_profiles = 3)` compares
every candidate with 3 profiles; `compare_models(grid)` compares them
all. Latent class models are compared the same way; their `n_profiles`
column counts classes. Fits from
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
[`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
can also be passed directly.

## References

Wagenmakers, E.-J., & Farrell, S. (2004). AIC model selection using
Akaike weights. *Psychonomic Bulletin & Review*, 11, 192–196.

## Examples

``` r
# Fixed trajectories against growth curves with a random intercept,
# for the first 60 students
few_students <- subset(growth_scores, student <= 60)
trajectories <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                                   id = "student", class_level = "group",
                                   seed = 1)
growth <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                             id = "student", class_level = "group",
                             random = "intercept", seed = 1)
compare_models(trajectories = trajectories, growth = growth)
#> Model comparison
#> 
#> Model         Classes  Random     k    LogLik      BIC    ΔBIC  Weight  Entropy
#> ------------  -------  ---------  -  --------  -------  ------  ------  -------
#> trajectories        2  none       7  -1277.94  2584.53  196.11   0.000     0.94
#> growth              2  intercept  9  -1175.79  2388.43    0.00   1.000     0.97  <- lowest BIC

# Candidates of an enumeration grid, by name, by class counts, or selected
grid <- enumerate_classes(subset(course_engagement, student <= 40),
                          c("browse", "lectures"), "student",
                          n_profiles = 2:3, n_group_classes = 1:2,
                          model = "VVI", n_starts = 2, seed = 1)
compare_models(grid, "profiles2_groups2", "profiles3_groups2")
#> Model comparison
#> 
#> Model              Profiles  Group classes   k    LogLik      BIC   ΔBIC  Weight      ICL  Entropy  Min AvePP  < 0.8  Group entropy  Group min AvePP
#> -----------------  --------  -------------  --  --------  -------  -----  ------  -------  -------  ---------  -----  -------------  ---------------
#> profiles2_groups2         2              2  11  -1458.86  2958.29   0.00   0.999  3233.59     0.67       0.90      0           0.96             0.98  <- lowest BIC
#> profiles3_groups2         3              2  17  -1454.49  2971.70  13.41   0.001  3409.04     0.67       0.72      1           0.96             0.97
compare_models(grid, c(2, 2), c(3, 2))
#> Model comparison
#> 
#> Model              Profiles  Group classes   k    LogLik      BIC   ΔBIC  Weight      ICL  Entropy  Min AvePP  < 0.8  Group entropy  Group min AvePP
#> -----------------  --------  -------------  --  --------  -------  -----  ------  -------  -------  ---------  -----  -------------  ---------------
#> profiles2_groups2         2              2  11  -1458.86  2958.29   0.00   0.999  3233.59     0.67       0.90      0           0.96             0.98  <- lowest BIC
#> profiles3_groups2         3              2  17  -1454.49  2971.70  13.41   0.001  3409.04     0.67       0.72      1           0.96             0.97
compare_models(grid, n_group_classes = 2)
#> Model comparison
#> 
#> Model              Profiles  Group classes   k    LogLik      BIC   ΔBIC  Weight      ICL  Entropy  Min AvePP  < 0.8  Group entropy  Group min AvePP
#> -----------------  --------  -------------  --  --------  -------  -----  ------  -------  -------  ---------  -----  -------------  ---------------
#> profiles2_groups2         2              2  11  -1458.86  2958.29   0.00   0.999  3233.59     0.67       0.90      0           0.96             0.98  <- lowest BIC
#> profiles3_groups2         3              2  17  -1454.49  2971.70  13.41   0.001  3409.04     0.67       0.72      1           0.96             0.97
```
