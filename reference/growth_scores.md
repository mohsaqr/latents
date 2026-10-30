# Simulated achievement scores with three growth trajectories

`growth_scores` is a simulated longitudinal dataset illustrating growth
mixture models. Each row represents an achievement score for one of 300
students, measured on up to seven occasions.

The data were generated from three trajectory classes, labelled
"improving", "stable" and "declining". Within each class, students
differ in their initial achievement and rate of change through random
intercepts and slopes.

The dataset illustrates how growth mixture models distinguish
differences between class mean trajectories from variation among
individuals within the same class. Generating class labels are retained
for assessing recovery, and randomly omitted observations produce an
unbalanced panel.

## Usage

``` r
growth_scores
```

## Format

A data frame with 1844 rows (300 students) and 5 columns:

- student:

  Integer student identifier, 1 to 300. Pass it as `id`.

- wave:

  Integer test occasion, 0 to 6. The time variable.

- score:

  Achievement score.

- motivation:

  Standardized motivation, constant within student.

- trajectory:

  Factor, `"improving"`, `"stable"` or `"declining"`: the student's
  generating class; retained for assessing class recovery.

## Source

Simulated by `data-raw/growth-scores.R` in the source repository.

## Details

Occasions are indexed by `wave`, from 0 to 6. The generating classes
have three linear mean trajectories: *improving* (intercept 48, slope
3.5 points per occasion), *stable* (intercept 62, slope 0.5), and
*declining* (intercept 58, slope -2.5). Within each class, individual
trajectories vary according to random intercepts and slopes with
standard deviations 4 and 0.8 and correlation -0.3. The residual
standard deviation is 3. Motivation is a time-invariant predictor of
class membership, with higher values associated with greater probability
of the improving class and lower probability of the declining class.
Approximately 12% of test occasions were removed at random, producing an
unbalanced panel.

A model without random effects may represent within-class variation by
additional trajectory classes.

## See also

[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).

## Examples

``` r
descriptives(growth_scores, c("score", "motivation"))
#>     variable    n n_missing        mean       sd   min   max n_distinct
#> 1      score 1844         0 57.57039046 8.909898 24.90 81.20        400
#> 2 motivation 1844         0  0.03240239 1.006345 -3.35  2.82        200
# \donttest{
fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 3,
                          id = "student", class_level = "group",
                          random = "wave", random_covariance = "equal",
                          n_starts = 4, seed = 1)
get_results(fit, "recovery", data = growth_scores, truth = "trajectory")
#> Recovery of a known classification
#> 
#> Assigned  trajectory  Persons  Share of assigned
#> --------  ----------  -------  -----------------
#> Class 1   improving       109              0.973
#> Class 2   improving         0              0.000
#> Class 3   improving         6              0.067
#> Class 1   stable            3              0.027
#> Class 2   stable            3              0.030
#> Class 3   stable           82              0.921
#> Class 1   declining         0              0.000
#> Class 2   declining        96              0.970
#> Class 3   declining         1              0.011
# }
```
