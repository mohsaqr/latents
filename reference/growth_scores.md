# Achievement growth of students with three kinds of trajectory

A simulated longitudinal dataset made for growth models. 300 students
sit an achievement test on up to seven occasions (`wave` 0 to 6). Three
kinds of student follow three trajectories: *improving* (starts near 48,
gains 3.5 points a wave), *stable* (starts near 62, gains 0.5) and
*declining* (starts near 58, loses 2.5). Within each kind, students vary
around the trajectory: each has their own start and rate of change
(random intercept SD 4, random slope SD 0.8, correlation -0.3, the same
spread in every kind), with residual SD 3. Motivation, measured once,
makes a student more likely to be improving and less likely to be
declining. About 12% of tests were missed at random, so the panel is
unbalanced.

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
  generating kind. Not for fitting; kept to check recovery.

## Source

Simulated by `data-raw/growth-scores.R` in the source repository.

## Details

Because students vary within their kind, a trajectory model without
random effects needs more classes than the three that generated the
data; a growth mixture model with random intercepts and slopes recovers
three.

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
