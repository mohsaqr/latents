# Reading growth of students nested in schools

A simulated multilevel longitudinal dataset made for the multilevel
growth mixture model. 600 students in 40 schools (15 per school) sit a
reading test on five occasions (`wave` 0 to 4). Two kinds of student
follow two trajectories: *improving* (starts near 48, gains 3 points a
wave) and *plateauing* (starts near 56, gains 0.3). Students vary around
their trajectory (random intercept SD 3, random slope SD 0.6,
correlation -0.2; residual SD 2). Schools are of two types: in
*supportive* schools 80% of students improve, in *struggling* schools
25% do, and a school's support programme (standardized, measured once)
makes it more likely to be supportive. About 8% of tests were missed at
random.

## Usage

``` r
growth_schools
```

## Format

A data frame with 2780 rows (600 students, 40 schools) and 7 columns:

- school:

  Integer school identifier, 1 to 40. Pass it as `cluster`.

- student:

  Integer student identifier, 1 to 600. Pass it as `id`.

- wave:

  Integer test occasion, 0 to 4. The time variable.

- score:

  Reading score.

- programme:

  Standardized school support programme, constant within school.

- school_type:

  Factor, `"supportive"` or `"struggling"`: the school's generating
  type. Not for fitting; kept to check recovery.

- trajectory:

  Factor, `"improving"` or `"plateauing"`: the student's generating
  kind. Not for fitting; kept to check recovery.

## Source

Simulated by `data-raw/growth-schools.R` in the source repository.

## Details

The school types are what a multilevel growth mixture model's group
classes recover:
`mixture_regression(..., cluster = "school", n_group_classes = 2)`, with
`group_membership = "programme"`.

## See also

[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
[`get_results.latents_growth_mixture()`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md).

## Examples

``` r
descriptives(growth_schools, c("score", "programme"))
#>    variable    n n_missing        mean        sd   min   max n_distinct
#> 1     score 2780         0 55.17003597 5.0187548 35.90 70.50        262
#> 2 programme 2780         0  0.04313669 0.9708406 -2.36  2.25         38
# \donttest{
fit <- mixture_regression(score ~ wave, growth_schools, n_classes = 2,
                          id = "student", class_level = "group",
                          random = "wave", random_covariance = "equal",
                          cluster = "school", n_group_classes = 2,
                          group_membership = "programme",
                          n_starts = 2, seed = 1)
get_results(fit, "group_classes")
#> Group classes: trajectory classes within each
#> 
#> Group class    Class    Probability     SE  Share of clusters  Clusters
#> -------------  -------  -----------  -----  -----------------  --------
#> Group class 1  Class 1        0.851  0.024              0.528        21
#> Group class 1  Class 2        0.149  0.024              0.528        21
#> Group class 2  Class 1        0.245  0.029              0.472        19
#> Group class 2  Class 2        0.755  0.029              0.472        19
get_results(fit, "recovery", data = growth_schools, truth = "trajectory")
#> Recovery of a known classification
#> 
#> Assigned  trajectory  Persons  Share of assigned
#> --------  ----------  -------  -----------------
#> Class 1   improving       321              0.950
#> Class 2   improving        14              0.053
#> Class 1   plateauing       17              0.050
#> Class 2   plateauing      248              0.947
# }
```
