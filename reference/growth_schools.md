# Reading growth of students nested in schools

`growth_schools` is a simulated dataset illustrating multilevel growth
mixture models. It contains repeated reading scores for 600 students
nested within 40 schools, with up to five observations per student.

The data were generated from two student trajectory classes, labelled
"improving" and "plateauing". Students vary in their initial reading
achievement and rate of change within each class. Two school-level
classes differ in the prevalence of these student trajectories.

The dataset illustrates how multilevel growth mixture models represent
individual variation in development alongside differences in trajectory
membership between schools. Generating class labels are retained at both
levels for assessing recovery.

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
  class; retained for assessing class recovery.

- trajectory:

  Factor, `"improving"` or `"plateauing"`: the student's generating
  class; retained for assessing class recovery.

## Source

Simulated by `data-raw/growth-schools.R` in the source repository.

## Details

There are 15 students per school, with occasions indexed by `wave` from
0 to 4. The generating student classes have two linear mean
trajectories: *improving* (intercept 48, slope 3 points per occasion)
and *plateauing* (intercept 56, slope 0.3). Random intercepts and slopes
have standard deviations 3 and 0.6 and correlation -0.2; the residual
standard deviation is 2. The generating school classes, labelled
*supportive* and *struggling*, have improving-class probabilities of 0.8
and 0.25, respectively. The standardized school-level `programme`
variable predicts school-class membership. Approximately 8% of test
occasions were removed at random.

The school types correspond to the group classes of
`mixture_regression(..., cluster = "school", n_group_classes = 2)`, with
`group_membership = "programme"`.

## See also

[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).

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
